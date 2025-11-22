##############################################################################
## Elastic net regularization to identify the gene signature for mAge acceleration
##############################################################################

## load R packages

library(tidyverse) 
library(caret)     
library(glmnet)    
library(pROC)      
library(themis)
library(glmnet)
library(GSEABase)
library(AUCell)


##############################################################################
## dataset pre-processing
##############################################################################

set.seed(123) 

## initial gene list, commonly up or down regulated genes in cohorts of Ruijin and Schmitz et al.
gset_acc <- getGmt("2025.08.Aging.Datasets.ACC.Dec.DEGs.gmt")
gset_acc <- subsetGeneSets(gset_acc, rownames(expression_data_nor)) 

## extract gene names
up_genes = geneIds(gset_acc)$RJ.NEJM18.Acc.Up
down_genes = geneIds(gset_acc)$RJ.NEJM18.Acc.Dn
candidate_genes <- union(up_genes, down_genes)

## metadata with mAge groups results
pheno = age_results_sub
pheno$group <- factor(pheno$mAge_group, levels = c("Dec","Syn","Acc"))

## extract expression data of the initial gene list
X_train <- expression_data_nor[candidate_genes, , drop=FALSE]

## data scaling
train_means <- apply(X_train, 2, mean)
train_sds <- apply(X_train, 2, sd)
X_train_scaled <- scale(X_train, center = train_means, scale = train_sds)

## align the samples in the metadata and expression matrix
samps <- if ("SampleID" %in% names(pheno)) {
  intersect(pheno$SampleID, colnames(X_train))
} else {
  intersect(rownames(pheno), colnames(X_train))
}
stopifnot(length(samps) > 100)  # set a minimal sample size

## reorder samples
if ("SampleID" %in% names(pheno)) {
  ph <- pheno[match(samps, pheno$SampleID), , drop=FALSE]
  rownames(ph) <- ph$SampleID
} else {
  ph <- pheno[match(samps, rownames(pheno)), , drop=FALSE]
}
X_train <- t(X_train[, samps, drop=FALSE])


##############################################################################
## Training for discovery model
##############################################################################


##### expression matrix #####
X_train<- as.data.frame(X_train) %>% rownames_to_column("SampleID")

## add mAge group information into the expression matrix
internal_model_data <- pheno[,c("SampleID","group")] %>% left_join(X_train, by = "SampleID")
rownames(internal_model_data) = internal_model_data$SampleID
internal_model_data <- subset(internal_model_data, select = -SampleID)

internal_model_data$group <- factor(internal_model_data$group, levels = c("Dec", "Syn", "Acc"))

cat("mAge group classification:\n")
print(table(internal_model_data$group))

##### assign weight to mAge groups #####
## The weights are inversely proportional to the class frequencies.

class_counts <- table(internal_model_data$group)
model_weights <- sum(class_counts) / (length(class_counts) * class_counts)  # normalize weight, so mean equals to 1
model_weights <- model_weights[internal_model_data$group]                   
mean(model_weights)

##### Train an Elastic Net model and perform gene selection #####

ctrl <- trainControl(method = "cv", 
                     number = 10,
                     classProbs = TRUE,
                     summaryFunction = multiClassSummary,
                     savePredictions = "final")

## optional: parameter grid for tuning
tune_grid <- expand.grid(alpha = seq(0.1, 1, length = 5),
                         lambda = 10^seq(-4, -2, length = 20))

##
cat("\nStart discovery model training...\n")
set.seed(123)
discovery_model <- train(
  group ~ .,
  data = internal_model_data,
  family = "multinomial",
  method = "glmnet",
  trControl = ctrl,
  tuneLength = 10,      # tuneGrid = tune_grid
  metric = "Accuracy",  
  preProcess = c("center", "scale"),
  weights = model_weights 
)
cat("Training completed\n")



##############################################################################
## Signature identification with discovery cohort
##############################################################################

## Extract and rank all potential genes (by importance/coefficient)
initial_coeffs <- coef(discovery_model$finalModel, s = discovery_model$bestTune$lambda)
accel_coeffs <- as.data.frame(as.matrix(initial_coeffs$Acc))
colnames(accel_coeffs) = "Coefficient"

ranked_genes_df <- accel_coeffs %>%
  rownames_to_column("Gene") %>%
  filter(Gene != "(Intercept)" & Coefficient != 0) %>%
  arrange(desc(abs(Coefficient)))

ranked_gene_list <- ranked_genes_df$Gene

## define the range for number of genes included in the training
## We tested gene numbers in the range of 20 to 60 at intervals of 5
gene_range_to_test <- seq(20, 60, by = 5) 
cat("Testing gene number:", paste(gene_range_to_test, collapse=", "), "\n")

## Create a Data Frame for storing iterative results
performance_summary <- data.frame()

## Commence the iterative testing.
cat("Start the iterative assessment with different numbers of genes...\n")
for (n_genes in gene_range_to_test) {
  
  # next if not enough candidate genes
  if (n_genes > length(ranked_gene_list)) next
  
  # select top n genes by their coefficients
  current_signature <- ranked_gene_list[1:n_genes]
  
  # subset with the top n genes
  internal_data_subset <- internal_model_data %>%
    dplyr::select(group, all_of(current_signature))
  
  # Evaluate the performance using cross-validation.
  set.seed(456) # use a different seed than the main function
  cv_fit <- train(
    group ~ ., data = internal_data_subset, method = "glmnet",
    trControl = ctrl, 
    tuneLength = 10,      # tuneGrid = tune_grid
    metric = "Accuracy", 
    preProcess = c("center", "scale"),
    weights = model_weights
  )
  
  performance_summary <- rbind(performance_summary, 
                               data.frame(Num_Genes = n_genes,
                                          Accuracy = cv_fit$results[which.max(cv_fit$results$Accuracy), "Accuracy"],
                                          Kappa = cv_fit$results[which.max(cv_fit$results$Accuracy), "Kappa"]))
  
  cat(sprintf("Completed: %d genes, CV Accuracy = %.4f\n", n_genes, performance_summary$Accuracy[nrow(performance_summary)]))
}

# visualize the result with elbow_plot
elbow_plot <- ggplot(performance_summary, aes(x = Num_Genes, y = Accuracy)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 3) +
  labs(title = "Gene numbers vs. Accuracy",
       x = "Number of Genes",
       y = "CV (Accuracy)") +
  scale_x_continuous(breaks = gene_range_to_test) +
  theme_bw(base_size = 14)

print(elbow_plot)

## Select the optimal number of genes based on the point with the highest accuracy.
best_performance <- performance_summary %>% filter(Accuracy == max(Accuracy)) %>% slice(1) 
best_n_genes <- best_performance$Num_Genes

final_signature_small <- ranked_gene_list[1:best_n_genes]



##############################################################################
## Validation of the final gene signature on discovery cohort
##############################################################################

###########################################
##### discovery cohort: Ruijin cohort #####

## extract expression data of the signature genes
internal_data_small <- internal_model_data %>%
  dplyr::select(group, all_of(final_signature_small))



##### Perform cross-validation of the signature on the discovery cohort #####
cat("Cross-validation of the signature on the discovery cohort...\n")
set.seed(123)
internal_model_final <- train(
              group ~ .,
              data = internal_data_small,
              method = "glmnet",
              trControl = ctrl,
              tuneLength = 10,      # tuneGrid = tune_grid
              metric = "Accuracy",  
              family = "multinomial",
              preProcess = c("center", "scale"),
              weights = model_weights
)

## mAge group prediction
pred_RJ <- predict(internal_model_final, newdata = internal_data_small)
prob_RJ <- predict(internal_model_final, newdata = internal_data_small, type = "prob")

## confusion matrix
confusionMatrix(pred_RJ, internal_data_small$group)



##### ROC curve #####
roc_list <- list()
for (class in levels(internal_data_small$group)) {
  binary_y <- as.factor(ifelse(internal_data_small$group == class, "Yes", "No"))
  # ROC curve
  roc_obj <- roc(response = binary_y, predictor = prob_RJ[, class])
  roc_list[[class]] <- roc_obj
}
## AUC of the ROC curve
roc_list

# Plot ROC curve
pdf("./datasets/Ruijin/RJ.Acc_Signature.ROC.pdf",height=4,width=4)
plot(roc_list$Acc, col = "#BF360C", main = "Ruijin cohort",lwd = 3)
plot(roc_list$Syn, col = "#0091EA", add = TRUE,lwd = 3)
plot(roc_list$Dec, col = "#4CAF50", add = TRUE,lwd = 3)
legend("bottomright", legend = c(paste("Dec: AUC",round(roc_list$Dec$auc,digits = 3)),
                                 paste("Syn: AUC",round(roc_list$Syn$auc,digits = 3)),
                                 paste("Acc: AUC",round(roc_list$Acc$auc,digits = 3))),
       col = c("#4CAF50", "#0091EA", "#BF360C"), lwd = 5, seg.len = 0.5)
dev.off()



##### plot confusion matrix #####

## function to plot the matrix
ggplot_confusion_matrix <- function(cfm) {
    mytitle <- paste("Accuracy", percent_format() (cfm$overall[1]),
                     "Kappa", percent_format() (cfm$overall[2]))
    p <-
        ggplot(data = as.data.frame(cfm$table),
               aes(x = Reference, y = Prediction)) +
        geom_tile(aes(fill = log(Freq)), colour = "white") +
        scale_fill_gradient(low = "white", high = "steelblue") +
        geom_text(aes(x = Reference, y = Prediction, label = Freq)) +
        theme_pubr(base_size = 12) +
        ggtitle(mytitle)
    return(p)
}

internal_data_small$group <- factor(internal_data_small$group, levels = c("Dec", "Syn", "Acc"))
levels(pred_RJ) = c("Dec", "Syn", "Acc")

p = ggplot_confusion_matrix(confusionMatrix(pred_RJ, internal_data_small$group))
ggsave("./datasets/Ruijin/2025.08.29.confusionMtx.RJ.pdf",p,width = 2.2,height = 3)


#####################################################
##### discovery cohort: Schmitz et al NEJM 2018 #####
## Expression matrix downloaded from: https://gdc.cancer.gov/about-data/publications/DLBCL-2018

## reading the normalized matrix
exprs_2018NEJM = read.table("./datasets/Schmitz_2018NEJM/RNAseq_gene_expression_562.formated.txt",header=T,sep = "\t")
rownames(exprs_2018NEJM) <- exprs_2018NEJM$Gene
exprs_2018NEJM <- exprs_2018NEJM[,c(-1,-2,-3)]

## metadata with mAge group information
age_info_2018NEJM = read.table("./datasets/Schmitz_2018NEJM/2018NEJM.age_results_sub.txt",header=T,sep = "\t")

## convertting the FPKM to TPM values
expr_fpkm <- 2^exprs_2018NEJM - 1

fpkm_to_tpm <- function(fpkm) {
  scaling <- colSums(fpkm, na.rm = TRUE) 
  tpm <- sweep(fpkm, 2, scaling, "/") * 1e6
  return(tpm)
}

expr_tpm <- fpkm_to_tpm(expr_fpkm)
expr_log2_tpm <- log2(expr_tpm + 1)

rm(expr_fpkm,exprs_2018NEJM)


## Ensure sample names match
common_samples <- intersect(colnames(expr_log2_tpm), age_info_2018NEJM$SampleID)
expr_log2_tpm <- expr_log2_tpm[,common_samples]
expr_log2_tpm = t(expr_log2_tpm)

## align the sample order
age_info_2018NEJM <- age_info_2018NEJM %>% filter(SampleID %in% common_samples)
age_info_2018NEJM <- age_info_2018NEJM[match(colnames(expr_log2_tpm), age_info_2018NEJM$SampleID), ]
age_info_2018NEJM$group = mAge_group

## extract expression matrix
X_test<- as.data.frame(expr_log2_tpm) %>% rownames_to_column("SampleID")

## add mAge group information into the expression matrix
external_2018NEJM <- age_info_2018NEJM[,c("SampleID","group")] %>% left_join(X_test, by = "SampleID")
rownames(external_2018NEJM) = external_2018NEJM$SampleID
external_2018NEJM <- subset(external_2018NEJM, select = -SampleID)

external_2018NEJM$group <- factor(external_2018NEJM$group, 
                                    levels = c("Dec", "Syn", "Acc"))

## assign weight to mAge groups
## ## The weights are inversely proportional to the class frequencies.

ext_counts <- table(external_2018NEJM$group)
ext_weights <- sum(ext_counts) / (length(ext_counts) * ext_counts)  # normalize weight, so mean equals to 1
ext_weights <- ext_weights[external_2018NEJM$group]                  
mean(ext_weights)

## extract expression for the final signature genes
external_2018NEJM <- external_2018NEJM %>%
  dplyr::select(group, all_of(final_signature_small))

## Cross-validation of the signature
cat("Cross-validation of the signature...\n")
ctrl_external <- trainControl(method = "cv", number = 10, classProbs = TRUE, 
                              summaryFunction = multiClassSummary, savePredictions = "final")

set.seed(123)
external_model_2018NEJM <- train(
          group ~ .,
          data = external_2018NEJM,
          method = "glmnet",
          family = "multinomial",
          preProcess = c("center", "scale"),
          trControl = ctrl_external,
          tuneLength = 10,      # tuneGrid = tune_grid
          weights = ext_weights, 
          metric = "Accuracy"
)


## mAge group prediction
pred_2018NEJM <- predict(external_model_2018NEJM, newdata = external_2018NEJM)
prob_2018NEJM <- predict(external_model_2018NEJM, newdata = external_2018NEJM, type = "prob")

## confusion Matrix
confusionMatrix(pred_2018NEJM, external_2018NEJM$group)

roc_list <- list()
for (class in levels(external_2018NEJM$group)) {
  binary_y <- as.factor(ifelse(external_2018NEJM$group == class, "Yes", "No"))
  # ROC curve
  roc_obj <- roc(response = binary_y, predictor = prob_2018NEJM[, class])
  roc_list[[class]] <- roc_obj
}

## Plot ROC curve
roc_list

pdf("./datasets/Schmitz_2018NEJM/2018NEJM.Acc_Signature.ROC.pdf",height=4,width=4)
plot(roc_list$Acc, col = "#BF360C", main = "Schmitz et al",lwd = 3)
plot(roc_list$Syn, col = "#0091EA", add = TRUE,lwd = 3)
plot(roc_list$Dec, col = "#4CAF50", add = TRUE,lwd = 3)
legend("bottomright", legend = c(paste("Dec: AUC",round(roc_list$Dec$auc,digits = 3)),
                                 paste("Syn: AUC",round(roc_list$Syn$auc,digits = 3)),
                                 paste("Acc: AUC",round(roc_list$Acc$auc,digits = 3))),
       col = c("#4CAF50", "#0091EA", "#BF360C"), lwd = 5, seg.len = 0.5)
dev.off()


##### plot confusion matrix #####
external_2018NEJM$group <- factor(external_2018NEJM$group, levels = c("Dec", "Syn", "Acc"))

p = ggplot_confusion_matrix(confusionMatrix(pred_2018NEJM, external_2018NEJM$group))
ggsave("./datasets/Schmitz_2018NEJM/2025.08.29.confusionMtx.2018NEJM.pdf",p,width = 2.2,height = 3)


##############################################################################
## Validation of the final gene signature on additional cohorts
##############################################################################

##################################################
##### DLBCL dataset: Sha et al JCO 2018 #####
## Expression matrix downloaded from: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117556/matrix/

## reading the normalized matrix
exprs_2018JCO = read.table("./datasets/Sha_2018JCO/GSE117556_series_matrix.nor.symbol.txt",header=T,sep = "\t")

## metadata with mAge group information
age_info_2018JCO = read.table("./datasets/Sha_2018JCO/2025.08.25.2018JCO.age_results_sub.txt",header=T,sep = "\t")

## Ensure sample names match
common_samples <- intersect(colnames(exprs_2018JCO), age_info_2018JCO$SampleID)
exprs_2018JCO <- exprs_2018JCO[,common_samples]
age_info_2018JCO <- age_info_2018JCO %>% filter(SampleID %in% common_samples)

## align the sample order
age_info_2018JCO <- age_info_2018JCO[match(colnames(exprs_2018JCO), age_info_2018JCO$SampleID), ]

## extract expression matrix
X_test<- as.data.frame(t(exprs_2018JCO)) %>% rownames_to_column("SampleID")
age_info_2018JCO$group = age_info_2018JCO$mAge_group

## add mAge group information into the expression matrix
external_2018JCO <- age_info_2018JCO[,c("SampleID","group")] %>% left_join(X_test, by = "SampleID")
rownames(external_2018JCO) = external_2018JCO$SampleID
external_2018JCO <- subset(external_2018JCO, select = -SampleID)

external_2018JCO$group <- factor(external_2018JCO$group, 
                                    levels = c("Dec", "Syn", "Acc"))

## assign weight to mAge groups
## The weights are inversely proportional to the class frequencies
ext_counts <- table(external_2018JCO$group)

ext_weights <- sum(ext_counts) / (length(ext_counts) * ext_counts)  # normalize weight, so mean equals to 1
ext_weights <- ext_weights[external_2018JCO$group]                   
mean(ext_weights)

## extract expression for the final signature genes
external_2018JCO <- external_2018JCO %>%
  dplyr::select(group, all_of(intersect(final_signature_small, colnames(external_2018JCO))))


## Cross-validation of the signature
cat("Cross-validation of the signature...\n")
ctrl_external <- trainControl(method = "cv", number = 10, classProbs = TRUE, 
                              summaryFunction = multiClassSummary, savePredictions = "final")

set.seed(123)
external_model_2018JCO <- train(
          group ~ .,
          data = external_2018JCO,
          method = "glmnet",
          family = "multinomial",
          # preProcess = c("center", "scale"),
          trControl = ctrl_external,
          tuneLength = 10,      # tuneGrid = tune_grid
          weights = ext_weights, 
          metric = "Accuracy"  # or "logLoss"
)


## mAge group prediction
pred_2018JCO <- predict(external_model_2018JCO, newdata = external_2018JCO)
prob_2018JCO <- predict(external_model_2018JCO, newdata = external_2018JCO, type = "prob")

## confusion Matrix
confusionMatrix(pred_2018JCO, external_2018JCO$group)

roc_list <- list()
for (class in levels(external_2018JCO$group)) {
  binary_y <- as.factor(ifelse(external_2018JCO$group == class, "Yes", "No"))
  roc_obj <- roc(response = binary_y, predictor = prob_2018JCO[, class])
  roc_list[[class]] <- roc_obj
}

## AUC of the ROC curve
roc_list

## Plot ROC curve
pdf("./datasets/Sha_2018JCO/2018JCO.Acc_Signature.ROC.pdf",height=4,width=4)
plot(roc_list$Acc, col = "#BF360C", main = "Sha et al",lwd = 3)
plot(roc_list$Syn, col = "#0091EA", add = TRUE,lwd = 3)
plot(roc_list$Dec, col = "#4CAF50", add = TRUE,lwd = 3)
legend("bottomright", legend = c(paste("Dec: AUC",round(roc_list$Dec$auc,digits = 3)),
                                 paste("Syn: AUC",round(roc_list$Syn$auc,digits = 3)),
                                 paste("Acc: AUC",round(roc_list$Acc$auc,digits = 3))),
       col = c("#4CAF50", "#0091EA", "#BF360C"), lwd = 5, seg.len = 0.5)
dev.off()

##### plot confusion matrix #####
external_2018JCO$group <- factor(external_2018JCO$group, levels = c("Dec", "Syn", "Acc"))
# levels(pred_RJ) = c("Dec", "Syn", "Acc")

p = ggplot_confusion_matrix(confusionMatrix(pred_2018JCO, external_2018JCO$group))
ggsave("./datasets/Sha_2018JCO/2025.08.29.confusionMtx.2018JCO.pdf",p,width = 2.2,height = 3)

##################################################
##### DLBCL dataset: Lacy et al Blood 2020 #####
## Expression matrix downloaded from: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE181nnn/GSE181063/matrix/

## reading the normalized matrix
exprs_2020Blood = read.table("./datasets/Lacy_2020Blood/GSE181063_series_matrix.nor.DLBCL.txt",header=T,sep = "\t")

## metadata with mAge group information
age_info_2020Blood = read.table("./datasets/Lacy_2020Blood/2025.07.25.2020Blood.age_results_sub.txt",header=T,sep = "\t")

## Ensure sample names match
common_samples <- intersect(colnames(exprs_2020Blood), age_info_2020Blood$SampleID)
exprs_2020Blood <- exprs_2020Blood[,common_samples]
age_info_2020Blood <- age_info_2020Blood %>% filter(SampleID %in% common_samples)

## align the sample order
age_info_2020Blood <- age_info_2020Blood[match(colnames(exprs_2020Blood), age_info_2020Blood$SampleID), ]

## extract expression matrix
X_test<- as.data.frame(t(exprs_2020Blood)) %>% rownames_to_column("SampleID")
age_info_2020Blood$group = age_info_2020Blood$mAge_group

## add mAge group information into the expression matrix
external_2020Blood <- age_info_2020Blood[,c("SampleID","group")] %>% left_join(X_test, by = "SampleID")
rownames(external_2020Blood) = external_2020Blood$SampleID
external_2020Blood <- subset(external_2020Blood, select = -SampleID)

external_2020Blood$group <- factor(external_2020Blood$group, 
                                    levels = c("Dec", "Syn", "Acc"))

## assign weight to mAge groups
## The weights are inversely proportional to the class frequencies.
ext_counts <- table(external_2020Blood$group)
ext_weights <- sum(ext_counts) / (length(ext_counts) * ext_counts)  # normalize weight, so mean equals to 1
ext_weights <- ext_weights[external_2020Blood$group]                   
mean(ext_weights)

## extract expression for the final signature genes
external_2020Blood <- external_2020Blood %>%
  dplyr::select(group, all_of(intersect(final_signature_small, colnames(external_2020Blood))))


## Cross-validation of the signature
cat("Cross-validation of the signature...\n")
ctrl_external <- trainControl(method = "cv", number = 10, classProbs = TRUE, 
                              summaryFunction = multiClassSummary, savePredictions = "final")

set.seed(123)
external_model_validation <- train(
          group ~ .,
          data = external_2020Blood,
          method = "glmnet",
          family = "multinomial",
          preProcess = c("center", "scale"),
          trControl = ctrl_external,
          tuneLength = 10,      # tuneGrid = tune_grid
          weights = ext_weights, 
          metric = "Accuracy"  # or "logLoss"
)

## mAge group prediction
pred_2020Blood <- predict(external_model_validation, newdata = external_2020Blood)
prob_2020Blood <- predict(external_model_validation, newdata = external_2020Blood, type = "prob")

## confusion Matrix
confusionMatrix(pred_2020Blood, external_2020Blood$group)

## AUC of the ROC curve
roc_list <- list()
for (class in levels(external_2020Blood$group)) {
  binary_y <- as.factor(ifelse(external_2020Blood$group == class, "Yes", "No"))
  roc_obj <- roc(response = binary_y, predictor = prob_2020Blood[, class])
  roc_list[[class]] <- roc_obj
}

roc_list

## Plot ROC curve
pdf("./datasets/Lacy_2020Blood/2020Blood.Acc_Signature.ROC.pdf",height=4,width=4)
plot(roc_list$Acc, col = "#BF360C", main = "Lacy et al",lwd = 3)
plot(roc_list$Syn, col = "#0091EA", add = TRUE,lwd = 3)
plot(roc_list$Dec, col = "#4CAF50", add = TRUE,lwd = 3)
legend("bottomright", legend = c(paste("Dec: AUC",round(roc_list$Dec$auc,digits = 3)),
                                 paste("Syn: AUC",round(roc_list$Syn$auc,digits = 3)),
                                 paste("Acc: AUC",round(roc_list$Acc$auc,digits = 3))),
       col = c("#4CAF50", "#0091EA", "#BF360C"), lwd = 5, seg.len = 0.5)
dev.off()

##### plot confusion matrix #####
external_2020Blood$group <- factor(external_2020Blood$group, levels = c("Dec", "Syn", "Acc"))

p = ggplot_confusion_matrix(confusionMatrix(pred_2020Blood, external_2020Blood$group))
ggsave("./datasets/Lacy_2020Blood/2025.08.29.confusionMtx.2020Blood.pdf",p,width = 2.2,height = 3)

