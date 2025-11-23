###############################################################
## load R package
###############################################################

setwd("./")

library(dplyr)
library(ggplot2)
library(limma)
library(edgeR)
library(scales)
library(ggrepel)
library(tidyr)

library(survival)
library(survminer)
library(paletteer)

library(viridis)
library(ggpointdensity)
library(RColorBrewer)
library(corrplot)
library(GSVA)
library(GSEABase)
library(AUCell)
library(edgeR)     
library(limma)      
library(splines)   
library(sva)


###############################################################
## expression matrix and patient data
###############################################################


###############################################################
##### Ruijin cohort #####

## expression matrix
expression_data <- read.table("DLBCL_759_RNAseq_rawcounts.txt",header=T)
rownames(expression_data) <- expression_data$gene
expression_data <- expression_data[,-1]

## match the expression matrix and patient data
group_data <- read.table("DLBCL_759.patient.info.txt",header=T,fill = TRUE,sep = "\t")

age_info <- group_data[match(colnames(expression_data), group_data$Seq_ID), ]

## correct batch effect with ComBat_seq
expression_data_adj <- ComBat_seq(as.matrix(expression_data), batch=group_data$batch, group=NULL)

## expression normalization
dge <- DGEList(counts = expression_data_adj)
dge <- calcNormFactors(dge)  
design0 <- model.matrix(~1, data = age_info)       
v <- voom(dge, design = design0)   

expression_data_nor = v$E



###############################################################
##### DLBCL cohort: Schmitz et al NEJM 2018 #####
## Expression matrix downloaded from: https://gdc.cancer.gov/about-data/publications/DLBCL-2018

## expression matrix
expression_data = read.table("./datasets/Schmitz_2018NEJM/RNAseq_gene_expression_562.formated.txt",header=T,sep = "\t")
rownames(expression_data) <- expression_data$Gene
expression_data <- expression_data[,c(-1,-2,-3)]

## patient metadata
group_data <- read.table("Supplementary_Appendix_2.TableS9_patients.txt",header=T,fill = TRUE,sep = "\t")
group_data$sampleID = gsub("-", ".", group_data$dbGaP.submitted.subject.ID)
age_info = subset(group_data, Included.in.Survival.Analysis == "Yes")


## Initial expression matrix was FPKM values, convert the FPKM to TPM values
expr_fpkm <- 2^expression_data_ag - 1

fpkm_to_tpm <- function(fpkm) {
  scaling <- colSums(fpkm, na.rm = TRUE)   # 每列 FPKM 总和
  tpm <- sweep(fpkm, 2, scaling, "/") * 1e6
  return(tpm)
}

expr_tpm <- fpkm_to_tpm(expr_fpkm)
expression_data_nor <- log2(expr_tpm + 1)

## Ensure sample names match
common_samples <- intersect(colnames(expression_data_nor), age_info$sampleID)
tpm_t <- expression_data_nor[,common_samples]
age_info <- age_info %>% filter(sampleID %in% common_samples)

## order the metadata by the expression matrix
age_info <- age_info[match(colnames(expression_data_nor), age_info$sampleID), ]



###############################################################
##### DLBCL cohort: Sha et al JCO 2018 #####
## Expression matrix downloaded from: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117556/matrix/

## expression matrix
expression_data = read.table("GSE117556_series_matrix.nor.txt",header=T,sep = "\t")
rownames(expression_data) <- expression_data$ID_REF
expression_data <- expression_data[,-1]

## match the expression matrix and patient data
group_data <- read.table("GSE117556.suppl2.patient-info.txt",header=T,fill = TRUE,sep = "\t")
age_info = subset(group_data, AGE != "n/a")

## gene annotation
library("AnnotationDbi")
library("illuminaHumanv4.db")
library("org.Hs.eg.db")                                          
columns(illuminaHumanv4.db)
library(dplyr)
library(tibble)

id_str <- rownames(expression_data)
annotable = AnnotationDbi::select(illuminaHumanv4.db,id_str,c("SYMBOL"), keytype="PROBEID")

expression_data_df <- expression_data %>%  tibble::rownames_to_column(var = "PROBEID")
expression_data_df <- expression_data_df %>% left_join(annotable, by = "PROBEID")

expression_data_nor <- expression_data_df %>% filter(!is.na(SYMBOL))

## collapse multiple probes with mean
expression_data_nor <- expression_data_nor %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), 
   ~ if(n() > 1) mean(.x, na.rm = TRUE) else .x[1])) %>% ungroup()

expression_data_nor <- expression_data_nor %>% column_to_rownames(var = "SYMBOL")

## Ensure sample names match
common_samples <- intersect(colnames(expression_data_nor), age_info$geo)
expression_data_nor <- expression_data_nor[,common_samples]
age_info <- age_info %>% filter(geo %in% common_samples)

## order the metadata by the expression matrix
age_info <- age_info[match(colnames(expression_data_nor), age_info$geo), ]


###############################################################
##### DLBCL cohort: Chapuy et al NatMed 2018 #####
## Expression matrix downloaded from: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98588/matrix/

## expression matrix
expression_data = read.table("./datasets/Chapuy_2018NatMed/GSE98588_series_matrix.nor.txt",header=T,sep = "\t")
rownames(expression_data) <- expression_data$ID_REF
expression_data <- expression_data[,-1]

## match the expression matrix and patient data
group_data <- read.table("GSE98588_series_matrix.sample.info.txt",header=T,fill = TRUE,sep = "\t")
age_info = subset(group_data, PFS != "na")

## gene annotation
library("AnnotationDbi")
library("illuminaHumanv4.db")
library("org.Hs.eg.db")                                          
columns(org.Hs.eg.db)

library(dplyr)
library(tibble)

id_str <- rownames(expression_data)
annotable = AnnotationDbi::select(org.Hs.eg.db,id_str,c("SYMBOL"), keytype="ENSEMBL")

expression_data_df <- expression_data %>%  tibble::rownames_to_column(var = "ENSEMBL")
expression_data_df <- expression_data_df %>% left_join(annotable, by = "ENSEMBL")

expression_data_nor <- expression_data_df %>% filter(!is.na(SYMBOL))

## collapse multiple probes with mean
expression_data_nor <- expression_data_nor %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), 
   ~ if(n() > 1) mean(.x, na.rm = TRUE) else .x[1])) %>% ungroup()

expression_data_nor <- expression_data_nor %>% column_to_rownames(var = "SYMBOL")

## Ensure sample names match
common_samples <- intersect(colnames(expression_data_nor), age_info$GSE_acc)
expression_data_nor <- expression_data_nor[,common_samples]
age_info <- age_info %>% filter(GSE_acc %in% common_samples)

## order the metadata by the expression matrix
age_info <- age_info[match(colnames(expression_data_nor), age_info$GSE_acc), ]



###############################################################
##### DLBCL cohort: Lacy et al Blood 2020 #####
## Expression matrix downloaded from: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE181nnn/GSE181063/matrix/

## expression matrix
expression_data = read.table("./datasets/Lacy_2020Blood/GSE181063_series_matrix.nor.DLBCL.txt",header=T,sep = "\t")
rownames(expression_data) <- expression_data$ID_REF
expression_data <- expression_data[,-1]

## match the expression matrix and patient data
group_data <- read.table("GSE181063_series_matrix.sample.info.fulltable.txt",header=T,fill = TRUE,sep = "\t")
age_info <- subset(group_data, disease == "DLBCL" & PFS_time != "NA")

## gene annotation
library("AnnotationDbi")
library("illuminaHumanv4.db")
library("org.Hs.eg.db")
columns(illuminaHumanv4.db)
library(dplyr)
library(tibble)

id_str <- rownames(expression_data)
annotable = AnnotationDbi::select(illuminaHumanv4.db,id_str,c("SYMBOL"), keytype="PROBEID")

expression_data_df <- expression_data %>%  tibble::rownames_to_column(var = "PROBEID")
expression_data_df <- expression_data_df %>% left_join(annotable, by = "PROBEID")

expression_data_nor <- expression_data_df %>% filter(!is.na(SYMBOL))

## collapse multiple probes with mean
expression_data_nor <- expression_data_nor %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), 
   ~ if(n() > 1) mean(.x, na.rm = TRUE) else .x[1])) %>% ungroup()

expression_data_nor <- expression_data_nor %>% column_to_rownames(var = "SYMBOL")

## Ensure sample names match
common_samples <- intersect(colnames(expression_data_nor), age_info$GSE_acc)
expression_data_nor <- expression_data_nor[,common_samples]
age_info <- age_info %>% filter(GSE_acc %in% common_samples)

## order the metadata by the expression matrix
age_info <- age_info[match(colnames(expression_data_nor), age_info$GSE_acc), ]


##############################################################################
## after dataset processing, codes for downstream analyses were the same across cohorts.

##############################################################################
## Natural Cubic spline to identify cAge-related genes
##############################################################################

##### ns regression to estimate significance p-values #####

## design with degree of freedom df = 3 
df_age <- 3                                        
design <- model.matrix(~ ns(Age, df = df_age), data = age_info)

fit <- lmFit(expression_data_nor, design)
fit <- eBayes(fit, trend = TRUE)

topF <- topTableF(fit, number = Inf, sort.by = "none")  


##### calculate the delta of gene expression #####

## calculate the delta of gene expression between patients with cAge at 85% and 15% percentile
ages_ref <- as.numeric(quantile(age_info$Age, probs = c(0.15, 0.85), na.rm = TRUE))

ref_row <- lapply(age_info, function(col) {
  if (is.numeric(col)) median(col, na.rm = TRUE)
  else levels(as.factor(col))[1]
})
newdata <- as.data.frame(ref_row, stringsAsFactors = FALSE)[rep(1, 2), ]
newdata$age <- ages_ref   

## design matrix for expression prediction.
X_pred <- model.matrix(~ ns(age, df = df_age), data = newdata)

## expression prediction and delta calculation
pred_expr <- tcrossprod(X_pred, fit$coefficients)   # 2 × gene
expr_Q1   <- pred_expr[1, ]                         # log2 CPM
expr_Q3   <- pred_expr[2, ]
delta     <- expr_Q3 - expr_Q1
direction <- ifelse(delta > 0, "positive", "negative")


##### results of ns regression #####
res <- data.frame(
  gene      = rownames(v$E),
  F_stat    = topF$F,
  P_value   = topF$P.Value,
  FDR       = topF$adj.P.Val,
  expr_Q1   = expr_Q1,
  expr_Q3   = expr_Q3,
  delta     = delta,
  direction = direction,
  row.names = NULL,
  check.names = FALSE
)

## ranking
res <- res %>% arrange(FDR, desc(delta))

## Ruijin cohort, Schmitz et al NEJM 2018, and Lacy et al Blood 2020: cut-off for differenitally expressed genes: |log2 delta| >1 and FDR <0.01
age_up_genes = subset(res, FDR < 0.01 & delta > 1)$gene
age_dn_genes = subset(res, FDR < 0.01 & delta < -1)$gene

## For microarray datasets, the delta values were relative lower, we used different cut-off to obtain similar number of differential genes
## Sha et al JCO 2018, Chapuy et al NatMed 2018, and : |log2 delta| >0.585 (1.5-fold) and FDR <0.01
age_up_genes = subset(res, FDR < 0.01 & delta >0.585)$gene
age_dn_genes = subset(res, FDR < 0.01 & delta <(-0.585))$gene


################################################################
# Identify PFS-related genes (Cox analysis)
################################################################

##### cox regression analysis #####

## cox regression across all expressed genes
cox_models <- apply(t(expression_data_nor), 2, function(g) {
  cox_model <- coxph(Surv(age_info$PFS2, age_info$PFSS2) ~ g)
  hr <- summary(cox_model)$coefficients[2]
  p <- summary(cox_model)$coefficients[5]
  return(c(HR = hr, P = p))
})
cox_df <- as.data.frame(t(cox_models))
cox_df$FDR <- p.adjust(cox_df$P, method = "fdr")
cox_df <- cox_df %>% arrange(P, desc(HR))

## upregulated genes <> Hazard ratio (HR) >1
pfs_genes <- row.names(subset(cox_df, cox_df$P < 0.05 & cox_df$HR > 1))
shared_genes <- intersect(pfs_genes, age_up_genes)

## downregulated genes <> Hazard ratio (HR) <1
pfs_genes_lt1 <- row.names(subset(cox_df, cox_df$P < 0.05 & cox_df$HR < 1))
shared_genes_lt1 <- intersect(pfs_genes_lt1,age_dn_genes)


##### data frame with ns and cox regression results #####

# Subset data to shared genes 
expression_selected <- t(expression_data_nor[shared_genes,])

# Save and merge with clinical info
age_results <- data.frame(
  SampleID = rownames(expression_selected),
  ActualAge = age_info$Age,   # cAge
  PFS_time = age_info$PFS,
  PFS_event = age_info$PFSS,
  Gender = age_info$Gender,
  IPI = age_info$IPI
)


################################################################
# visualize the gene expression changes with cAge
################################################################

## extract expression of selected genes, up-regulated genes
expression_selected_ag <- as.data.frame(t(expression_data_nor[age_up_genes,]))
expression_selected_ag$SampleID = rownames(expression_selected_ag)

expression_selected_ag = left_join(expression_selected_ag, age_results[, c("SampleID", "ActualAge")], , by = c("SampleID" = "SampleID"))
expression_selected_ag$SampleID <- NULL


library(tidyr)
data_long <- expression_selected_ag %>% pivot_longer(
    cols = -ActualAge,
    names_to = "Gene",
    values_to = "Expression"
  )

data_scaled <- data_long %>%
  group_by(Gene) %>% 
  mutate( Expression_scaled = (Expression - mean(Expression)) / sd(Expression)
  ) %>% ungroup() 

pdf("cAge-gene-up.corr.pdf",height=4,width=4)
ggplot(data_scaled, aes(x = ActualAge, y = Expression_scaled)) +
  geom_smooth(
    aes(group = Gene), 
    # method = "gam",     
    method = "loess",    
    se = FALSE,        
    color = "grey",    
    size = 0.5,        
    alpha = 0.5        
  ) +
  geom_smooth(
    method = "loess",
    se = FALSE,        
    color = "blue",    
    size = 1.5        
  ) +
  theme_pubr() + labs(
    title = "",
    x = "cAge",
    y = "Z-score of expression"
  ) + theme(legend.position = "none")
dev.off()


## extract expression of selected genes, down-regulated genes

expression_selected_ag <- as.data.frame(t(expression_data_nor[age_dn_genes,]))
expression_selected_ag$SampleID = rownames(expression_selected_ag)

expression_selected_ag = left_join(expression_selected_ag, age_results[, c("SampleID", "ActualAge")], , by = c("SampleID" = "SampleID"))
expression_selected_ag$SampleID <- NULL

## transform the table
library(tidyr)
data_long <- expression_selected_ag %>% pivot_longer(
    cols = -ActualAge,
    names_to = "Gene",
    values_to = "Expression"
  )
data_scaled <- data_long %>%
  group_by(Gene) %>% 
  mutate( Expression_scaled = (Expression - mean(Expression)) / sd(Expression)
  ) %>% ungroup() 

## ploting 
pdf("cAge-gene-dn.corr.pdf",height=4,width=4)
ggplot(data_scaled, aes(x = ActualAge, y = Expression_scaled)) +
  geom_smooth(
    aes(group = Gene), 
    # method = "gam",     
    method = "loess",    
    se = FALSE,        
    color = "grey",    
    size = 0.5,        
    alpha = 0.5        
  ) +
  geom_smooth(
    method = "loess",
    se = FALSE,        
    color = "blue",    
    size = 1.5        
  ) +
  theme_pubr() + labs(
    title = "",
    x = "cAge",
    y = "Z-score of expression"
  ) + theme(legend.position = "none")
dev.off()


################################################################
## gene set scoring using ssGSEA method
################################################################

##### GSEA analysis ##### 

## Prepare expression matrix: gene x sample
expr_gsva <- as.matrix(expression_data_nor)

## gene sets
aging_gene_set_list <- list(AgingShared = shared_genes, AgingShared_lt1 = shared_genes_lt1, Aging_Up = age_up_genes, Aging_Dn = age_dn_genes)

## ssGSEA using the gsva package
ssgseaPar <- ssgseaParam(as.matrix(expr_gsva), aging_gene_set_list, normalize = T)
ssgsea_result <- gsva(ssgseaPar, verbose=FALSE)

## enrichment score for cAge-related genes
age_results$Aging_Up <- as.numeric(ssgsea_result["Aging_Up", ])
age_results$Aging_Dn <- as.numeric(ssgsea_result["Aging_Dn", ])

age_results$ssGSEA_Score <- as.numeric(ssgsea_result["AgingShared", ])
age_results$ssGSEA_Score_lt1 <- as.numeric(ssgsea_result["AgingShared_lt1", ])

## calculate composite score with cox weight
cox2 <- coxph(Surv(PFS_time, PFS_event) ~ Aging_Up + Aging_Dn, data = age_results)
beta  <- coef(cox2)[c("Aging_Up","Aging_Dn")]
age_results$Aging_cScore <- beta["Aging_Up"]*age_results$Aging_Up + beta["Aging_Dn"]*age_results$Aging_Dn

cox2 <- coxph(Surv(PFS_time, PFS_event) ~ ssGSEA_Score + ssGSEA_Score_lt1, data = age_results)
beta  <- coef(cox2)[c("ssGSEA_Score","ssGSEA_Score_lt1")]
age_results$ssGSEA_cScore <- beta["ssGSEA_Score"]*age_results$ssGSEA_Score + beta["ssGSEA_Score_lt1"]*age_results$ssGSEA_Score_lt1



##### KM curve on composite scores ##### 

## cAge-related genes 

age_results$score_quartile <- cut(age_results$Aging_cScore,
                                   breaks = quantile(age_results$Aging_cScore, probs = seq(0, 1, 0.25)),
                                   include.lowest = TRUE,
                                   labels = c("Q1", "Q2", "Q3", "Q4"))
surv_obj <- Surv(age_results$PFS_time, age_results$PFS_event)
fit_km4 <- survfit(surv_obj ~ score_quartile, data = age_results)
ggsurvplot(fit_km4, data = age_results, pval = TRUE, title = "KM: Aging Risk Score (Quartiles)", risk.table = TRUE)

pdf("Aging.ns-Aged-related-gene-cScore.PFS.pdf",height=5.3,width=4.2)
plot.new()
p = ggsurvplot(fit_km4, data = age_results, pval = TRUE, risk.table = TRUE,
           title = "PFS by ns-Aging_Score",
           xlab = "Time (Years)", ylab = "PFS (%)",
           legend.title = "Group",
           legend.labs = c("Q1","Q2","Q3", "Q4"),
           legend = c(0.8, 0.6), font.legend=12,
           font.x = c(12),
           font.y = c(12),
           font.tickslab = c(12),
           risk.table.height = 0.25,
           risk.table.fontsize = 12/.pt,
           palette = "Dark2",
          tables.theme = theme_cleantable() + 
          theme(plot.title = element_text(size = 12))
        )
p$plot <- p$plot +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

print(p, newpage = FALSE)
dev.off()


## cAge-PFS-related genes

age_results$score_quartile <- cut(age_results$ssGSEA_cScore,
                                   breaks = quantile(age_results$ssGSEA_cScore, probs = seq(0, 1, 0.25)),
                                   include.lowest = TRUE,
                                   labels = c("Q1", "Q2", "Q3", "Q4"))
surv_obj <- Surv(age_results$PFS_time, age_results$PFS_event)
fit_km4 <- survfit(surv_obj ~ score_quartile, data = age_results)
ggsurvplot(fit_km4, data = age_results, pval = TRUE, title = "KM: Aging Risk Score (Quartiles)", risk.table = TRUE)

pdf("Aging.GSEA.composite_Score.PFS.pdf",height=5.3,width=4.2)
plot.new()
p = ggsurvplot(fit_km4, data = age_results, pval = TRUE, risk.table = TRUE,
           title = "PFS by cAge composite Score",
           xlab = "Time (Years)", ylab = "PFS (%)",
           legend.title = "Group",
           legend.labs = c("Q1","Q2","Q3", "Q4"),
           legend = c(0.8, 0.6), font.legend=12,
           font.x = c(12),
           font.y = c(12),
           font.tickslab = c(12),
           risk.table.height = 0.25,
           risk.table.fontsize = 12/.pt,
           palette = "Dark2",
          tables.theme = theme_cleantable() + 
          theme(plot.title = element_text(size = 12))
        )
p$plot <- p$plot +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

print(p, newpage = FALSE)
dev.off()



################################################################
## molecular age estimation
################################################################

##### R package ##### 

library(glmnet);               
library(survival)    
library(rms)          
library(timeROC)      
library(survIDINRI)   


##### beta coefficient ##### 

# only include cases with PFS data
age_results_sub = subset(age_results, !is.na(PFS_time))
y <- with(age_results_sub, Surv(PFS_time, PFS_event))

cox_age   <- coxph(y ~ ActualAge,  data = age_results_sub)
cox_score <- coxph(y ~ ssGSEA_cScore, data = age_results_sub)
beta_age   <- coef(cox_age)["ActualAge"]
beta_score <- coef(cox_score)["ssGSEA_cScore"]


##### molecular age estimation ##### 


## shift the score to all great than zero as representative of mAge
score_shift   <- age_results_sub$ssGSEA_cScore - min(age_results_sub$ssGSEA_cScore)
beta_ratio    <- beta_score / beta_age             
age_results_sub$mAge <- score_shift * beta_ratio

## median-anchored the mAge based on actual cAge
scale_factor <- IQR(age_results_sub$ActualAge) / IQR(age_results_sub$mAge)
shift_factor <- median(age_results_sub$ActualAge) - median(age_results_sub$mAge * scale_factor)
age_results_sub$mAge_aligned <- age_results_sub$mAge * scale_factor + shift_factor


##### KM curves with cAge and mAge ##### 

## KM curve on actual cAge group
age_results_sub$Age_group <- as.factor(cut(age_results_sub$ActualAge,
                             # breaks = c(0, 40, 50, 60, 70, 80, 100),
                             breaks = c(0, 44, 59, 74, 150),
                             labels = c("lt_45","45-59","60-74", "gt_75"),
                             right = FALSE))
table(age_results_sub$Age_group)

surv_obj <- Surv(age_results_sub$PFS_time, age_results_sub$PFS_event)
fit_km <- survfit(surv_obj ~ Age_group, data = age_results_sub)
ggsurvplot(fit_km, data = age_results_sub, pval = TRUE, title = "KM: Aging Risk Score", risk.table = TRUE)

pdf("Aging.cAge-groups.PFS.pdf",height=5.3,width=4.2)
plot.new()
p = ggsurvplot(fit_km, data = age_results_sub, pval = TRUE, risk.table = TRUE,
           title = "PFS by cAge group",
           xlab = "Time (Years)", ylab = "PFS (%)",
           legend.title = "Group",
           legend.labs = c("lt_45","45-59","60-74", "gt_75"),
           legend = c(0.8, 0.6), font.legend=12,
           font.x = c(12),
           font.y = c(12),
           font.tickslab = c(12),
           risk.table.height = 0.25,
           risk.table.fontsize = 12/.pt,
           palette = "Dark2",
          tables.theme = theme_cleantable() + 
          theme(plot.title = element_text(size = 12))
        )
p$plot <- p$plot +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

print(p, newpage = FALSE)
dev.off()


## KM curves on actual mAge group
## PFS
age_results_sub$Age_group <- as.factor(cut(age_results_sub$mAge_aligned,
                             # breaks = c(0, 45, 55, 65, 75, 150),
                             breaks = c(0, 44, 59, 74, 150),
                             labels = c("lt_45","45-59","60-74", "gt_75"),
                             right = FALSE))
table(age_results_sub$Age_group)

surv_obj <- Surv(age_results_sub$PFS_time, age_results_sub$PFS_event)
fit_km <- survfit(surv_obj ~ Age_group, data = age_results_sub)
ggsurvplot(fit_km, data = age_results_sub, pval = TRUE, title = "KM: Aging Risk Score", risk.table = TRUE)

pdf("Aging.mAge-groups.PFS.pdf",height=5.3,width=4.2)
plot.new()
p = ggsurvplot(fit_km, data = age_results_sub, pval = TRUE, risk.table = TRUE,
           title = "PFS by cAge group",
           xlab = "Time (Years)", ylab = "PFS (%)",
           legend.title = "Group",
           legend.labs = c("lt_45","45-59","60-74", "gt_75"),
           legend = c(0.8, 0.6), font.legend=12,
           font.x = c(12),
           font.y = c(12),
           font.tickslab = c(12),
           risk.table.height = 0.25,
           risk.table.fontsize = 12/.pt,
           palette = "Dark2",
          tables.theme = theme_cleantable() + 
          theme(plot.title = element_text(size = 12))
        )
p$plot <- p$plot +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

print(p, newpage = FALSE)
dev.off()

## OS

age_results_sub$Age_group <- as.factor(cut(age_results_sub$mAge_aligned,
                             # breaks = c(0, 45, 55, 65, 75, 150),
                             breaks = c(0, 44, 59, 74, 150),
                             labels = c("lt_45","45-59","60-74", "gt_75"),
                             right = FALSE))
table(age_results_sub$Age_group)

surv_obj <- Surv(age_results_sub$OS, age_results_sub$OSS)
fit_km <- survfit(surv_obj ~ Age_group, data = age_results_sub)
ggsurvplot(fit_km, data = age_results_sub, pval = TRUE, title = "KM: Aging Risk Score", risk.table = TRUE)

pdf("Aging.mAge-groups.OS.pdf",height=5.3,width=4.2)
plot.new()
p = ggsurvplot(fit_km, data = age_results_sub, pval = TRUE, risk.table = TRUE,
           title = "OS by cAge group",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "Group",
           legend.labs = c("lt_45","45-59","60-74", "gt_75"),
           legend = c(0.8, 0.6), font.legend=12,
           font.x = c(12),
           font.y = c(12),
           font.tickslab = c(12),
           risk.table.height = 0.25,
           risk.table.fontsize = 12/.pt,
           palette = "Dark2",
          tables.theme = theme_cleantable() + 
          theme(plot.title = element_text(size = 12))
        )
p$plot <- p$plot +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

print(p, newpage = FALSE)
dev.off()


################################################################
## Analysis of mAge acceleration relative to cAge
################################################################

library(dplyr)
library(ggplot2)
library(broom)
library(splines)


##### mAge acceleration relative to cAge ##### 

## baseline fitted using natural splines

fit <- lm(mAge_aligned ~ ns(ActualAge, df = 3), data = age_results_sub)
age_results_sub <- age_results_sub %>% 
         mutate(pred_mAge = predict(fit, se.fit = TRUE)$fit,
                pred_mAge.se = predict(fit, se.fit = TRUE)$se.fit,
                residual  = (mAge_aligned - pred_mAge)/pred_mAge.se,
                deltaAge  = mAge_aligned - ActualAge)
hist(age_results_sub$residual,breaks = 50)


## studentized residual 
age_results_sub$stud_res <- rstudent(fit)
hist(age_results_sub$stud_res,breaks = 50)


## define mAge acceleration and deceleration using a residual of 1 and -1, respectively
age_results_sub$mAge_group = ifelse(age_results_sub$stud_res > 1, "Acc", ifelse(age_results_sub$stud_res < -1, "Dec", "Syn"))

age_results_sub$mAge_group = factor(age_results_sub$mAge_group, levels = c("Dec","Syn","Acc"))
table(age_results_sub$mAge_group)



##### merge the results from DLBCL classification analysis ##### 

## COO
age_results_sub = left_join(age_results_sub, res_coo[, c("Sample","COO_est")], , by = c("SampleID" = "Sample"))
## LME
age_results_sub = left_join(age_results_sub, lme_predicted[, c("Sample","LME")], , by = c("SampleID" = "Sample"))



##### visualize mAge acceleration ##### 

library(broom)  

## mAge versus cAge in a scatter plot
p_scatter <- ggplot(age_results_sub,
             aes(x = ActualAge, y = mAge_aligned, colour = mAge_group)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", aes(group = 1),
              formula = y ~ splines::ns(x, df = 3),
              se = FALSE, linetype = "dashed",color="black") +
  scale_colour_manual(values = c(Acc = "#BF360C90",
                                 Dec = "#4CAF5090",
                                 Syn = "#0091EA15")) +
  labs(title = "Chronological vs Molecular age",
       x = "Chronological age (years)",
       y = "Molecular age (years)",
       colour = NULL) +
  theme_pubr(base_size = 12)

pdf("cAge-mAge.model.correlation.pdf", height=5, width=4.5)
p_scatter 
dev.off()


## frequency of mAge_group by cAge groups
age_results_sub$cAge_pct <- as.factor(cut(age_results_sub$ActualAge,
                             # breaks = c(0, 45, 55, 65, 75, 150),
                             breaks = c(0, 44, 59, 74, 150),
                             labels = c("<45","45-59","60-74", "≥75"),
                             right = FALSE))
table(age_results_sub$cAge_pct)


prop_tbl <- age_results_sub %>%
    group_by(cAge_pct, mAge_group) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(cAge_pct) %>%
    mutate(prop = n / sum(n))   
p <- ggplot(prop_tbl,
                     aes(x = cAge_pct, y = prop, fill = mAge_group)) +
    geom_col(position = "fill", width = 0.9) +   # fill = 100%，条上下堆叠
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c(Acc = "#BF360C98",
                                 Dec = "#4CAF5098",
                                 Syn = "#0091EA98")) +
    labs(title = "",
         x = "",
         y = "Proportion of patients",
         fill = NULL) +
    theme_pubr(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("cAge_pct.Acc_group.pdf",p,width = 2.8,height = 5)


## Age density plot

p1 <- ggplot(age_results_sub, aes(x=ActualAge)) + geom_density() +
    geom_vline(aes(xintercept=median(ActualAge)),
               color="blue", linetype="dashed", size=1) +
    theme_pubr(base_size = 12) +
    labs(title = "Chronological age", x = "cAge (Years)") + xlim(20,100)

p2 <- ggplot(age_results_sub, aes(x=mAge_aligned)) + geom_density() +
    geom_vline(aes(xintercept=median(mAge_aligned)),
               color="blue", linetype="dashed", size=1) +
    theme_pubr(base_size = 12) +
    labs(title = "Molecular age", x = "mAge (Years)") + xlim(20,100)

pdf("cAge-mAge.density.pdf",height=6,width=4)
p1 / p2
dev.off()



##### cAge distributions vs mAge distributions ##### 

age_results_sub$mAge_pct <- as.factor(cut(age_results_sub$mAge_aligned,
                             # breaks = c(0, 45, 55, 65, 75, 150),
                             breaks = c(0, 44, 59, 74, 150),
                             labels = c("lt_45","45-59","60-74", "gt_75"),
                             right = FALSE))
table(age_results_sub$mAge_pct)

age_results_sub$cAge_pct <- as.factor(cut(age_results_sub$ActualAge,
                             # breaks = c(0, 45, 55, 65, 75, 150),
                             breaks = c(0, 44, 59, 74, 150),
                             labels = c("lt_45","45-59","60-74", "gt_75"),
                             right = FALSE))
table(age_results_sub$cAge_pct)

allu <- age_results_sub %>% group_by(mAge_pct, cAge_pct) %>% summarise(Freq = n())  

p = ggplot(data = allu, aes(axis1 = cAge_pct, axis2 = mAge_pct, y = Freq)) +
    geom_alluvium(aes(fill = mAge_pct)) +
    geom_stratum(aes(fill = cAge_pct)) +
    scale_fill_brewer(palette="Dark2") +
    geom_text(stat = "stratum",aes(label = after_stat(stratum))) +
    theme_pubr(base_size = 12) + theme(legend.position = "none") + labs(y = "Number of cases") +
    scale_x_discrete(limits = c("cAge_pct", "mAge_pct"), labels=c("mAge_pct" = "mAge", "cAge_pct" = "cAge"), expand = c(0.15, 0.05))

ggsave("cAge-mAge.alluplot.Age_pct.pdf",p,width = 3,height = 4.5)


##### KM curve based on mAge group ##### 

## PFS
surv_obj <- Surv(age_results_sub$PFS_time, age_results_sub$PFS_event)
fit_km <- survfit(surv_obj ~ mAge_group, data = age_results_sub)
ggsurvplot(fit_km, data = age_results_sub, pval = TRUE, title = "KM: Aging Risk Score", risk.table = TRUE)

pdf("mAge.Acceleration.group.PFS.pdf",height=5,width=4.2)
plot.new()
p = ggsurvplot(fit_km, data = age_results_sub, pval = TRUE, risk.table = TRUE,
           title = "PFS by mAge_group",
           xlab = "Time (Years)", ylab = "PFS (%)",
           legend.title = "Group",
           legend.labs = c("Dec", "Syn","Acc"),
           legend = c(0.8, 0.2), font.legend=12,
           font.x = c(12),
           font.y = c(12),
           font.tickslab = c(12),
           risk.table.height = 0.22,
           risk.table.fontsize = 12/.pt,
           palette = c("#4CAF50", "#0091EA", "#BF360C"),
          tables.theme = theme_cleantable() + 
          theme(plot.title = element_text(size = 12))
        )
p$plot <- p$plot +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

print(p, newpage = FALSE)
dev.off()

## OS
surv_obj <- Surv((age_results_sub$OS2)/12, age_results_sub$OSS2)
fit_km <- survfit(surv_obj ~ mAge_group, data = age_results_sub)
ggsurvplot(fit_km, data = age_results_sub, pval = TRUE, title = "KM: Aging Risk Score", risk.table = TRUE)

pdf("mAge.Acceleration.group.OS.pdf",height=5,width=4.2)
plot.new()
p = ggsurvplot(fit_km, data = age_results_sub, pval = TRUE, risk.table = TRUE,
           title = "OS by mAge_group",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "Group",
           legend.labs = c("Dec", "Syn","Acc"),
           legend = c(0.8, 0.2), font.legend=12,
           font.x = c(12),
           font.y = c(12),
           font.tickslab = c(12),
           risk.table.height = 0.22,
           risk.table.fontsize = 12/.pt,
           palette = c("#4CAF50", "#0091EA", "#BF360C"),
          tables.theme = theme_cleantable() + 
          theme(plot.title = element_text(size = 12))
        )
p$plot <- p$plot +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

print(p, newpage = FALSE)
dev.off()



##### KM curve based on mAge group, within each COO subtype ##### 

## ABC
age_results_sub2 = subset(age_results_sub, COO_est == "ABC")
surv_obj = Surv(age_results_sub2$PFS_time, age_results_sub2$PFS_event)

fit_km <- survfit(surv_obj ~ mAge_group, data = age_results_sub2)
ggsurvplot(fit_km, data = age_results_sub2, pval = TRUE, title = "KM: Aging Risk Score", risk.table = TRUE)

pdf("mAge.Acceleration.group.PFS.ABC.pdf",height=4,width=3.2)
plot.new()
p = ggsurvplot(fit_km, data = age_results_sub2, pval = TRUE, risk.table = TRUE,
           title = "PFS by mAge_group",
           xlab = "Time (Years)", ylab = "PFS (%)",
           legend.title = "Group",
           legend.labs = c("Dec", "Syn","Acc"),
           legend = c(0.8, 0.2), font.legend=12,
           font.x = c(12),
           font.y = c(12),
           font.tickslab = c(12),
           risk.table.height = 0.22,
           risk.table.fontsize = 12/.pt,
           palette = c("#4CAF50", "#0091EA", "#BF360C"),
          tables.theme = theme_cleantable() + 
          theme(plot.title = element_text(size = 12))
        )
p$plot <- p$plot +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

print(p, newpage = FALSE)
dev.off()

## GCB
age_results_sub2 = subset(age_results_sub, COO_est == "GCB")
surv_obj = Surv(age_results_sub2$PFS_time, age_results_sub2$PFS_event)

fit_km <- survfit(surv_obj ~ mAge_group, data = age_results_sub2)
ggsurvplot(fit_km, data = age_results_sub2, pval = TRUE, title = "KM: Aging Risk Score", risk.table = TRUE)

pdf("mAge.Acceleration.group.GCB.pdf",height=4,width=3.2)
plot.new()
p = ggsurvplot(fit_km, data = age_results_sub2, pval = TRUE, risk.table = TRUE,
           title = "PFS by mAge_group",
           xlab = "Time (Years)", ylab = "PFS (%)",
           legend.title = "Group",
           legend.labs = c("Dec", "Syn","Acc"),
           legend = c(0.8, 0.2), font.legend=12,
           font.x = c(12),
           font.y = c(12),
           font.tickslab = c(12),
           risk.table.height = 0.22,
           risk.table.fontsize = 12/.pt,
           palette = c("#4CAF50", "#0091EA", "#BF360C"),
          tables.theme = theme_cleantable() + 
          theme(plot.title = element_text(size = 12))
        )
p$plot <- p$plot +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

print(p, newpage = FALSE)
dev.off()


##### KM curve based on mAge group, within each LME-DE(DP) or LME-ME(MS) subtype ##### 

## LME-DE
age_results_sub2 = subset(lme_coldata, LME.y == "LME-DE")
surv_obj = Surv(age_results_sub2$PFS_time, age_results_sub2$PFS_event)

fit_km <- survfit(surv_obj ~ mAge_group, data = age_results_sub2)
ggsurvplot(fit_km, data = age_results_sub2, pval = TRUE, title = "KM: Aging Risk Score", risk.table = TRUE)

pdf("mAge.Acceleration.group.PFS.LME-DE.pdf",height=4,width=3.2)
plot.new()
p = ggsurvplot(fit_km, data = age_results_sub2, pval = TRUE, risk.table = TRUE,
           title = "PFS by mAge_group",
           xlab = "Time (Years)", ylab = "PFS (%)",
           legend.title = "Group",
           legend.labs = c("Dec", "Syn","Acc"),
           legend = c(0.8, 0.2), font.legend=12,
           font.x = c(12),
           font.y = c(12),
           font.tickslab = c(12),
           risk.table.height = 0.22,
           risk.table.fontsize = 12/.pt,
           palette = c("#4CAF50", "#0091EA", "#BF360C"),
          tables.theme = theme_cleantable() + 
          theme(plot.title = element_text(size = 12))
        )
p$plot <- p$plot +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

print(p, newpage = FALSE)
dev.off()


## LME-ME
age_results_sub2 = subset(lme_coldata, LME.y == "LME-ME")
surv_obj = Surv(age_results_sub2$PFS_time, age_results_sub2$PFS_event)

fit_km <- survfit(surv_obj ~ mAge_group, data = age_results_sub2)
ggsurvplot(fit_km, data = age_results_sub2, pval = TRUE, title = "KM: Aging Risk Score", risk.table = TRUE)

pdf("mAge.Acceleration.group.PFS.LME-ME.pdf",height=4,width=3.2)
plot.new()
p = ggsurvplot(fit_km, data = age_results_sub2, pval = TRUE, risk.table = TRUE,
           title = "PFS by mAge_group",
           xlab = "Time (Years)", ylab = "PFS (%)",
           legend.title = "Group",
           legend.labs = c("Dec", "Syn","Acc"),
           legend = c(0.8, 0.2), font.legend=12,
           font.x = c(12),
           font.y = c(12),
           font.tickslab = c(12),
           risk.table.height = 0.22,
           risk.table.fontsize = 12/.pt,
           palette = c("#4CAF50", "#0091EA", "#BF360C"),
          tables.theme = theme_cleantable() + 
          theme(plot.title = element_text(size = 12))
        )
p$plot <- p$plot +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

print(p, newpage = FALSE)
dev.off()




##############################################################################
## differential expression between acc and dec patients
##############################################################################

library(edgeR)      
library(limma)     


## sample IDs for each mAge group
pt_acc = subset(age_results_sub, mAge_group == "Acc" )$SampleID
pt_dec = subset(age_results_sub, mAge_group == "Dec" )$SampleID
pt_syn = subset(age_results_sub, mAge_group == "Syn" )$SampleID

## create the design matrix for limma
group_info$mAge_group = ifelse(group_info$SampleID %in% pt_acc,"Acc",
                        ifelse(group_info$SampleID %in% pt_dec,"Dec",
                          ifelse(group_info$SampleID %in% pt_syn,"Syn",
                            "ohters"
                            )
                          )
                        )


group = factor(as.character((group_info[colnames(expression_data_nor),])$mAge_group))

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit = lmFit(expression_data_nor, design)

## limma contrasts.fit, Accelerated vs Synchronized
contrast.matrix <- makeContrasts(Acc_vs_Syn = Acc - Syn, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# differential expression analysis
res_Acc_vs_Syn <- topTable(fit2, coef = "Acc_vs_Syn", adjust = "fdr", number = Inf)
res_Acc_vs_Syn$logP <- -log10(res_Acc_vs_Syn$P.Value)
head(res_Acc_vs_Syn)


## limma contrasts.fit, Decelerated vs Synchronized
contrast.matrix <- makeContrasts(Dec_vs_Syn = Dec - Syn, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# differential expression analysis
res_Dec_vs_Syn <- topTable(fit2, coef = "Dec_vs_Syn", adjust = "fdr", number = Inf)
res_Dec_vs_Syn$logP <- -log10(res_Dec_vs_Syn$P.Value)
head(res_Dec_vs_Syn)


## volcano plot

## Accelerated vs Synchronized
res_Acc_vs_Syn$logQ <- -log10(res_Acc_vs_Syn$adj.P.Val)

res_Acc_vs_Syn$Significance <- "Not Significant"
res_Acc_vs_Syn$Significance[res_Acc_vs_Syn$logFC > 1 & res_Acc_vs_Syn$adj.P.Val < 0.05] <- "Upregulated"
res_Acc_vs_Syn$Significance[res_Acc_vs_Syn$logFC < -1 & res_Acc_vs_Syn$adj.P.Val < 0.05] <- "Downregulated"


p <- ggplot(res_Acc_vs_Syn, aes(x = logFC, y = logQ)) +
  geom_point(aes(color = Significance), alpha = 0.6, show.legend = FALSE) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme_bw(base_size=12) +
  theme(axis.text=element_text(colour="black"), aspect.ratio = 1) +
  labs(title = "Ruijin cohort (Acc vs Syn)", x = "Log2 Fold Change", y = "-Log10 Adj.P")


pdf("mAge.Acceleration.group.DE.Accel_vs_Sync.limma.volcano-plot.pdf", width = 4, height = 4)
p + geom_hline(yintercept=1.3, linetype="dashed", color = "GREY", size=1) + geom_vline(xintercept=1, linetype="dashed", color = "GREY", size=1) + geom_vline(xintercept=-1, linetype="dashed", color = "GREY", size=1)
dev.off()


## Decelerated vs Synchronized

res_Dec_vs_Syn$logQ <- -log10(res_Dec_vs_Syn$adj.P.Val)

res_Dec_vs_Syn$Significance <- "Not Significant"
res_Dec_vs_Syn$Significance[res_Dec_vs_Syn$logFC > 1 & res_Dec_vs_Syn$adj.P.Val < 0.05] <- "Upregulated"
res_Dec_vs_Syn$Significance[res_Dec_vs_Syn$logFC < -1 & res_Dec_vs_Syn$adj.P.Val < 0.05] <- "Downregulated"


p <- ggplot(res_Dec_vs_Syn, aes(x = logFC, y = logQ)) +
  geom_point(aes(color = Significance), alpha = 0.6, show.legend = FALSE) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme_bw(base_size=12) +
  theme(axis.text=element_text(colour="black"), aspect.ratio = 1) +
  labs(title = "Ruijin cohort (Acc vs Syn)", x = "Log2 Fold Change", y = "-Log10 Adj.P")


pdf("mAge.Acceleration.group.DE.Decel_vs_Sync.limma.volcano-plot.pdf", width = 4, height = 4)
p + geom_hline(yintercept=1.3, linetype="dashed", color = "GREY", size=1) + geom_vline(xintercept=1, linetype="dashed", color = "GREY", size=1) + geom_vline(xintercept=-1, linetype="dashed", color = "GREY", size=1)
dev.off()



##############################################################################
## correlation analysis, mAge groups with DLBCL classifications
##############################################################################

library(ggalluvial)

## COO
age_results_sub2 = subset(age_results_sub, !is.na(COO_est))
chisq.test(age_results_sub2$COO_est, age_results_sub2$mAge_group)
# X-squared = 53.49, df = 4, p-value = 6.73e-11

allu <- age_results_sub2 %>% group_by(mAge_group, COO_est) %>% summarise(Freq = n())  
allu_percent <- allu %>% group_by(mAge_group) %>%            
  mutate(total = sum(Freq)) %>% 
  ungroup() %>% mutate(percent = Freq / total * 100) 

p = ggplot(allu_percent, aes(x = mAge_group, stratum = COO_est, alluvium = COO_est, y = percent, fill = COO_est)) +
    geom_flow(width = 0.3) +
    geom_stratum(width = 0.5) +
    scale_fill_manual(values=c("#EE5C42","#7CCD7C","grey")) +
    theme_pubr(base_size = 12) +
    labs(title = "",
         x = "",
         y = "Percentage of cases")+ 
    theme(legend.position=c("right")) + guides(fill=guide_legend(title="COO"))

ggsave("Alluplot.mAge-Acceleration-group.COO_est.pdf",p,width = 4.2,height = 3)


## LME subtypes

age_results_sub2 = subset(age_results_sub, !is.na(LME))
chisq.test(age_results_sub2$LME, age_results_sub2$mAge_group)
# X-squared = 195.57, df = 6, p-value < 2.2e-16

allu <- age_results_sub2 %>% group_by(mAge_group, LME) %>% summarise(Freq = n())  
allu_percent <- allu %>% group_by(mAge_group) %>%            
  mutate(total = sum(Freq)) %>% 
  ungroup() %>% mutate(percent = Freq / total * 100) 


p = ggplot(allu_percent, aes(x = mAge_group, stratum = LME, alluvium = LME, y = percent, fill = LME)) +
    geom_flow(width = 0.3) +
    geom_stratum(width = 0.5) +
    scale_fill_manual(values = c("#F39B7FB2","#76afda","#F1CC2F","#00A087B2")) +
    theme_pubr(base_size = 12) +
    labs(title = "",
         x = "",
         y = "Percentage of cases")+ 
    theme(legend.position=c("right")) + guides(fill=guide_legend(title="LME"))

ggsave("Alluplot.mAge-Acceleration-group.LME.pdf",p,width = 4.2,height = 3)



##############################################################################
## correlation analysis, mAge groups with PFS24 event incidences
##############################################################################

## PFS24 event incidence by mAge groups, in all patients

pfs24 = subset(age_results_sub, PFS_time <= 24 & PFS_event == 1)$SampleID
age_results_sub$PFS24 = ifelse(age_results_sub$SampleID %in% pfs24, "Event", "No Event")

chisq.test(age_results_sub$PFS24, age_results_sub$mAge_group)
# X-squared = 41.085, df = 2, p-value = 1.198e-09

allu <- age_results_sub %>% group_by(mAge_group, PFS24) %>% summarise(Freq = n())  
allu_percent <- allu %>% group_by(mAge_group) %>%            
  mutate(total = sum(Freq)) %>% 
  ungroup() %>% mutate(percent = Freq / total * 100) 


p = ggplot(allu_percent, aes(x = mAge_group, stratum = PFS24, alluvium = PFS24, y = percent, fill = PFS24)) +
    geom_flow(width = 0.3) +
    geom_stratum(width = 0.4) +
    scale_fill_brewer(palette = "Set1") +
    theme_pubr(base_size = 12) +
    labs(title = "",
         x = "",
         y = "Number of cases")+ 
    theme(legend.position=c("right")) + guides(fill=guide_legend(title="PFS24"))

ggsave("Alluplot.mAge-Acceleration-group.PFS24.pdf",p,width = 4.2,height = 4)


## PFS24 event incidence by mAge groups, in patients with cAge <60

age_results_sub2 = subset(age_results_sub, ActualAge < 60)
age_results_sub2$PFS24 = ifelse(age_results_sub2$SampleID %in% pfs24, "Event", "No Event")
chisq.test(age_results_sub2$PFS24, age_results_sub2$mAge_group)
# X-squared = 14.684, df = 2, p-value = 0.0006478

allu <- age_results_sub2 %>% group_by(mAge_group, PFS24) %>% summarise(Freq = n())  
allu_percent <- allu %>% group_by(mAge_group) %>%            
  mutate(total = sum(Freq)) %>% 
  ungroup() %>% mutate(percent = Freq / total * 100) 

p = ggplot(allu_percent, aes(x = mAge_group, stratum = PFS24, alluvium = PFS24, y = percent, fill = PFS24)) +
    geom_flow(width = 0.3) +
    geom_stratum(width = 0.4) +
    scale_fill_brewer(palette = "Set1") +
    theme_pubr(base_size = 12) +
    labs(title = "",
         x = "",
         y = "Number of cases")+ 
    theme(legend.position=c("right")) + guides(fill=guide_legend(title="PFS24"))

ggsave("Alluplot.mAge-Acceleration-group.PFS24.lt60.pdf",p,width = 4.2,height = 4)


## PFS24 event incidence by mAge groups, in patients with cAge >= 60

age_results_sub2 = subset(age_results_sub, ActualAge >= 60)
age_results_sub2$PFS24 = ifelse(age_results_sub2$SampleID %in% pfs24, "Event", "No Event")
chisq.test(age_results_sub2$PFS24, age_results_sub2$mAge_group)
# X-squared = 26.63, df = 2, p-value = 1.65e-06

allu <- age_results_sub2 %>% group_by(mAge_group, PFS24) %>% summarise(Freq = n())  
allu_percent <- allu %>% group_by(mAge_group) %>%            
  mutate(total = sum(Freq)) %>% 
  ungroup() %>% mutate(percent = Freq / total * 100) 

p = ggplot(allu_percent, aes(x = mAge_group, stratum = PFS24, alluvium = PFS24, y = percent, fill = PFS24)) +
    geom_flow(width = 0.3) +
    geom_stratum(width = 0.4) +
    scale_fill_brewer(palette = "Set1") +
    theme_pubr(base_size = 12) +
    labs(title = "",
         x = "",
         y = "Number of cases")+ 
    theme(legend.position=c("right")) + guides(fill=guide_legend(title="PFS24"))

ggsave(".Alluplot.mAge-Acceleration-group.PFS24.gt60.pdf",p,width = 4.2,height = 4)



##############################################################################
## multivariate cox regression analysis, with PFS
##############################################################################

## overall PFS

pt_fisher = age_results_sub

pt_fisher$cAge_gt_60 <- ifelse(pt_fisher$ActualAge == "",NA, ifelse(pt_fisher$ActualAge >= 60, "1", "0"))
pt_fisher$mAge_gt_60 <- ifelse(pt_fisher$mAge_aligned == "",NA, ifelse(pt_fisher$mAge_aligned >= 60, "1", "0"))
pt_fisher$Gender_gp <- ifelse(pt_fisher$Gender == "",NA, ifelse(pt_fisher$Gender == "Male", "1", "0"))
pt_fisher$IPI_3to5 <- ifelse(pt_fisher$IPI == "",NA, ifelse(pt_fisher$IPI >2,"1", "0"))
pt_fisher$COO_est = factor(pt_fisher$COO_est, level = c("GCB","ABC","Unclassified"))
pt_fisher$mAge_group = factor(pt_fisher$mAge_group, level = c("Syn","Acc","Dec"))

cox_model <- coxph(Surv(PFS_time, PFS_event) ~ mAge_group +cAge_gt_60 +Gender_gp +COO_est +IPI_3to5, data = pt_fisher)
cox_summary = summary(cox_model)
cox_summary


## PFS24
pt_fisher$PFS24_status = as.numeric(ifelse(pt_fisher$PFS24 == "Event", "1", "0"))
pt_fisher$PFS24_time = ifelse(pt_fisher$PFS_time <= 2, pt_fisher$PFS_time, 2)

cox_model <- coxph(Surv(PFS24_time, PFS24_status) ~ mAge_group +cAge_gt_60 +Gender_gp +COO_est +IPI_3to5, data = pt_fisher)
cox_summary = summary(cox_model)
cox_summary

## PFS12
pfs12 = subset(pt_fisher, PFS_time <= 1 & PFS_event == 1)$SampleID
pt_fisher$PFS12 = ifelse(pt_fisher$SampleID %in% pfs12, "Event", "No Event")
pt_fisher$PFS12_status = as.numeric(ifelse(pt_fisher$PFS12 == "Event", "1", "0"))
pt_fisher$PFS12_time = ifelse(pt_fisher$PFS_time <= 1, pt_fisher$PFS_time, 1)

cox_model <- coxph(Surv(PFS12_time, PFS12_status) ~ mAge_group +cAge_gt_60 +Gender_gp +COO_est +IPI_3to5, data = pt_fisher)
cox_summary = summary(cox_model)
cox_summary


## correlation analysis with clinical baseline variables
selected_columns <- c("cAge_gt_60", "Gender_gp", "COO_est","ECOG_gt2","AA_34","LDH","Extranodal_gt2","IPI_3to5")

## Acc vs Syn
pt_fisher_sub = subset(pt_fisher, mAge_group %in% c("Syn","Acc"))
pt_fisher_sub$mAge_group =factor(pt_fisher_sub$mAge_group, levels=c("Syn","Acc"))

results <- lapply(selected_columns, function(col) {
  pt_fisher_sub[[col]] = as.factor(pt_fisher_sub[[col]])
  contingency_table <- table(pt_fisher_sub$mAge_group, pt_fisher_sub[[col]])
  #test <- chisq.test(contingency_table)
  test <- fisher.test(contingency_table)

  data.frame(
    Column = col,
    P_Value = test$p.value,
    Odd_ratio = as.numeric(test$estimate),
    OR_low = test$conf.int[1],
    OR_high = test$conf.int[2]
  )
})

print(do.call(rbind, results))

## Dec vs Syn
pt_fisher_sub = subset(pt_fisher, mAge_group %in% c("Syn","Dec"))
pt_fisher_sub$mAge_group =factor(pt_fisher_sub$mAge_group, levels=c("Syn","Dec"))

results <- lapply(selected_columns, function(col) {
  pt_fisher_sub[[col]] = as.factor(pt_fisher_sub[[col]])
  contingency_table <- table(pt_fisher_sub$mAge_group, pt_fisher_sub[[col]])
  #test <- chisq.test(contingency_table)
  test <- fisher.test(contingency_table)

  data.frame(
    Column = col,
    P_Value = test$p.value,
    Odd_ratio = as.numeric(test$estimate),
    OR_low = test$conf.int[1],
    OR_high = test$conf.int[2]
  )
})
print(do.call(rbind, results))



##############################################################################
## ecotyper results
##############################################################################

## read ecotyper outputs, and create table with metadata and ecotypes
ecotyper_types = read.table("./ecotyper_output/Ecotypes/ecotype_assignment.txt",header=T,sep="\t")
coldata = merge(age_results_sub, ecotyper_types, by.x = "SampleID", by.y = "ID")


##### correlation analysis #####

chisq.test(table(coldata$mAge_group, coldata$Ecotype))
# X-squared = 355.67, df = 16, p-value < 2.2e-16

allu <- coldata %>% group_by(mAge_group, Ecotype) %>% summarise(Freq = n())  
allu_percent <- allu %>% group_by(mAge_group) %>%            
  mutate(total = sum(Freq)) %>% 
  ungroup() %>% mutate(percent = Freq / total * 100) 

p = ggplot(allu_percent, aes(x = mAge_group, stratum = Ecotype, alluvium = Ecotype, y = percent, fill = Ecotype)) +
    geom_flow(width = 0.3) +
    geom_stratum(width = 0.5) +
    scale_fill_brewer(palette = "Spectral") +
    theme_pubr(base_size = 12) +
    labs(title = "",
         x = "",
         y = "Percentage of cases")+ 
    theme(legend.position=c("right")) + guides(fill=guide_legend(title="Ecotypes"))

ggsave("Alluplot.mAge-Acceleration-group.Ecotype.pdf",p,width = 3.6,height = 4)


##### cox with ecotypes and LME subtypes #####

pt_fisher_eco = left_join(pt_fisher, ecotyper_types[, c("ID", "Ecotype")], , by = c("SampleID" = "ID"))
pt_fisher_eco = left_join(pt_fisher_eco, age_results_sub[, c("SampleID", "LME")], , by = c("SampleID" = "SampleID"))


cox_model <- coxph(Surv(PFS_time, PFS_event) ~ mAge_group +cAge_gt_60 +Gender_gp +COO_est +IPI_3to5 +Ecotype +LME.y, data = pt_fisher_eco)
cox_summary = summary(cox_model)
cox_summary

cox_model <- coxph(Surv(PFS24_time, PFS24_status) ~ mAge_group +cAge_gt_60 +Gender_gp +COO_est +IPI_3to5 +Ecotype +LME.y, data = pt_fisher_eco)
cox_summary = summary(cox_model)
cox_summary



##############################################################################
## Survival analysis with the signature for mAge acceleration
##############################################################################

## Signature score using ssGSEA analysis
gset_acc <- getGmt("./2025.08.Aging.Datasets.ACC.Dec.DEGs.gmt")
gset_acc <- subsetGeneSets(gset_acc, rownames(expression_data_nor)) 

ssgseaPar <- ssgseaParam(as.matrix(expression_data_nor), gset_acc)
gset_stem_ssgsea <- gsva(ssgseaPar, verbose=FALSE)
gset_stem_ssgsea = as.data.frame(gset_stem_ssgsea)


##
transposed_data <- as.data.frame(t(gset_stem_ssgsea))
transposed_data <- data.frame(ACC = rownames(transposed_data), transposed_data)
transposed_data$EN_Acc_Sig <- scale(transposed_data$EN_Acc_Sig_Up) - scale(transposed_data$EN_Acc_Sig_Dn) 

## merge into the new table with metadata
merged_data_sur <- merge(age_results_sub, transposed_data, by.x = "SampleID", by.y = "ACC")

## define high vs low by 75% percentile
merged_data_sur$EN_Acc_group <- ifelse(merged_data_sur$EN_Acc_Sig >= quantile(merged_data_sur$EN_Acc_Sig,0.75), "High", "Low")
merged_data_sur$EN_Acc_group <- factor(merged_data_sur$EN_Acc_group, levels = c("Low", "High"))


# PFS
surv_object <- Surv(time = merged_data_sur$PFS_time, event = merged_data_sur$PFS_event)
fit <- survfit(surv_object ~ EN_Acc_group, data = merged_data_sur)

pdf("mAge-Acceleration.EN_Acc_Sig.q75.PFS.pdf",height=4.8,width=4.2)
plot.new()
p = ggsurvplot(fit, data = merged_data_sur, pval = TRUE, risk.table = TRUE,
           title = "PFS by EN_Acc_Sig",
           xlab = "Time (Years)", ylab = "PFS (%)",
           legend.title = "Group",
           legend.labs = c("Low","High"),
           legend = c(0.8, 0.6), font.legend=12,
           font.x = c(12),
           font.y = c(12),
           font.tickslab = c(12),
           risk.table.height = 0.18,
           risk.table.fontsize = 12/.pt,
           palette = c("#2E9FDF", "#FF2400"),
          tables.theme = theme_cleantable() + 
          theme(plot.title = element_text(size = 12))
        )
p$plot <- p$plot +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

print(p, newpage = FALSE)
dev.off()


##############################################################################
## Cell type deconvolution analysis using xCell
##############################################################################

## xCell analysis was implemented with the immunedeconv package
library(immunedeconv)

##### xCell analysis #####

## cell type scoring
result_xcell <- deconvolute(expression_data_nor, method = "xcell")

## reformat the results, and merge into metadata table
xcell_result_numeric <- result_xcell
xcell_result_numeric[ , -1] <- lapply(result_xcell[ , -1], as.numeric)

transposed_data <- as.data.frame(t(as.data.frame(result_xcell)))
colnames(transposed_data) = transposed_data[1,]
transposed_data = transposed_data[-1,]
rn <- rownames(transposed_data)

transposed_data <- as.data.frame(lapply(transposed_data, as.numeric))
rownames(transposed_data) <- rn

transposed_data <- data.frame(ACC = rownames(transposed_data), transposed_data)
immunedeconv_xcell <- merge(age_results_sub, transposed_data, by.x = "SampleID", by.y = "ACC")


##### t.test with results of xCell analysis #####

## immune.score
t.test(subset(immunedeconv_xcell, mAge_group == "Acc")$immune.score, subset(immunedeconv_xcell, mAge_group == "Dec")$immune.score)
# 1.042910  0.774302, p-value = 3.878e-10
t.test(subset(immunedeconv_xcell, mAge_group == "Syn")$immune.score, subset(immunedeconv_xcell, mAge_group == "Dec")$immune.score)
# 1.018421  0.774302, p-value = 6.744e-09
t.test(subset(immunedeconv_xcell, mAge_group == "Acc")$immune.score, subset(immunedeconv_xcell, mAge_group == "Syn")$immune.score)
# 1.042910  1.018421, p-value = 0.143

## stroma.score
t.test(subset(immunedeconv_xcell, mAge_group == "Acc")$stroma.score, subset(immunedeconv_xcell, mAge_group == "Dec")$stroma.score)
# 0.03245126 0.08778818, p-value = 1.553e-12
t.test(subset(immunedeconv_xcell, mAge_group == "Syn")$stroma.score, subset(immunedeconv_xcell, mAge_group == "Dec")$stroma.score)
# 0.06109018 0.08778818, p-value = 0.0004647
t.test(subset(immunedeconv_xcell, mAge_group == "Acc")$stroma.score, subset(immunedeconv_xcell, mAge_group == "Syn")$stroma.score)
# 0.03245126 0.06109018, p-value = 1.878e-12



##### visualize xCell result by mAge groups #####

library(patchwork)

df_long <- immunedeconv_xcell %>%
  dplyr::select("mAge_group","immune.score","stroma.score") %>%   # 这里加上 dplyr::
  tidyr::pivot_longer(
    cols = -mAge_group,
    names_to = "variable",
    values_to = "xCell score"
  )

df_long$variable = factor(df_long$variable, levels = c("immune.score","stroma.score"))


pB <- df_long %>% filter(variable=="immune.score") %>%
  ggplot(aes(mAge_group, `xCell score`, fill = mAge_group)) +
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = c(Acc = "#BF360C98",
                               Dec = "#4CAF5098",
                               Syn = "#0091EA98")) +
  theme_pubr(base_size = 12) + theme(legend.position="none", axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x="", y="xCell immune score") + ylim(0,1.5)

pC <- df_long %>% filter(variable=="stroma.score") %>%
  ggplot(aes(mAge_group, `xCell score`, fill = mAge_group)) +
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = c(Acc = "#BF360C98",
                               Dec = "#4CAF5098",
                               Syn = "#0091EA98")) +
  theme_pubr(base_size = 12) + theme(legend.position="none", axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x="", y="xCell stroma score") + ylim(0,0.3)

p = (pB | pC)

ggsave("mAge-Acceleration-group.xcell.overall.scores.pdf",p,width = 3.3,height = 2.6)
