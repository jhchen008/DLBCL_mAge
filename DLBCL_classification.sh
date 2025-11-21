
#########################################################################
## COO classification by ABC/GCB gene sets
#########################################################################
## R codes

## define the scoring function
coo_by_mean_zscore <- function(expr,
                               abc_genes,
                               gcb_genes,
                               diff_thr = 0.25,
                               mean_cap = 0.75) {
  stopifnot(is.matrix(expr) || is.data.frame(expr))
  expr <- as.matrix(expr)
  
  # convert to upper cases
  rownames(expr) <- toupper(rownames(expr))
  abc_genes <- toupper(abc_genes)
  gcb_genes <- toupper(gcb_genes)
  
  abc_use <- intersect(abc_genes, rownames(expr))
  gcb_use <- intersect(gcb_genes, rownames(expr))
  if (length(abc_use) < 3) warning("Too few ABC genes：", length(abc_use))
  if (length(gcb_use) < 3) warning("Too few GCB genes：", length(gcb_use))
  
  # data scaling
  z_row <- function(v){
    s <- stats::sd(v, na.rm = TRUE); m <- mean(v, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(v)))  # 方差为0则设为0
    (v - m) / s
  }
  zexpr <- t(apply(expr, 1, z_row))
  
  # calculate the mean values
  ABC_mean <- colMeans(zexpr[abc_use, , drop = FALSE], na.rm = TRUE)
  GCB_mean <- colMeans(zexpr[gcb_use, , drop = FALSE], na.rm = TRUE)
  diff_val <- ABC_mean - GCB_mean
  
  # classification
  call <- ifelse(diff_val >=  diff_thr & GCB_mean <= mean_cap, "ABC",
          ifelse(diff_val <= -diff_thr & ABC_mean <= mean_cap, "GCB", "Unclassified"))
  
  # result table
  out <- data.frame(
    Sample = colnames(expr),
    ABC_mean = ABC_mean,
    GCB_mean = GCB_mean,
    Diff = diff_val,
    COO = call,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  attr(out, "used_genes") <- list(ABC = abc_use, GCB = gcb_use)
  attr(out, "params") <- list(diff_thr = diff_thr, mean_cap = mean_cap)
  return(out)
}


## ABC and GCB gene sets from Reddy et al. PMID: 28985567 
ABC_set <- c("PIM1","ENTPD1","BLNK","CCND2","ETV6","FUT8","BMF","IL16","PTPN1","SH3BP5","IRF4")
GCB_set <- c("ITPKB","MME","BCL6","MYBL1","DENND3","NEK6","LMO2","LRMP","SERPINA9")

## applying the 
res_coo <- coo_by_mean_zscore(expression_data_nor, ABC_set, GCB_set, diff_thr = 0.25, mean_cap = 0.75)
res_coo$COO_est = res_coo$COO



#########################################################################
## LME classification
#########################################################################
## Python codes

## https://github.com/BostonGene/LME/blob/main/LME_Classification.ipynb

##### load modules #####
%load_ext autoreload
%autoreload 2
%matplotlib inline
%config IPCompleter.use_jedi = False
import matplotlib.pyplot as plt
import os
import pandas as pd
import numpy as np
import seaborn as sns
from scipy import stats

%config InlineBackend.figure_format = 'png'
plt.rcParams['pdf.fonttype'] = 'truetype'
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['figure.dpi'] = 120
sns.set_style('white')

import warnings
warnings.filterwarnings("ignore")

from lme.utils import *
from lme.plotting import *
from lme.classification import *
from lme.pathway_scoring import *



##### loading reference for LME classification #####
## this should be found within the LME folder
REFERENCE_COHORT_ANNOTATION = './datasets/pan_cohort_annotation.tsv'
REFERENCE_COHORT_EXPRESSION = './datasets/pan_cohort_signatures.tsv.gz'
GENE_SIGNATURES = 'databases/signatures.gmt'

## list of selected/needed signatures for further calssification.
signatures_selected = [
 'Lymphatic_endothelium',
 'Angiogenesis',
 'CAF',
 'Fibroblastic_reticular_cells',
 'Matrix',
 'Matrix_remodeling',
 'Granulocyte_traffic',
 'Protumor_cytokines',
 'Follicular_dendritic_cells',
 'Macrophages',
 'M1_signature',
 'T_cell_traffic',
 'MHCII',
 'MHCI',
 'Follicular_B_helper_T_cells',
 'Treg',
 'T_cells',
 'Checkpoint_inhibition',
 'NK_cells',
 'B_cells_traffic',
 'B_cells',
 'Proliferation_rate']
len(signatures_selected)



##### Data preparation, reference cohort #####

# Load the reference cohort
cohort_signatures_scaled = read_dataset(REFERENCE_COHORT_EXPRESSION).T
cohort_signatures_scaled.shape
cohort_signatures_scaled.head()

# Read reference cohort annotation file
cohort_ann = read_dataset(REFERENCE_COHORT_ANNOTATION)
cohort_ann.shape
cohort_ann.head()

len(cohort_signatures_scaled.index & cohort_ann.index)

# Load gene signatures
signatures = read_gene_sets(GENE_SIGNATURES)

# Create a KNN model based on reference cohort
cohort_ann_filtered = cohort_ann[(cohort_ann.Diagnosis=='Diffuse_Large_B_Cell_Lymphoma') & 
                                 (~cohort_ann.LME.isna())]
cohort_ann_filtered.shape

LME_MODEL = KNeighborsClusterClassifier(norm=False, scale=False, clip=3, k=35).fit(
    *to_common_samples([cohort_signatures_scaled[signatures_selected + progeny_selected], cohort_ann_filtered.LME]))



##################################################
##### DLBCL dataset: Schmitz et al NEJM 2018 #####

SAMPLE_EXPRESSION = './datasets/Schmitz_2018NEJM/RNAseq_gene_expression_562.formated.txt'

expression_matrix = read_dataset(SAMPLE_EXPRESSION)
expression_matrix_t = expression_matrix.T

# Check if expression matrix is normalized if not log2 transform it.
if all(0<=sample<=18 for sample in expression_matrix_t.mean()):
    print(expression_matrix_t.head())
else:
    expression_matrix_t = np.log2(1+expression_matrix_t)
    print(expression_matrix_t.head())

# Calculate ssGSEA and PROGENy scores of signatures
ssgsea_scores = ssgsea_formula(expression_matrix_t, signatures)
progeny_scores = run_progeny(expression_matrix_t)

signatures_calculated = pd.concat([ssgsea_scores, progeny_scores], axis=1)
signatures_calculated.shape
signatures_calculated.head()

# median scaling
signatures_scaled = median_scale_genes(signatures_calculated)
signatures_scaled.shape
signatures_scaled.head()

lme_predicted = LME_MODEL.predict(signatures_scaled[LME_MODEL.X.columns]).rename('LME')
lme_predicted.to_csv('Schmitz_2018NEJM/lme_predicted.sample.csv', index=True)



##################################################
##### DLBCL dataset: Lacy et al Blood 2020 #####

SAMPLE_EXPRESSION = './datasets/Lacy_2020Blood/GSE181063_series_matrix.nor.DLBCL.txt'

expression_matrix = read_dataset(SAMPLE_EXPRESSION)
expression_matrix_t = expression_matrix.T

# Check if expression matrix is normalized if not log2 transform it.
if all(0<=sample<=18 for sample in expression_matrix_t.mean()):
    print(expression_matrix_t.head())
else:
    expression_matrix_t = np.log2(1+expression_matrix_t)
    print(expression_matrix_t.head())

# Calculate ssGSEA and PROGENy scores of signatures
ssgsea_scores = ssgsea_formula(expression_matrix_t, signatures)
progeny_scores = run_progeny(expression_matrix_t)

signatures_calculated = pd.concat([ssgsea_scores, progeny_scores], axis=1)
signatures_calculated.shape
signatures_calculated.head()

# median scaling
signatures_scaled = median_scale_genes(signatures_calculated)
signatures_scaled.shape
signatures_scaled.head()

lme_predicted = LME_MODEL.predict(signatures_scaled[LME_MODEL.X.columns]).rename('LME')
lme_predicted.to_csv('Lacy_2020Blood/lme_predicted.sample.csv', index=True)



##################################################
##### DLBCL dataset: Sha et al JCO 2018 #####

SAMPLE_EXPRESSION = './datasets/Sha_2018JCO/GSE117556_series_matrix.nor.symbol.txt'

expression_matrix = read_dataset(SAMPLE_EXPRESSION)
expression_matrix_t = expression_matrix.T

# Check if expression matrix is normalized if not log2 transform it.
if all(0<=sample<=18 for sample in expression_matrix_t.mean()):
    print(expression_matrix_t.head())
else:
    expression_matrix_t = np.log2(1+expression_matrix_t)
    print(expression_matrix_t.head())

# Calculate ssGSEA and PROGENy scores of signatures
ssgsea_scores = ssgsea_formula(expression_matrix_t, signatures)
progeny_scores = run_progeny(expression_matrix_t)

signatures_calculated = pd.concat([ssgsea_scores, progeny_scores], axis=1)
signatures_calculated.shape
signatures_calculated.head()

# median scaling
signatures_scaled = median_scale_genes(signatures_calculated)
signatures_scaled.shape
signatures_scaled.head()

lme_predicted = LME_MODEL.predict(signatures_scaled[LME_MODEL.X.columns]).rename('LME')
lme_predicted.to_csv('Sha_2018JCO/lme_predicted.sample.csv', index=True)



##################################################
##### DLBCL dataset: Chapuy et al NatMed 2018 #####

SAMPLE_EXPRESSION = './datasets/Chapuy_2018NatMed/GSE98588_series_matrix.nor.symbol.txt'

expression_matrix = read_dataset(SAMPLE_EXPRESSION)
expression_matrix_t = expression_matrix.T

# Check if expression matrix is normalized if not log2 transform it.
if all(0<=sample<=18 for sample in expression_matrix_t.mean()):
    print(expression_matrix_t.head())
else:
    expression_matrix_t = np.log2(1+expression_matrix_t)
    print(expression_matrix_t.head())

# Calculate ssGSEA and PROGENy scores of signatures
ssgsea_scores = ssgsea_formula(expression_matrix_t, signatures)
progeny_scores = run_progeny(expression_matrix_t)

signatures_calculated = pd.concat([ssgsea_scores, progeny_scores], axis=1)
signatures_calculated.shape
signatures_calculated.head()

# median scaling
signatures_scaled = median_scale_genes(signatures_calculated)
signatures_scaled.shape
signatures_scaled.head()

lme_predicted = LME_MODEL.predict(signatures_scaled[LME_MODEL.X.columns]).rename('LME')
lme_predicted.to_csv('Chapuy_2018NatMed/lme_predicted.sample.csv', index=True)



#########################################################################
## Ecotype classification
#########################################################################
## shell codes with EcoTyper Rscript obtained from: https://github.com/digitalcytometry/ecotyper

Rscript EcoTyper_recovery_bulk.R -d Lymphoma -m ./datasets/Schmitz_2018NEJM/RNAseq_gene_expression_562.formated.txt -o ./datasets/Schmitz_2018NEJM/ecotyper

Rscript EcoTyper_recovery_bulk.R -d Lymphoma -m ./datasets/Lacy_2020Blood/GSE181063_series_matrix.nor.DLBCL.txt -o ./datasets/Lacy_2020Blood/ecotyper

Rscript EcoTyper_recovery_bulk.R -d Lymphoma -m ./datasets/Sha_2018JCO/GSE117556_series_matrix.nor.symbol.txt -o ./datasets/Sha_2018JCO/ecotyper

Rscript EcoTyper_recovery_bulk.R -d Lymphoma -m ./datasets/Chapuy_2018NatMed/GSE98588_series_matrix.nor.symbol.txt -o ./datasets/Chapuy_2018NatMed/ecotyper
