#Clear the environment
rm(list=ls())
gc()

# Install the R package on bioconductor
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
if(!require("BiocManager")) install.packages("BiocManager")
if(!require("TCGAbiolinks")) BiocManager::install("TCGAbiolinks")
if(!require("SummarizedExperiment")) BiocManager::install("SummarizedExperiment")
if(!require("DESeq2")) BiocManager::install("DESeq2")
if(!require("edgeR")) BiocManager::install("edgeR")
if(!require("limma")) BiocManager::install("limma")

# Install the R package above cran
if(!require("survival")) install.packages("survival")
if(!require("broom")) install.packages("broom")
if(!require("devtools")) install.packages("devtools")
if(!require("cli")) install.packages("cli")

devtools::install_github("ayueme/easyTCGA")

library(easyTCGA)
help(package="easyTCGA")


# TCGA pan-cancer，Download TCGA expression matrix and sample information：
# https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_gene_tpm.gz
# https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp

# The pan-cancer data of GTEx is also organized in the same way. First, the expression matrix file and sample information file are downloaded：
# https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/gtex_RSEM_gene_tpm.gz
# https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/GTEX_phenotype.gz

# File import
tcga_expr_file <- "tcga_RSEM_gene_tpm.gz"
tcga_clin_file <- "Survival_SupplementalTable_S1_20171025_xena_sp"
gtex_expr_file <- "gtex_RSEM_gene_tpm.gz"
gtex_pheno_file <- "GTEX_phenotype.gz"

# Generate the RDS file
getpancancer_xena(tcga_expr_file = tcga_expr_file,
                  tcga_clin_file = tcga_clin_file,
                  gtex_expr_file = gtex_expr_file,
                  gtex_pheno_file = gtex_pheno_file,
                  type = "tcga+gtex")



# Then wait, the RDS file will generate 10 rdata files, stored in the output_pancancer_xena folder in the current working directory:

# TCGA_pancancer_expr.rdata: TCGA expression matrix, rows are genes, columns are samples
# TCGA_pancancer_clin.rdata: The clinical information, sample quantity and sequence of TCGA samples are exactly the same as the sample quantity and sequence of the expression matrix above
# TCGA_pancancer_lncrna_clin.rdata: data that is integrated with lncRNA and sample information. Rows are samples, columns are genes, and the first 34 columns are clinical information, including survival data
# TCGA_pancancer_mrna_clin.rdata: data that integrates mRNA and sample information together. Rows are samples, columns are genes, and the first 34 columns are clinical information, including survival data.
# GTEx_pancancer_expr.rdata: GTEx expression matrix, rows are genes, columns are samples
# GTEx_pancancer_pheno.rdata: The sample information, sample quantity and order of GTEx are exactly the same as the sample quantity and order of the expression matrix above
# GTEx_pancancer_lncrna_pheno.rdata: data that is integrated with lncRNA and sample information. The row is the sample, the column is the gene, and the first two columns are sample_id and sample_type
# GTEx_pancancer_mrna_pheno.rdata: data that integrates mRNA and sample information. The row is the sample, the column is the gene, and the first two columns are sample_id and sample_type
# TCGA_GTEx_pancancer_lncRNA_pheno.rdata: lncRNA expression matrix and sample information integrated by TCGA and GTEx. Note: rows are samples! The first four columns are sample information, and the last column is lncRNA
# TCGA_GTEx_pancancer_mRNA_pheno.rdata: The mRNA expression matrix and sample information integrated by TCGA and GTEx. Note: The rows are samples! The first four columns are the sample information, and the following columns are the mRNA

