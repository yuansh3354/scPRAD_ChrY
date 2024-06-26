# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
# Load RDS files
mutation_data = readRDS('mPC_scomatic.rds')

# Summarize mutation data by SampleID and Cell_types
result <- mutation_data %>%
  group_by(SampleID, Cell_types) %>%
  summarise(Count = log(n(), 2))
result
# Create a ggplot object for mutation burden
x = ggplot(result, aes(x = as.numeric(Cell_types), y = Count)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.1, cex = 1) +
  stat_summary(geom = 'point', fun = 'mean', aes(fill = Cell_types), size = 4, pch = 21) +
  stat_summary(fun = "mean", geom = "line", cex = 1) +
  yuansh_theme +
  scale_x_continuous(breaks = 1:8, labels = mysplit('MetaCs,OstcCs,ProlifCs,EMTCs,OstbCs,NEndoCs,ImmuCs,InvCs')) +
  labs(x = "", y = "") +
  ggtitle("Log2(Mut burden)") + 
  scale_fill_manual(values = Ca.colors)
x

# Merge data with Seurat metadata and write to Excel
df = x$data %>% as.data.frame() %>% na.omit()
df = merge(df, unique(sce@meta.data[, mysplit('SampleID,ID')]))

# Filter mutation data for chrY and summarize by Cell_types and Gene
chrY = mutation_data[mutation_data$CHROM == 'chrY',]
Gene = unique(chrY$Gene)
Gene <- Gene[!grepl("-", Gene)]
result <- mutation_data %>%
  group_by(Cell_types, Gene) %>%
  summarise(Count = log(n(), 2))
result = result[which(result$Gene %in% Gene),]
genes = names(which(table(result$Gene) > 7))
result = result[which(result$Gene %in% genes),]

# Summarize mutation data by Cell_types, Gene, and CHROM
result <- mutation_data %>%
  group_by(Cell_types, Gene, CHROM) %>%
  summarise(Count = log(n(), 2))
genes = names(which(table(result$Gene) > 7))
result = result[which(result$Gene %in% genes),]
result = result[order(result$Count, decreasing = TRUE),]
result <- result[!grepl("-", result$Gene),]
Genes = unique(result$Gene)[1:10]
result = result[which(result$Gene %in% Genes),]
result = na.omit(result)

# Create a ggplot object for gene mutation burden
x = ggplot(result, aes(x = Cell_types, y = Count, group = Gene, color = Gene)) +
  geom_line() +
  geom_point() +
  labs(x = "", y = "", title = "") +
  theme_minimal() +
  scale_color_simpsons()
x


# Filter mutation data for chrY and summarize by Gene
chrY = mutation_data[mutation_data$CHROM == 'chrY',]
Gene = unique(chrY$Gene)
Gene <- Gene[!grepl("-", Gene)]
result <- mutation_data %>%
  group_by(Gene) %>%
  summarise(Count = log(n(), 2))
result = result[which(result$Gene %in% Gene),]
result = result[order(result$Count, decreasing = TRUE),]
result = result[1:10,]

# Create a ggplot object for chrY gene mutation burden
x = ggplot(result, aes(reorder(Gene, Count, decreasing = TRUE), Count)) + 
  geom_bar(stat = 'identity', fill = '#00A087B2') + 
  yuansh_theme
x

# Load necessary libraries
library(maftools)
library(readxl)
library(ggplot2)

# Read MAF and clinical data
laml <- read.maf(maf = "maf_data.txt", clinicalData = "clinical_data.txt", vc_nonSyn = NULL)

# Perform TiTv analysis and plot results
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)

# Convert MAF data to data frame
df = as.data.frame(laml@data)

# Update SampleID in Seurat object
sce$SampleID = paste0(sce$SubtypeAnnotation, "(", sce$SampleID, ")")

# Extract genes of interest from mutation data
genes_of_interest <- unique(mutation_data[which(mutation_data$CHROM == 'chrY'), 'Gene'])

# Read TCGA MAF and clinical data
laml.tcga = read.maf('mc3.v0.2.8.PUBLIC.maf.gz')
phe = readxl::read_excel('TCGA-CDR-SupplementalTableS1.xlsx', sheet = 1) %>% as.data.frame()
phe = phe[, -1]  # Remove first column
colnames(phe)[1] = 'Tumor_Sample_Barcode'
df = as.data.frame(laml.tcga@data)
df$Tumor_Sample_Barcode = substring(df$Tumor_Sample_Barcode, 1, 12)
cg = as.character(phe$Tumor_Sample_Barcode)[phe$gender == 'MALE']
laml <- read.maf(df[df$Tumor_Sample_Barcode %in% cg, ], clinicalData = phe[phe$gender == 'MALE', ], vc_nonSyn = NULL)

# Define pathways and genes
pathways = data.frame(
  Genes = c(
    "GNMT", "ACSM1", "MMP26", "PCAT18", "TFF3", "EPHA6", "FABP5", "PCA3", "NUDT8", "ALOX15B", 
    "OR51E2", "AMACR", "CAMKK2", "TRPM4", "TMEFF2", "SPON2", "HPN", "ANPEP", "HSD17B6", "SMS", 
    "FAM13C", "MYBPC1", "SMIM22", "REPS2", "MCCC2", "LUZP2", "FASN", "PCAT4", "GOLM1", "TMSB15A", 
    "SDK1", "MESP1", "PDLIM5", "ALDH1A3", "H2AFJ", "MIPEP", "TRPM8", "PNKP", "SPOCK1", "SLC2A12", 
    "TRGC1", "LRRC26", "RAMP1", "RAB3B", "PPP1R14B", "DUS1L", "ABCC4", "CORO1B", "TSPAN1", "ENTPD5",
    mysplit('KLK2,KLK3,KLK4,APOO,DHRS7,AZGP1,BMPR1B,FOXA1,PTOV1'),
    mysplit('NPAS2,ADORA2A,GABRA2,AVP,PER1, PER2, PER3,SIK3,CRY1,NFAT5'),
    "SRY", "TSPY1", "USP9Y", "UTY", "ZFY", "KDM5D", "NLGN4Y"
  ),
  Pathway = rep(c(
    "Metabolism", "Signal Transduction", "Extracellular Matrix", "Transcription Regulation", 
    "Transport", "Prostate Cancer", "Circadian rhythm", "Y Chromosome"
  ), c(6, 15, 10, 9, 10, 9, 10, 7)),
  stringsAsFactors = FALSE
)

# Define variant classification colors
vc_cols = c(Multi_Hit = '#DC9F64', Missense_Mutation = '#DE5613', Nonsense_Mutation = '#4B4F5B')

# Generate oncoplot and save as PDF
pdf(file = 'plts/SubtypeAnnotation-Mut.pdf', height = 12, width = 12)
oncoplot(
  maf = laml, pathways = pathways,
  clinicalFeatures = c('SubtypeAnnotation'),
  sortByAnnotation = TRUE,
  annotationColor = list(SubtypeAnnotation = Ca.colors, SITE = Site.colors),
  colors = vc_cols,
  bgCol = "#F0F0F0"
)
dev.off()

# Subset MAF data by genes in pathways
maf_subset <- subsetMaf(laml, genes = pathways$Genes)
sce$Group <- ifelse(sce$SampleID %in% maf_subset@clinical.data$Tumor_Sample_Barcode, "mutated", "non-mutated")
sce$Group = paste0(sce$Group, '(', sce$SubtypeAnnotation, ')')

# Generate dot plot for pathways genes
gg = jjDotPlot(object = sce,
               gene = pathways$Genes,
               xtree = TRUE,
               ytree = FALSE, lwd = 0.2, bar.width = 3,
               id = 'Group', rescale = TRUE,
               point.geom = FALSE,
               dot.col = Strength.cls,
               legend.position = 'bottom',
               tile.geom = TRUE) + 
  scale_size(limits = c(0, 100), range = c(1, 5))
gg

# Save survival plot as PDF
pdf(file = 'plts/TCGA-MALE-Mut-Survival.pdf', width = 6, height = 4)
mafSurvival(maf = laml, genes = pathways$Genes, time = "PFI.time", Status = "PFI", isTCGA = TRUE)
dev.off()

# Additional survival analysis
mafSurvival(maf = laml, genes = mysplit('TSPAN8,ALDOA,HSP90AA1,SDHD,CLCN3,SRSF5,PTMA'), time = "PFI.time", Status = "PFI", isTCGA = TRUE)

# Subset Seurat object by unique sample IDs
maf_subset <- subsetMaf(laml, genes = pathways$Genes)
Idents(sce) = sce$SampleID 
sce = subset(sce, ident = unique(laml@clinical.data$Tumor_Sample_Barcode))
sce$Group <- ifelse(sce$SampleID %in% maf_subset@clinical.data$Tumor_Sample_Barcode, "mutated", "non-mutated")
table(sce$Group)

# Dimensional plot for groups
DimPlot(sce, group.by = 'Group')

# Differential expression analysis
Idents(sce) = 'Group'
degs = FindMarkers(sce, ident.1 = 'mutated', ident.2 = 'non-mutated', min.pct = 0.25)

# Visualize cell type distribution
x = visualizeCellTypeDistribution(sce@meta.data, 'SubtypeAnnotation', 'Group') + 
  scale_fill_manual(values = alpha(c("maroon", "royalblue"), 0.8))
x
