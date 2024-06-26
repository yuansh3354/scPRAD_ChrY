# Remove all objects from the environment
rm(list = ls())
# Perform garbage collection to free up memory
gc(reset = TRUE)

# Define the file path for output files
file.path = 'result/inferCNV/'

# Define paths for expression data, group information, and gene information files
expFile = paste0(file.path, 'expFile.txt')
groupFiles = paste0(file.path, 'groupFiles.txt')
geneFile = paste0(file.path, 'geneFile.txt')

# Convert RNA counts data to a data frame
dfcount = as.data.frame(sce@assays$RNA$counts)
# Create a data frame for group information with cell IDs
groupinfo = data.frame(cellId = colnames(dfcount))
# Add cell type information to the group info data frame
groupinfo$cellType = sce$celltype
# Get gene annotation information, sorted by chromosome and start position, and remove duplicate genes
geneInfor = annoGene(rownames(dfcount), "SYMBOL", 'human')
geneInfor = geneInfor[with(geneInfor, order(chr, start)), c(1, 4:6)]
geneInfor = geneInfor[!duplicated(geneInfor[, 1]), ]
# Filter the count data to include only genes present in the gene information
dfcount = dfcount[rownames(dfcount) %in% geneInfor[, 1], ]
# Reorder the count data to match the gene information order
dfcount = dfcount[match(geneInfor[, 1], rownames(dfcount)), ]

# Load the data.table library
library(data.table)
# Write the count data, group info, and gene info to their respective files
write(dfcount, file = expFile, sep = '\t', row.names = TRUE, quote = FALSE)
write(groupinfo, file = groupFiles, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
write(geneInfor, file = geneFile, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)

# Create an InferCNV object using the raw counts, annotations, and gene order files
my.infercnv_obj = CreateInfercnvObject(
  raw_counts_matrix = expFile,
  annotations_file = groupFiles,
  delim = "\t",
  gene_order_file = geneFile,
  ref_group_names = c("Normal")
)
# Save the InferCNV object to an RDS file
saveRDS(my.infercnv_obj, file = 'result/inferCNV/inferCNV.rds')

# Set the maximum allowed size for global variables in future operations
options(future.globals.maxSize = 500 * 1024^3)
# Remove all objects from the environment again
rm(list = ls())
# Perform garbage collection to free up memory again
gc(reset = TRUE)

# Run the InferCNV analysis using the previously saved InferCNV object
my.infercnv_obj = infercnv::run(
  readRDS('result/inferCNV/inferCNV.rds'),
  cutoff = 0.1,  # Use 1 for smart-seq, 0.1 for 10x-genomics
  out_dir = 'result/inferCNV/',  # Directory is auto-created for storing outputs
  num_threads = 32,
  HMM = FALSE,
  denoise = TRUE
)


# Load the InferCNV object from the RDS file
infercnv_obj = readRDS('result/inferCNV/run.final.infercnv_obj')
# Load the ComplexHeatmap library
library(ComplexHeatmap)

# Extract expression data from the InferCNV object
expr <- infercnv_obj@expr.data
# Get gene position information, annotated by symbol and human genome
gene_pos = annoGene(rownames(expr), "SYMBOL", 'human')
# Split genes by chromosome
chrCNV = split(gene_pos$SYMBOL, gene_pos$chr)
# Calculate the column means of (expression - 1)^2 for each chromosome
chrCNV = lapply(chrCNV, function(x){
  df = as.data.frame(colMeans((expr[x, ] - 1)^2))
})
# Assign column names to the chromosome data frames
for (i in names(chrCNV)) { colnames(chrCNV[[i]]) = i }
# Combine all chromosome data frames into one
chrCNV1 = do.call(cbind, chrCNV)

# Calculate the CNV score as the column means of (expression - 1)^2
CNV_score = as.data.frame(colMeans((expr - 1)^2))
colnames(CNV_score) = "CNV_score"
CNV_score$cell_id = rownames(CNV_score)

# Get metadata and add cell IDs
meta = sce@meta.data
meta$cell_id = rownames(meta)
# Merge CNV score with metadata
CNV_score = merge(meta, CNV_score)

# Extract reference and observation data from the InferCNV object
ref = infercnv_obj@expr.data[, infercnv_obj@reference_grouped_cell_indices$Normal]
obs = infercnv_obj@expr.data[, infercnv_obj@observation_grouped_cell_indices$Epithelium]

# Define a function to estimate CNV based on observations and references
estimateCNV <- function(obs, ref, score_threshold, cor_threshold) {
  cnv_obs <- colMeans((obs - 1)^2)
  cnv_ref <- colMeans((ref - 1)^2)
  cell_top <- names(sort(cnv_obs, decreasing = TRUE))[1:round(length(cnv_obs) * 0.5)]
  cnv_top <- rowMeans(obs[, cell_top])
  cor_obs <- apply(obs, 2, function(x) cor(x, cnv_top))
  cor_ref <- apply(ref, 2, function(x) cor(x, cnv_top))
  
  cnv <- data.frame(score = c(cnv_obs, cnv_ref), 
                    cor = c(cor_obs, cor_ref), 
                    barcode = c(colnames(obs), colnames(ref)))
  
  cnv$type <- 'Other'
  cnv$type[cnv$score > score_threshold & cnv$cor > cor_threshold] <- 'Malignant'
  cnv$type[cnv$score < score_threshold & cnv$cor < cor_threshold] <- 'Not Malignant'
  return(cnv)
}

# Estimate CNV scores and correlations with specified thresholds
CNV_score_cor <- estimateCNV(obs, ref, 
                             score_threshold = 0.001, 
                             cor_threshold = 0.25)
# Rename the third column to "cell_id"
colnames(CNV_score_cor)[3] <- "cell_id"

# Merge CNV score correlations with metadata
CNV_score_cor = merge(CNV_score_cor, meta)

# Load the ggplot2 library for plotting
library(ggplot2)

# Create a ggplot object for CNV score and correlation 
gg = ggplot(CNV_score_cor, 
            aes(x = score, y = cor, color = SubtypeAnnotation)) +
  geom_point() +  # Add points to the plot
  theme_bw() +  # Apply a white background theme
  facet_wrap(~ID, ncol = 10) +  # Create facets based on ID, with 10 columns
  theme(axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.title = element_text(color = 'black', size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "CNV Score",
       y = "CNV Correlation") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +  # Add a dashed horizontal line at y = 0.5
  scale_color_manual(values = Ca.colors) +  # Set manual color scale
  yuansh_theme  # Apply custom theme

# Print the ggplot object
gg

# Save the plot as a PNG file
ggsave('plts/CNSscoreforID.png', width = 30, height = 20, dpi = 300)

# Create a ggplot object for violin plot of CNV scores across different cell types
gg = ggplot(CNV_score_cor, aes(celltype, score)) +
  geom_violin(aes(fill = celltype), color = "black") +  # Add violin plot with black border
  theme_bw() +  # Apply a white background theme
  theme(axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title = element_text(color = 'black', size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1),
               geom = "pointrange", color = "black", size = 0.5) +  # Add summary statistics
  scale_fill_manual(values = c(Ca.colors, Normal = '#BDD8C5'))  # Set manual fill color scale

# Print the ggplot object
gg
