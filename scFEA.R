library(ComplexHeatmap)
library(Seurat)
library(dplyr)
# Read the adj_flux.csv file
adj_flux <- read.csv("./scFEA/output/adj_flux.csv", row.names = 1)
rownames(adj_flux) <- gsub("\\.", "-", rownames(adj_flux))
ids = rownames(adj_flux)
# Fetch metadata from the Seurat object
cell_anno = FetchData(sce, vars = c('SubtypeAnnotation', 'SITE', 'SUBTYPE', 'DISEASE', 'TIME_Type'))
cell_anno = cell_anno[ids, , drop = FALSE]
cell_anno <- cell_anno[order(cell_anno$SubtypeAnnotation), , drop = FALSE]
ids = rownames(cell_anno)
adj_flux <- adj_flux[ids, ]
print(dim(adj_flux))
cell_anno$group = paste0(cell_anno$SUBTYPE, '(', cell_anno$TIME_Type, ')')
Attribute = 'SUBTYPE'

# Calculate average metabolic flux for each group
df_averages <- adj_flux %>%
  group_by(group = cell_anno[[Attribute]]) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  select(-group)
rownames(df_averages) <- unique(cell_anno[[Attribute]])
df_averages <- t(df_averages) %>% as.data.frame()

df_flux <- df_averages
# Remove rows with zero standard deviation
df_flux = df_flux[apply(df_flux, 1, function(x) sd(x) != 0), ]
human_moduleInfo <- read.csv("./scFEA.human.moduleinfo.csv", header = TRUE, row.names = 1)

# Annotate modules and replace flux row names
human_moduleInfo$module_name <- paste0(human_moduleInfo$Module_id, ": ",
                                       human_moduleInfo$Compound_IN_name,
                                       "_", human_moduleInfo$Compound_OUT_name)
select_moduleInfo = human_moduleInfo[rownames(df_flux), ]
df_flux_new <- df_flux
rownames(df_flux_new) <- select_moduleInfo$module_name

# Indicate key modules of interest
modules = human_moduleInfo[human_moduleInfo$Module_id %in% mysplit('M_49,M_34,M_35,M_105,M_78,M_65,M_66,M_167,M_48,M_51,M_92,M_169'), 'module_name'] %>% unique()
modules <- as.data.frame(modules)

# Initialize a list to store top 10 genes per column
top_10_genes_per_column <- list()
# Iterate over each column to find top 10 genes
for (col in colnames(df_flux_new)) {
  col_values <- df_flux_new[, col]
  top_genes <- rownames(df_flux_new)[order(col_values, decreasing = TRUE)[1:10]]
  top_10_genes_per_column[[col]] <- top_genes
}

ids = unique(c(unlist(top_10_genes_per_column), modules$modules))
ids = modules$modules
pheatmap(df_flux_new, scale = "row",
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_cols = TRUE, cluster_rows = TRUE,
         col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(512),
         fontsize_col = 8,
         treeheight_row = 20,
         heatmap_legend_param = list(title = "Flux")) + 
  rowAnnotation(link = anno_mark(at = which(rownames(df_flux_new) %in% ids), 
                                 labels = ids, labels_gp = gpar(fontsize = 8)))

# Subset Seurat object based on adj_flux cell IDs
adj_scRNA = subset(sce, cells = rownames(adj_flux))
predFlux <- t(data.matrix(adj_flux))
adj_scRNA[["FLUX"]] <- CreateAssayObject(counts = predFlux)
adj_scRNA$group = paste0(adj_scRNA$SUBTYPE, '(', adj_scRNA$SITE, ')')
DefaultAssay(adj_scRNA) <- 'FLUX'
human_moduleInfo$SM_anno %>% unique()
myModule = split(gsub('_', '-', human_moduleInfo$Module_id), human_moduleInfo$SM_anno)

adj_scRNA = AddModuleScore(adj_scRNA, features = myModule, ctrl = 5)
colnames(adj_scRNA@meta.data)[grep('Cluster', colnames(adj_scRNA@meta.data))] = names(myModule)
adj_scRNA$group = paste0(gsub(' Comrade', '', sce$TIME_Type), '(', adj_scRNA$SubtypeAnnotation, ')')

x = jjDotPlot(object = adj_scRNA,
              gene = names(myModule), 
              xtree = FALSE,
              ytree = TRUE, lwd = 0.2, bar.width = 3,
              id = 'group',
              dot.col = Strength.cls,
              point.geom = TRUE,
              rescale = TRUE, legend.position = 'bottom',
              tile.geom = TRUE) + 
  scale_size(limits = c(2, 100), range = c(2, 5))
x
# Differential analysis with volcano plot
Idents(adj_scRNA) = 'celltype'
cells <- mysplit('Epithelium')
DEG_FLUX <- list()
for (i in 1:length(cells)) {
  obj = subset(adj_scRNA, idents = cells[i])
  Idents(obj) = 'TIME_Type'
  module_deg = FindMarkers(obj, 
                           ident.1 = "BM Comrade", 
                           ident.2 = "PC Comrade",
                           min.cells.group = 1,
                           min.cells.feature = 1,
                           min.pct = 0,
                           logfc.threshold = 0)
  
  rownames(module_deg) <- gsub("-", "_", rownames(module_deg))
  module_deg$MN <- rownames(module_deg)
  module_deg <- module_deg[human_moduleInfo$Module_id, ]
  module_deg <- cbind(module_deg, human_moduleInfo)
  module_deg$M_name <- paste(module_deg$Compound_IN_name, "→", module_deg$Compound_OUT_name)
  DEG_FLUX[[i]] <- module_deg
  names(DEG_FLUX)[i] <- cells[i]
}

# Add cell type and group information to adj_flux
adj_flux$celltype <- adj_scRNA$celltype
adj_flux$group <- adj_scRNA$SITE

# Calculate Cohen's d for each flux
for (i in 1:length(DEG_FLUX)) {
  cells_flux <- adj_flux[adj_flux$celltype == names(DEG_FLUX)[i], ]
  cells_SD <- rownames(cells_flux)[cells_flux$group == "Bone"]
  cells_HC <- rownames(cells_flux)[cells_flux$group == "Prostate"]
  
  df <- as.data.frame(t(cells_flux))
  df <- df[-c(169, 170), ]
  
  df_SD <- df[, cells_SD]
  df_HC <- df[, cells_HC]
  
  for (flux_id in rownames(DEG_FLUX[[i]])) {
    A <- as.numeric(df_SD[flux_id, ])
    B <- as.numeric(df_HC[flux_id, ])
    c_d <- effsize::cohen.d(A, B)$estimate
    DEG_FLUX[[i]][flux_id, 'cohens_d'] <- c_d
  }
}

# Plot volcano plot
library(ggplot2)
library(ggrepel)
library(cowplot)

df = DEG_FLUX[['Epithelium']] %>% na.omit()
mymodule = df[order(abs(df$p_val_adj), decreasing = FALSE), ] %>% rownames()
ids = mysplit('M_49,M_34,M_35,M_105,M_78,M_65,M_66,M_167,M_48,M_51,M_92,M_169')
p = ggplot(df, aes(x = cohens_d, y = -log10(p_val))) +
  geom_hline(aes(yintercept = 1.3), color = "#999999", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = 0), color = "#999999", linetype = "dashed", size = 1) +
  geom_point(size = 0.75, color = "grey80") +
  geom_point(data = df[df$Module_id %in% ids, ], stroke = 0.5, size = 2, shape = 16, aes(color = ifelse(cohens_d > 0, "Bone", "Prostate"))) +
  labs(x = "Cohen's D", y = "-Log10(P-value)", title = "") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1, colour = "black"),
        axis.title = element_text(size = 14), axis.text = element_text(size = 8, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  theme(plot.margin = unit(c(0, 1, 2, 1), 'cm')) +
  geom_text_repel(data = df[df$Module_id %in% ids, ], 
                  aes(label = M_name), color = "black", size = 4, 
                  arrow = arrow(ends = "first", length = unit(0.01, "npc")), box.padding = 0.2,
                  point.padding = 0.3, segment.color = 'black', segment.size = 0.35, force = 1, max.iter = 3e3,
                  max.overlaps = Inf) +
  scale_color_manual(values = c("#B51F29", "#006699")) + 
  NoLegend()
p

# Use cowplot to combine plots and other elements
x = ggdraw(xlim = c(0, 1), ylim = c(0, 1.1)) + 
  draw_plot(p, x = 0, y = 0) +  
  draw_line(x = c(0.55, 0.75), 
            y = c(0.1, 0.1),
            lineend = "round",
            size = 1, col = "#B51F29",  
            arrow = arrow(angle = 15, 
                          length = unit(0.1, "inches"),
                          type = "closed")) +
  draw_line(x = c(0.25, 0.45), 
            y = c(0.1, 0.1),
            lineend = "round",
            size = 1, col = "#006699",  
            arrow = arrow(angle = 15, 
                          length = unit(0.1, "inches"),
                          type = "closed",
                          ends = "first")) +
  draw_text(text = "Activate in Bone", size = 12,
            x = 0.88, y = 0.1,
            color = "black", fontface = "italic") +
  draw_text(text = "Activate in Prostate", size = 12,
            x = 0.15, y = 0.1,
            color = "black", fontface = "italic")
x