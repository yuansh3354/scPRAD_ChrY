# Load the CellChat library
library(CellChat)

# Select cells with TIME_Type 'BM Comrade' and subset the original object
use.cells = sce@meta.data[which(sce$TIME_Type == 'BM Comrade'),] %>% rownames()
Bone = subset(sce, cells = use.cells)

# Select cells with TIME_Type 'PC Comrade' and subset the original object
use.cells = sce@meta.data[which(sce$TIME_Type == 'PC Comrade'),] %>% rownames()
Prostate = subset(sce, cells = use.cells)

# Remove the original object and perform garbage collection
rm(list = c('sce'))
gc(reset = TRUE)

# Find variable features and scale data for Prostate and Bone objects
Prostate = FindVariableFeatures(Prostate) %>% ScaleData()
Bone = FindVariableFeatures(Bone) %>% ScaleData()

# Create a list of CellChat objects for Prostate and Bone
cell.chat.list = list(
  Prostate = createCellChat(object = Prostate, meta = Prostate@meta.data, group.by = "SubtypeAnnotation", assay = "RNA"),
  Bone = createCellChat(object = Bone, meta = Bone@meta.data, group.by = "SubtypeAnnotation", assay = "RNA")
)

# Remove intermediate objects and perform garbage collection
rm(list = c('Bone', 'Prostate'))
gc(reset = TRUE)

# Print the levels of cell identities for each CellChat object in the list
lapply(cell.chat.list, function(x) { print(levels(x@idents)) })

# Use human CellChat database
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB)

# Set maximum allowed size for global variables in future operations
options(future.globals.maxSize = 500 * 1024^3)

# Process each CellChat object in the list with parallel computation
cell.chat.list = lapply(cell.chat.list, function(cellchat) {
  cellchat@DB <- CellChatDB.use
  future::plan("multisession", workers = 8) # Use parallel computation
  cellchat <- subsetData(cellchat) # Subset data for processing
  cellchat <- identifyOverExpressedGenes(cellchat) # Identify overexpressed genes
  cellchat <- identifyOverExpressedInteractions(cellchat) # Identify overexpressed interactions
})

# Increase the maximum allowed size for global variables
options(future.globals.maxSize = 5000 * 1024^3)

# Further process each CellChat object in the list
cell.chat.list = lapply(cell.chat.list, function(cellchat) {
  cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.05) # Compute communication probability
  cellchat <- filterCommunication(cellchat, min.cells = 10) # Filter communications
  cellchat <- computeCommunProbPathway(cellchat, thresh = 0.05) # Compute communication probability pathways
  cellchat <- aggregateNet(cellchat) # Aggregate networks
})

# Reset the maximum allowed size for global variables
options(future.globals.maxSize = 500 * 1024^3)

# Final processing for each CellChat object in the list
cell.chat.list = lapply(cell.chat.list, function(cellchat) {
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # Compute centrality
  cellchat <- computeNetSimilarity(cellchat, type = "functional") # Compute network similarity (functional)
  cellchat <- netEmbedding(cellchat, type = "functional", umap.method = 'uwot') # Network embedding (functional)
  cellchat <- netClustering(cellchat, type = "functional") # Network clustering (functional)
})

# Save the list of CellChat objects to an RDS file
saveRDS(cell.chat.list, 'result/CellChatList_for_all.RDS')

object.list = cell.chat.list
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = P.bm, title = names(object.list)[i], width = 15, height = 7)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = P.pc, title = names(object.list)[i+1], width = 15, height = 7)
draw(ht1 %v% ht2, ht_gap = unit(0.5, "cm"))
Ca.colors = c("ProlifCs"="#7E6148B2", "NEndoCs"="#DC0559b9",  "MetaCs"="#91D1C2B2" ,"InvCs" = "#9970acb2",   "ImmuCs" ="#FFDC91FF",  "EMTCs" ="#4DBBD5B2", "OstbCs" ="#5591B4B2", "OstcCs"="#3C5488B2")
x1 = netAnalysis_signalingChanges_scatter(cellchat, idents.use = names(Ca.colors)[1])
x2 = netAnalysis_signalingChanges_scatter(cellchat, idents.use = names(Ca.colors)[2])
x3 = netAnalysis_signalingChanges_scatter(cellchat, idents.use = names(Ca.colors)[3])
x4 = netAnalysis_signalingChanges_scatter(cellchat, idents.use = names(Ca.colors)[4])
x5 = netAnalysis_signalingChanges_scatter(cellchat, idents.use = names(Ca.colors)[5])
x6 = netAnalysis_signalingChanges_scatter(cellchat, idents.use = names(Ca.colors)[6])
x7 = netAnalysis_signalingChanges_scatter(cellchat, idents.use = names(Ca.colors)[7])
x8 = netAnalysis_signalingChanges_scatter(cellchat, idents.use = names(Ca.colors)[8])
x = (x1+x2+x3+x4)|(x5+x6+x7+x8)
print(x)

x = compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
print(x)
