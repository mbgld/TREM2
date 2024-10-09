### 0. Clear R environments and get libraries
outdir = './Results/CellChat/'
library(CellChat); library(patchwork); library(Seurat)
library(future); options(stringsAsFactors = FALSE); future::plan("multicore", workers = 4)

###############################################
##### Part A. Preprocessing dataset and DB ####
###############################################

### 1. Get dataset
MY.sfj <- readRDS('./Results/MY.sfj')
MY.sfj <- SetIdent(MY.sfj, value = 'subcell.id')
TREM2.cell <- WhichCells(MY.sfj, expression = TREM2 > 2.5)
Idents(object=MY.sfj, cells = TREM2.cell) <- "TREM2(+):Mono-mc"
MY.sfj$subcell.id <- Idents(MY.sfj)

### 2. Prep metadata
meta = MY.sfj@meta.data
meta$subcell.id = droplevels(meta$subcell.id, exclude = setdiff(levels(meta$subcell.id),unique(meta$subcell.id)))

### 3. Create Cellchat object
cellchat.obj <- createCellChat(MY.sfj, meta = meta, group.by = 'subcell.id')
cellchat.obj@idents <- factor(cellchat.obj@idents, 
                              levels=c("TREM2(+):Mono-mc", 
                                       "Mo-lineage", "Mono-Mc","CD163/LGMN", "FCN1-mono",
                                       "mo-DC","cDC1","cDC2","pDC", 
                                       "Alv-Mc1", "Alv-Mc2", "Alv-Mc3", "Alv-Mc4",
                                       "Alv-Mc5", "Alv-Mc6", "Alv-Mc7", "Alv-Mc8",
                                       "Prol-Mc"))
### 4. Set database
# 4.1 Check CellChatDB
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

# 4.2 Select a subset of cellchatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB

# 4.3 Set the selected DB in the object
cellchat.obj@DB <- CellChatDB.use

### 5. Preprocessing expression dataset
cellchat.obj <- subsetData(cellchat.obj); future::plan("multicore", workers = 4)
cellchat.obj <- identifyOverExpressedGenes(cellchat.obj)
cellchat.obj <- identifyOverExpressedInteractions(cellchat.obj)

###############################################################
#### Part B. Inference of cell-cell communications network ####
###############################################################

### 1. Compute communication prob & infer cellular communication network
computeAveExpr(cellchat.obj, features = 'TREM2', type='truncatedMean', trim=0.1)
cellchat.obj <- computeCommunProb(cellchat.obj)
cellchat.obj <- filterCommunication(cellchat.obj, min.cells = 10)

### 2. Extract inferred cellular communication network as a df
df.net <- subsetCommunication(cellchat.obj)

### 3. Infer the cell-cell communication at a signaling pathway level
cellchat.obj <- computeCommunProbPathway(cellchat.obj)

### 4. Calculate aggregated cell-cell communication network
## 4.1 Select cell cluster
# 4.1.1 Using all cell cluster
cellchat.obj <- aggregateNet(cellchat.obj)

# 4.1.2 Using TREM2(+):Mono-mc as key cluster
source.cluster = c("TREM2(+):Mono-mc")
target.cluster = c("CD163/LGMN", "Mo-lineage", "mo-DC", "Alv-Mc2", "Mono-Mc", "cDC1", "Alv-Mc7", "Alv-Mc3", "FCN1-mono", "Alv-Mc5", "Alv-Mc1", "pDC", "Alv-Mc4", "cDC2", "Prol-Mc", "Alv-Mc6")
cellchat.obj <- aggregateNet(cellchat.obj,  sources.use = source.cluster, targets.use = target.cluster, remove.isolate = FALSE)

## 4.2 Visualize
# 4.2.1
groupSize <- as.numeric(table(cellchat.obj@idents))

# 4.2.1.1 First graph
svg(paste0(outdir,"No_interaction_TREM2_netVisual.svg"), width = 16, height = 10, pointsize = 16)
netVisual_circle(cellchat.obj@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
dev.off()

# 4.2.1.2 Second graph
svg(paste0(outdir,"Weight_interaction_TREM2_netVisual.svg"), width = 16, height = 10, pointsize = 16)
netVisual_circle(cellchat.obj@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")
dev.off()

##############################################################
#### Part C. Visualize of cell-cell communication network ####
##############################################################
###  1. Visualize cell-cell communication mediated by multiple L-R or signaling pathways

# 1.1 Bubble plot
# 1.1.1 Significant interactions of all signaling 
svg(paste0(outdir,"netVisual_bubble_TREM2.svg"), width = 8, height = 16, pointsize = 12)
netVisual_bubble(cellchat.obj, sources.use = source.cluster, targets.use = target.cluster, remove.isolate = FALSE)
dev.off()

