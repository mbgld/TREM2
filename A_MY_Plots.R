### 1.Set-up
library(pheatmap); library(Seurat); library(dplyr); library(ggpubr); library(SeuratObject)
outdir = './Results/'
MY.sfj <- readRDS('./Results/MY.sfj')
MY.sfj <- SetIdent(MY.sfj, value = 'subcell.id')

### 2. Feature Plots
# 2.1 Check TREM2 expressing cells
jpeg(filename = paste0(outdir, "TREM2_positive_cells.jpeg"), width = 4000, height =3200, res = 600)
FeaturePlot(MY.sfj, features = 'TREM2', min.cutoff = 1.5, repel = T,
            pt.size = 0.5, label.size = 4, label = T, label.color = "black")
dev.off()

# 2.2 TREM2 expression by tissue
jpeg(filename = paste0(outdir, "TREM2_Positive_by_tissue.jpeg"), width = 4000, height =3200, res = 600)
FeaturePlot(MY.sfj, features = 'TREM2', min.cutoff = 1.5, repel = T,
            pt.size = 0.5, label.size = 4, label = F, label.color = "black", split.by = 'tissue.id')
dev.off()

# 3. Select TREM2 overexpressing cells
TREM2.cell <- WhichCells(MY.sfj, expression = TREM2 > 2.5)
TREM2.mk <- FindMarkers(MY.sfj, ident.1 = TREM2.cell)
TREM2.mk <- TREM2.mk[order(TREM2.mk$p_val),]
write.csv(TREM2.mk, file = paste0(outdir, "/TREM2_vs_other_Mc_mk.csv"))
TREM2.mk <- read.csv(file = paste0(outdir, "/TREM2_vs_other_Mc_mk.csv"), row.names = 1)

# 4 TREM2 cells Volcano plot (1) with annotation
detach("package:Seurat", unload=TRUE); library(EnhancedVolcano)

tiff(paste0(outdir, "Enhanced_Volcano_TREM2_cells_Selected_gene.tiff"), width = 8, height = 7, units = 'in', res = 300)
EnhancedVolcano(TREM2.mk, 
                lab = rownames(TREM2.mk), 
                x='avg_log2FC', y='p_val_adj',
                selectLab = c('A2M', 'LIPA', 'SPP1', 'MMP9',
                            'TREM2','PIGR', 'PLA2G7',
                            'PPARG', 'CEBPB', 'FABP4'),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'Padj'),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 1.0,
                labSize = 4.0,
                labCol = 'black',
                labFace = 'italic',
                boxedLabels = FALSE,
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')
dev.off()

