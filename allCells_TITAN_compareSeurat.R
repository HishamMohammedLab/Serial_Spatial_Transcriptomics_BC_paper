library(Seurat)
library(TITAN)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(bluster)
library(rstatix)
library(effectsize)

SO <- readRDS("data/cosmx_masterSO.rds")

SO <- subset(SO, nCount_RNA > 200) # filtering out low count cells

## All instances of TITAN are commented out because they take a long time to run, but just uncomment if you would like to run it

#runLDA(SO, ntopics = seq(10,100, by=10), parallel = T, outDir = "data/allCosmx_Models", cores = 10)

#pdf("allCosmx_Models/TITAN_elbowPlot.pdf")
#LDAelbowPlot("allCosmx_Models/", SO) + theme_bw() + theme(panel.border = element_blank())
#dev.off()

model <- readRDS("data/allCosmx_Models/Model_50topics.rds")

top_genes <- TopTopicGenes(model, ngenes = 50)

SO <- addTopicsToSeuratObject(model, SO)

anno_SO <- readRDS("data/cosmx_allCells_newAnno.rds")

SO$new_celltype <- anno_SO@meta.data[rownames(SO@meta.data), "predicted.id"]

SO <- RunUMAP(SO, reduction = "lda", dims = 1:50, reduction.name = "TITANumap")

patient_colors <- list(Patient_A = "#C1121F", Patient_B = "#2D6A4F", Patient_C = "#0077B6", Patient_D = "#E6AC00")

#Figure 1B All cells by Patient

png("plots/allCosmx_50topics_UMAP_byPatient.png")
DimPlot(SO, reduction = "TITANumap", group.by = "Patient", cols = patient_colors)
dev.off()

png("plots/allCosmx_50topics_UMAP_byCelltypes.png")
DimPlot(SO, reduction = "TITANumap", group.by = "new_celltype")
dev.off()

SO$grouped_epi_celltype <- SO$new_celltype
SO$grouped_epi_celltype[SO$grouped_epi_celltype %in% c("basal_myoepithelial", "cancer", "luminal_asp", "luminal_hs")] <- "epithelial"

retro_muted_colors <- c(
  "Cancer"                 = "#82B0AF",
  "epithelial"             = "#E69138",
  "tcell"           = "#E9C44D",
  "bcell"           = "#D95F5F",
  "myeloid"                = "#B7B0A8",
  "plasma"                 = "#8B426F",
  "fibroblast"             = "#A47C66",
  "endothelial_vascular" = "#5D85AD",
  "endothelial_lymphatic"= "#F198A4",
  "adipocyte"              = "#5D8E4D",
  "pericyte"               = "#A5799D"
)

# Figure 1B All Cells by Cell types

png("allCosmx_50topics_UMAP_grouped_Epi_celltypes.png")
DimPlot(SO, reduction = "TITANumap", group.by = "grouped_epi_celltype", cols = retro_muted_colors)
dev.off()

################################################################################
########## Compare with Seurat #################################################
################################################################################

SO <- NormalizeData(SO)
all.genes <- rownames(SO)
SO <- ScaleData(object = SO, features = all.genes)

SO <- RunPCA(SO, features = all.genes)
ElbowPlot(SO, ndims = 50)

SO <- FindNeighbors(object=SO, dims=1:20)
SO <- FindClusters(object=SO, resolution = 0.8)
SO <- RunUMAP(object=SO, dims=1:20, n.neighbors = 30, min.dist = 0.3)

#Extended Figure 1D

pdf("plots/allCosmx_TITANumap_byCoverage.pdf")
FeaturePlot(SO, features = "nCount_RNA", reduction = "TITANumap", max.cutoff = 1000, min.cutoff = 200, cols = c("grey90", "red"))
dev.off()

pdf("plots/allCosmx_SeuratUMAP_byCoverage.pdf")
FeaturePlot(SO, features = "nCount_RNA", reduction = "umap", max.cutoff = 1000, min.cutoff = 200, cols = c("grey90", "red"))
dev.off()

VlnPlot(SO, features = "nCount_RNA", pt.size = 0)

SO$covBin <- "unassigned"

SO$covBin[SO$nCount_RNA <= quantile(SO$nCount_RNA, 0.2)] <- "1"
SO$covBin[SO$nCount_RNA > quantile(SO$nCount_RNA, 0.2) & SO$nCount_RNA <= quantile(SO$nCount_RNA, 0.4)] <- "2"
SO$covBin[SO$nCount_RNA > quantile(SO$nCount_RNA, 0.4) & SO$nCount_RNA <= quantile(SO$nCount_RNA, 0.6)] <- "3"
SO$covBin[SO$nCount_RNA > quantile(SO$nCount_RNA, 0.6) & SO$nCount_RNA <= quantile(SO$nCount_RNA, 0.8)] <- "4"
SO$covBin[SO$nCount_RNA > quantile(SO$nCount_RNA, 0.8) & SO$nCount_RNA <= quantile(SO$nCount_RNA, 1.0)] <- "5"

# Calculating Coverage Purity

Seurat_purity <- neighborPurity(SO@reductions$pca@cell.embeddings, clusters = SO$covBin)
saveRDS(Seurat_purity, "data/allCosmx_Seurat_purity.rds")

TITAN_purity <- neighborPurity(SO@reductions$lda@cell.embeddings, clusters = SO$covBin)
saveRDS(TITAN_purity, "data/allCosmx_TITAN_purity.rds")

Seurat_purity <- readRDS("allCosmx_Seurat_purity.rds")
TITAN_purity  <- readRDS("allCosmx_TITAN_purity.rds")

purity_plot_df <- data.frame(Purity = c(TITAN_purity$purity, Seurat_purity$purity), 
                             Method = c(rep("TITAN", nrow(TITAN_purity)), rep("Seurat", nrow(Seurat_purity))),
                             celltype = c(TITAN_purity$maximum, Seurat_purity$maximum))

# Extended Figure 1G                             
                                                                                                                                                                                   
pdf("plots/allCosmx_coverage_purity_clean_noFill.pdf")
ggplot(data=purity_plot_df, mapping=aes(x=Method,y=Purity,color=Method)) + geom_boxplot(linewidth=1.5) + theme_classic() + scale_fill_manual(values=c("#2C7BB6", "#F26B22")) + NoLegend()
dev.off()
                                                                                                                                                                                   
SO <- FindNeighbors(object=SO, dims=1:20, reduction = "pca", graph.name = "pca.graph")
SO <- FindClusters(object=SO, resolution = 1.0, graph.name = "pca.graph", cluster.name = "pca.clusters")
                                                                                                                                                                                   
SO <- FindNeighbors(object=SO, dims=1:50, reduction = "lda", graph.name = "TITAN.graph")
SO <- FindClusters(object=SO, resolution = 1.0, graph.name = "TITAN.graph", cluster.name = "TITAN.clusters")
                                                                                                                                                                                   
SO <- subset(SO, TITAN.clusters %in% c("0","1","2","3","4","5","6","7","8","9","10","11","12","13", "14"))
                                                                                                                                                                                  
#SO <- subset(SO, pca.clusters != "14")
                                                                                                                                                                                   
Idents(SO) <- "pca.clusters"

# Extended Figure 1E and 1F
                                                                                                                                                                                   
pdf("allCosmx_seuratPCAclusters_byCoverage_VlnPlot.pdf")
VlnPlot(SO, features = "nCount_RNA", pt.size = 0) + geom_boxplot()
dev.off()
                                                                                                                                                                                   
Idents(SO) <- "TITAN.clusters"
                                                                                                                                                                                   
pdf("allCosmx_TITANclusters_byCoverage_VlnPlot.pdf")
VlnPlot(SO, features = "nCount_RNA", pt.size = 0) + geom_boxplot()
dev.off()
                                                                                                                                                                                   
SO$pca.clusters <- factor(SO$pca.clusters, levels = c("0","1","2","3","4","5","6","7","8","9","10","11","12"))
                                                                                                                                                                                   
Idents(SO) <- "pca.clusters"
                                                                                                                                                                                   
blues <- colorRampPalette(brewer.pal(9, "Blues"))
                                                                                                                                                                                   
SO$logCount_RNA <- log(SO$nCount_RNA)

# Extended Figure 1H
                                                                                                                                                                                   
pdf("allCosmx_seuratPCAclusters_byLogCoverage_BoxPlot_clean.pdf")
ggplot(SO@meta.data, aes(x=pca.clusters, y=logCount_RNA, fill=pca.clusters)) + geom_boxplot() + scale_fill_manual(values=blues(14)) + theme_classic() + NoLegend()
dev.off()
                                                                                                                                                                                   
SO$TITAN.clusters <- factor(SO$TITAN.clusters, levels = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14"))
                                                                                                                                                                                   
Idents(SO) <- "TITAN.clusters"
                                                                                                                                                                                   
oranges <- colorRampPalette(brewer.pal(9, "Oranges"))
                                                                                                                                                                                   
pdf("allCosmx_TITANclusters_byLogCoverage_BoxPlot_clean.pdf")
ggplot(SO@meta.data, aes(x=TITAN.clusters, y=logCount_RNA, fill = TITAN.clusters)) + geom_boxplot() + scale_fill_manual(values=oranges(15)) + theme_classic() + NoLegend()
dev.off()
   
# Extended Figure 1C
                                                                                                                                                                                
pdf("allCosmx_PCAclusters_UMAP.pdf")
DimPlot(SO, reduction = "umap", group.by = "pca.clusters")
dev.off()
                                                                                                                                                                                   
pdf("allCosmx_TITANclusters_TITANumap.pdf")
DimPlot(SO, reduction = "TITANumap", group.by = "TITAN.clusters")
dev.off()
                                                                                                                                                                                   
Idents(SO) <- "Patient"
                                                                                                                                                                                   
pdf("Coverage_boxPlot_byPatient.pdf")
VlnPlot(SO, features = "nCount_RNA", pt.size = 0) + geom_boxplot()
dev.off()
                                                                                                                                                                                   
pdf("UniqueFeature_boxPlot_byPatient.pdf")
VlnPlot(SO, features = "nFeature_RNA", pt.size = 0) + geom_boxplot()
dev.off()
                                                                                                                                                                                   