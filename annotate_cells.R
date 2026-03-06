library(Seurat)
library(dplyr)
library(ggplot2)
library(Signac)

SO <- readRDS("../cosmx_masterSO.rds")

anno_SO <- readRDS(annotation_object) # The annotation object could not be included at this time but is based on the annotation signatures from Gray et al (2025)

DefaultAssay(anno_SO) <- "RNA"

options(future.globals.maxSize = 60*1024^3)

transfer.anchors <- FindTransferAnchors(reference = Ryan_SO, query = SO, dims = 1:20, reference.reduction = "pca")
predictions <- TransferData(anchorset = transfer.anchors, refdata = Ryan_SO$assigned_celltype, dims = 1:20)
SO <- AddMetaData(SO, metadata = predictions)

saveRDS(SO, "cosmx_allCells_newAnno.rds")