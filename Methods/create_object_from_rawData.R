library(Giotto)
library(TITAN)
library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)
library(lda)
library(pheatmap)
library(RColorBrewer)
library(data.table)

load("complete_giotto_object.RData") # raw data was put into Giotto object by core and then extracted by us

gem@expression[["rna"]][["raw"]] <- list(as.matrix(gem@expression[["rna"]][["raw"]]))
names(gem@expression[["rna"]][["raw"]]) <- "one"
gem@expression[["rna"]][["normalized"]] <- list(as.matrix(gem@expression[["rna"]][["normalized"]]))
names(gem@expression[["rna"]][["normalized"]]) <- "one"

gem@expression[["negprobes"]][["raw"]] <- list(as.matrix(gem@expression[["negprobes"]][["raw"]]))
names(gem@expression[["negprobes"]][["raw"]]) <- "one"

gem@expression[["falsecode"]][["raw"]] <- list(as.matrix(gem@expression[["falsecode"]][["raw"]]))
names(gem@expression[["falsecode"]][["raw"]]) <- "one"

test <- gem@expression$rna$raw

SO <- CreateSeuratObject(counts=test, min.cells=0, min.features=0)

SO@meta.data$nb_clus <- gem@cell_metadata[["rna"]][["nb_clus"]]
SO@meta.data$leiden_clus <- gem@cell_metadata[["rna"]][["leiden_clus"]]
coords <- gem@dimension_reduction[["cells"]][["umap"]][["umap"]][["coordinates"]]
SO[["giotto_UMAP"]] <- CreateDimReducObject(embeddings=coords, key = "gUMAP_", assay = "RNA", global = T)

SO$Patient <- "Unknown"

SO$Patient[which(gem@cell_metadata[["rna"]][["Slide_name"]] %in% c("R1134_S3","R1124_S1"))] <- "Patient_A"
SO$Patient[which(gem@cell_metadata[["rna"]][["Slide_name"]] %in% c("R1124_S3"))] <- "Patient_B"
SO$Patient[which(gem@cell_metadata[["rna"]][["Slide_name"]] %in% c("R1134_S1"))] <- "Patient_C"
SO$Patient[which(gem@cell_metadata[["rna"]][["Slide_name"]] %in%
                   c("R1134_S2","R1124_S2"))] <- "Patient_D"

saveRDS(SO, "cosmx_masterSO.rds")

spatial_locs_raw <- as.data.frame(gem@spatial_locs[["raw"]])
spatial_locs_viz <- as.data.frame(gem@spatial_locs[["viz"]])

PatientA_raw_spatial_locs <- subset(spatial_locs_raw, cell_ID %in% colnames(Patient_A))
PatientA_viz_spatial_locs <- subset(spatial_locs_viz, cell_ID %in% colnames(Patient_A))

PatientB_raw_spatial_locs <- subset(spatial_locs_raw, cell_ID %in% colnames(Patient_B))
PatientB_viz_spatial_locs <- subset(spatial_locs_viz, cell_ID %in% colnames(Patient_B))

PatientC_raw_spatial_locs <- subset(spatial_locs_raw, cell_ID %in% colnames(Patient_C))
PatientC_viz_spatial_locs <- subset(spatial_locs_viz, cell_ID %in% colnames(Patient_C))

PatientD_raw_spatial_locs <- subset(spatial_locs_raw, cell_ID %in% colnames(Patient_D))
PatientD_viz_spatial_locs <- subset(spatial_locs_viz, cell_ID %in% colnames(Patient_D))

##############################################################################
############## Per patient objects ###########################################
##############################################################################

Patient_A <- subset(SO, Patient == "Patient_A")

Patient_A <- NormalizeData(Patient_A)
all.genes <- rownames(Patient_A)
Patient_A <- ScaleData(object = Patient_A, features = all.genes)

Patient_A <- RunPCA(Patient_A, features = all.genes)
ElbowPlot(Patient_A, ndims = 50)

Patient_A <- FindNeighbors(object=Patient_A, dims=1:20)
Patient_A <- FindClusters(object=Patient_A, resolution = 0.3)
Patient_A <- RunUMAP(object=Patient_A, dims=1:20, n.neighbors = 30, min.dist = 0.3)

bx3_cells <- PatientA_raw_spatial_locs$cell_ID[PatientA_raw_spatial_locs$sdimy < -15]
bx2_cells <- PatientA_raw_spatial_locs$cell_ID[PatientA_raw_spatial_locs$sdimy > -15 & PatientA_raw_spatial_locs$sdimx < 10]
bx4_cells <- PatientA_raw_spatial_locs$cell_ID[PatientA_raw_spatial_locs$sdimy > -15 & PatientA_raw_spatial_locs$sdimx > 10]

Patient_A$Timepoint <- "Unassigned"
Patient_A$Timepoint[rownames(Patient_A@meta.data) %in% bx2_cells] <- "Bx2"
Patient_A$Timepoint[rownames(Patient_A@meta.data) %in% bx4_cells] <- "Bx4"
Patient_A$Timepoint[rownames(Patient_A@meta.data) %in% bx3_cells] <- "Bx3"

saveRDS(Patient_A, "Patient_A_SOwTimepoints.rds")

Patient_B <- subset(SO, Patient == "Patient_B")

Patient_B <- NormalizeData(Patient_B)
all.genes <- rownames(Patient_B)
Patient_B <- ScaleData(object = Patient_B, features = all.genes)

Patient_B <- RunPCA(Patient_B, features = all.genes)
ElbowPlot(Patient_B, ndims = 50)

Patient_B <- FindNeighbors(object=Patient_B, dims=1:20)
Patient_B <- FindClusters(object=Patient_B, resolution = 0.4)
Patient_B <- RunUMAP(object=Patient_B, dims=1:20, n.neighbors = 30, min.dist = 0.3)

bx4_cells <- PatientB_raw_spatial_locs$cell_ID[PatientB_raw_spatial_locs$sdimx < 30]
bx1_cells <- PatientB_raw_spatial_locs$cell_ID[PatientB_raw_spatial_locs$sdimx > 30]

Patient_B$Timepoint <- "Unassigned"
Patient_B$Timepoint[rownames(Patient_B@meta.data) %in% bx1_cells] <- "Bx1"
Patient_B$Timepoint[rownames(Patient_B@meta.data) %in% bx4_cells] <- "Bx4"

saveRDS(Patient_B, "Patient_B_SOwTimepoints.rds")

Patient_C <- subset(SO, Patient == "Patient_C")

Patient_C <- NormalizeData(Patient_C)
all.genes <- rownames(Patient_C)
Patient_C <- ScaleData(object = Patient_C, features = all.genes)

Patient_C <- RunPCA(Patient_C, features = all.genes)
ElbowPlot(Patient_C, ndims = 50)

Patient_C <- FindNeighbors(object=Patient_C, dims=1:20)
Patient_C <- FindClusters(object=Patient_C, resolution = 0.5)
Patient_C <- RunUMAP(object=Patient_C, dims=1:20, n.neighbors = 100, min.dist = 0.3)

bx2_cells <- PatientC_raw_spatial_locs$cell_ID[PatientC_raw_spatial_locs$sdimx < 10]
bx1_cells <- PatientC_raw_spatial_locs$cell_ID[PatientC_raw_spatial_locs$sdimx > 10]

Patient_C$Timepoint <- "Unassigned"
Patient_C$Timepoint[rownames(Patient_C@meta.data) %in% bx1_cells] <- "Bx1"
Patient_C$Timepoint[rownames(Patient_C@meta.data) %in% bx2_cells] <- "Bx2"

saveRDS(Patient_C, "Patient_C_SOwTimepoints.rds")

Patient_D <- subset(SO, Patient == "Patient_D")

Patient_D <- NormalizeData(Patient_D)
all.genes <- rownames(Patient_D)
Patient_D <- ScaleData(object = Patient_D, features = all.genes)

Patient_D <- RunPCA(Patient_D, features = all.genes)
ElbowPlot(Patient_D, ndims = 50)

Patient_D <- FindNeighbors(object=Patient_D, dims=1:20)
Patient_D <- FindClusters(object=Patient_D, resolution = 0.3)
Patient_D <- RunUMAP(object=Patient_D, dims=1:20, n.neighbors = 30, min.dist = 0.3)

bx2_cells <- PatientD_raw_spatial_locs$cell_ID[PatientD_raw_spatial_locs$sdimy > -30 & PatientD_raw_spatial_locs$sdimx > 31]
bx1_cells <- PatientD_raw_spatial_locs$cell_ID[PatientD_raw_spatial_locs$sdimy > -30 & PatientD_raw_spatial_locs$sdimx < 31]
bx3_cells <- PatientD_raw_spatial_locs$cell_ID[PatientD_raw_spatial_locs$sdimy < -30]

Patient_D$Timepoint <- "Unassigned"
Patient_D$Timepoint[rownames(Patient_D@meta.data) %in% bx1_cells] <- "Bx1"
Patient_D$Timepoint[rownames(Patient_D@meta.data) %in% bx2_cells] <- "Bx2"
Patient_D$Timepoint[rownames(Patient_D@meta.data) %in% bx3_cells] <- "Bx3"

saveRDS(Patient_D, "Patient_D_SOwTimepoints.rds")


