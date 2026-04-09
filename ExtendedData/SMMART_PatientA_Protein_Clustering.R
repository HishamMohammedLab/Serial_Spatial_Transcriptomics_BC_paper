# ==============================================================================
# SMMART PROTEIN - Patient A Clustering (Corrected Mapping)
# Cluster cells based on protein markers
# 
# CONFIRMED SLIDE MAPPING:
#   Slide 1 = Patient A Bx2
#   Slide 2 = Patient A Bx4
# ==============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(viridis)

# ------------------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("PATIENT A PROTEIN CLUSTERING\n")
cat("============================================================\n\n")

output_dir <- "PatientA_Figures/Protein_Analysis/Clustering/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load protein data
cat("Loading protein Seurat object...\n")
obj_protein <- readRDS("data/SMMART_protein_cosmx_SO.RDS")
obj_protein <- UpdateSeuratObject(obj_protein)

cat("Total cells in object:", ncol(obj_protein), "\n\n")

# ------------------------------------------------------------------------------
# CONFIRMED SAMPLE MAPPING
# ------------------------------------------------------------------------------

cat("CONFIRMED SLIDE MAPPING:\n")
cat("  Slide 1 = Patient A Bx2\n")
cat("  Slide 2 = Patient A Bx4\n\n")

# Subset to Patient A slides 1 and 2 only
obj <- subset(obj_protein, subset = slide_ID_numeric %in% c(1, 2))

# Add corrected timepoint labels
obj$Timepoint <- ifelse(obj$slide_ID_numeric == 1, "Bx2", "Bx4")
obj$Patient <- "A"
obj$Sample <- paste0("A_", obj$Timepoint)

cat("Patient A cells (Bx2 + Bx4):", ncol(obj), "\n")
cat("  Bx2 (slide 1):", sum(obj$Timepoint == "Bx2"), "cells\n")
cat("  Bx4 (slide 2):", sum(obj$Timepoint == "Bx4"), "cells\n\n")

# ------------------------------------------------------------------------------
# CANCER FOV ANNOTATION
# ------------------------------------------------------------------------------

cat("Adding Cancer FOV annotations...\n")

# User-identified cancer FOVs (determined via CNV analysis of RNA)
cancer_fovs_Bx2 <- c(4, 5, 6, 13, 14, 15, 37, 38, 39, 40, 41, 42, 43)
cancer_fovs_Bx4 <- c(15, 16, 18, 19, 20, 23, 24, 26, 27)

# Create Cancer_FOV annotation
obj$Cancer_FOV <- ifelse(
  (obj$Timepoint == "Bx2" & obj$fov %in% cancer_fovs_Bx2) |
  (obj$Timepoint == "Bx4" & obj$fov %in% cancer_fovs_Bx4),
  "Cancer_FOV",
  "Other_FOV"
)

# Create combined FOV label for more detail
obj$FOV_Type <- paste0(obj$Timepoint, "_", obj$Cancer_FOV)

cat("  Cancer FOV cells:", sum(obj$Cancer_FOV == "Cancer_FOV"), "\n")
cat("  Other FOV cells:", sum(obj$Cancer_FOV == "Other_FOV"), "\n\n")

# Breakdown by timepoint
cat("  Bx2 Cancer FOVs:", sum(obj$Timepoint == "Bx2" & obj$Cancer_FOV == "Cancer_FOV"), "\n")
cat("  Bx2 Other FOVs:", sum(obj$Timepoint == "Bx2" & obj$Cancer_FOV == "Other_FOV"), "\n")
cat("  Bx4 Cancer FOVs:", sum(obj$Timepoint == "Bx4" & obj$Cancer_FOV == "Cancer_FOV"), "\n")
cat("  Bx4 Other FOVs:", sum(obj$Timepoint == "Bx4" & obj$Cancer_FOV == "Other_FOV"), "\n\n")

# ------------------------------------------------------------------------------
# DEFINE MARKERS FOR CLUSTERING
# ------------------------------------------------------------------------------

# Get all protein features (exclude channels and negprobes)
all_features <- rownames(obj)
proteins <- all_features[!grepl("^Channel-|^Ms IgG|^Rb IgG", all_features)]

cat("Proteins for clustering:", length(proteins), "\n\n")

# Key lineage markers for annotation
lineage_markers <- list(
  T_cell = c("CD3", "CD45"),
  CD4_T = c("CD4"),
  CD8_T = c("CD8"),
  Treg = c("FOXP3"),
  T_exhaustion = c("PD-1", "LAG3", "Tim-3", "CD39"),
  T_activation = c("4-1BB", "ICOS", "GITR"),
  T_cytotoxic = c("GZMA", "GZMB"),
  T_memory = c("CD45RA", "CCR7", "TCF7"),
  
  Myeloid = c("CD68", "CD14", "CD11b"),
  DC = c("CD11c", "HLA-DR"),
  M2_mac = c("CD163"),
  M1_mac = c("iNOS"),
  
  B_cell = c("CD19", "CD20"),
  Plasma = c("CD138", "CD38"),
  
  NK = c("CD56", "CD16"),
  
  Epithelial = c("EpCAM"),
  Stromal = c("Vimentin", "SMA", "Fibronectin"),
  Endothelial = c("CD31", "CD34"),
  
  Checkpoint_ligands = c("PD-L1", "PD-L2", "B7-H3", "VISTA"),
  Proliferation = c("Ki-67"),
  Signaling = c("STING", "NF-kB p65", "Beta-catenin")
)

# Flatten for heatmap
all_lineage_markers <- unique(unlist(lineage_markers))
all_lineage_markers <- intersect(all_lineage_markers, proteins)

# ------------------------------------------------------------------------------
# PREPROCESSING & CLUSTERING
# ------------------------------------------------------------------------------

cat("Running clustering pipeline...\n\n")

# Normalize data
cat("1. Normalizing (CLR)...\n")
obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2, verbose = FALSE)

# Scale data
cat("2. Scaling...\n")
obj <- ScaleData(obj, features = proteins, verbose = FALSE)

# PCA
cat("3. Running PCA...\n")
obj <- RunPCA(obj, features = proteins, npcs = 30, verbose = FALSE)

# Check variance explained
pct_var <- obj@reductions$pca@stdev^2 / sum(obj@reductions$pca@stdev^2) * 100
cum_var <- cumsum(pct_var)
n_pcs <- min(which(cum_var > 80))
cat("   PCs for 80% variance:", n_pcs, "\n")

n_pcs_use <- min(20, n_pcs + 5)
cat("   Using", n_pcs_use, "PCs\n")

# Find neighbors
cat("4. Finding neighbors...\n")
obj <- FindNeighbors(obj, dims = 1:n_pcs_use, verbose = FALSE)

# Cluster at multiple resolutions
cat("5. Clustering at multiple resolutions...\n")
resolutions <- c(0.3, 0.5, 0.8, 1.0, 1.5)

for (res in resolutions) {
  obj <- FindClusters(obj, resolution = res, verbose = FALSE)
  cat("   Resolution", res, ":", length(unique(obj@meta.data[[paste0("RNA_snn_res.", res)]])), "clusters\n")
}

# Set default to resolution 0.5
Idents(obj) <- "RNA_snn_res.0.5"
obj$seurat_clusters <- obj$RNA_snn_res.0.5

n_clusters <- length(unique(obj$seurat_clusters))
cat("\nUsing resolution 0.5 with", n_clusters, "clusters\n")

# UMAP
cat("6. Running UMAP...\n")
obj <- RunUMAP(obj, dims = 1:n_pcs_use, verbose = FALSE)

# ------------------------------------------------------------------------------
# COLOR PALETTES (DARKER)
# ------------------------------------------------------------------------------

# Darker cluster colors
cluster_colors <- c(
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
  "#A6761D", "#666666", "#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3",
  "#FDB462", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#FFED6F",
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628"
)[1:n_clusters]
names(cluster_colors) <- levels(obj$seurat_clusters)

# Timepoint colors (darker)
timepoint_colors <- c("Bx2" = "#2166AC", "Bx4" = "#B2182B")

# Cancer FOV colors (distinct)
cancer_fov_colors <- c("Cancer_FOV" = "#E41A1C", "Other_FOV" = "#4DAF4A")

# FOV Type colors (4 categories)
fov_type_colors <- c(
  "Bx2_Cancer_FOV" = "#B2182B",   # Dark red
  "Bx2_Other_FOV" = "#2166AC",    # Dark blue
  "Bx4_Cancer_FOV" = "#D6604D",   # Light red
  "Bx4_Other_FOV" = "#4393C3"     # Light blue
)

# ------------------------------------------------------------------------------
# VISUALIZATIONS
# ------------------------------------------------------------------------------

cat("\n7. Creating visualizations...\n\n")

# 7.1 UMAP by cluster
p_umap_cluster <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters",
                           cols = cluster_colors, pt.size = 0.3) +
  ggtitle("Clusters (res=0.5)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 7.2 UMAP by timepoint
p_umap_timepoint <- DimPlot(obj, reduction = "umap", group.by = "Timepoint",
                             cols = timepoint_colors, pt.size = 0.3) +
  ggtitle("Timepoint") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 7.3 UMAP by Cancer FOV
p_umap_cancer <- DimPlot(obj, reduction = "umap", group.by = "Cancer_FOV",
                          cols = cancer_fov_colors, pt.size = 0.3) +
  ggtitle("Cancer vs Other FOVs") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 7.4 UMAP by FOV Type (4 categories)
p_umap_fovtype <- DimPlot(obj, reduction = "umap", group.by = "FOV_Type",
                           cols = fov_type_colors, pt.size = 0.3) +
  ggtitle("FOV Type by Timepoint") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Combined UMAP plot - 2x2 grid
p_umap_combined <- (p_umap_cluster | p_umap_timepoint) / (p_umap_cancer | p_umap_fovtype) +
  plot_annotation(
    title = "Patient A Protein Clustering",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

ggsave(file.path(output_dir, "PatientA_UMAP_Overview.pdf"),
       p_umap_combined, width = 14, height = 12)
cat("Saved: PatientA_UMAP_Overview.pdf\n")

# 7.5 Larger Cancer FOV UMAP
p_cancer_large <- DimPlot(obj, reduction = "umap", group.by = "Cancer_FOV",
                           cols = cancer_fov_colors, pt.size = 0.5) +
  ggtitle("Cancer FOVs vs Other FOVs") +
  labs(subtitle = paste0("Cancer FOV: ", sum(obj$Cancer_FOV == "Cancer_FOV"), 
                         " cells | Other: ", sum(obj$Cancer_FOV == "Other_FOV"), " cells")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

ggsave(file.path(output_dir, "PatientA_UMAP_CancerFOVs.pdf"),
       p_cancer_large, width = 10, height = 8)
cat("Saved: PatientA_UMAP_CancerFOVs.pdf\n")

# 7.6 Split by timepoint
p_split_timepoint <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters",
                              split.by = "Timepoint", cols = cluster_colors, pt.size = 0.3) +
  ggtitle("Clusters Split by Timepoint") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "PatientA_UMAP_SplitTimepoint.pdf"),
       p_split_timepoint, width = 14, height = 6)
cat("Saved: PatientA_UMAP_SplitTimepoint.pdf\n")

# 7.7 Split by Cancer FOV
p_split_cancer <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters",
                           split.by = "Cancer_FOV", cols = cluster_colors, pt.size = 0.3) +
  ggtitle("Clusters Split by Cancer vs Other FOVs") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "PatientA_UMAP_SplitCancerFOV.pdf"),
       p_split_cancer, width = 14, height = 6)
cat("Saved: PatientA_UMAP_SplitCancerFOV.pdf\n")

# 7.8 Feature plots for key markers (darker gradient)
key_markers <- c("CD3", "CD4", "CD8", "CD68", "CD163", "EpCAM", 
                 "CD20", "CD56", "Vimentin", "PD-L1", "PD-1", "Ki-67")
key_markers <- intersect(key_markers, proteins)

p_features <- FeaturePlot(obj, features = key_markers, 
                           ncol = 4, pt.size = 0.2,
                           cols = c("grey85", "darkred")) &
  theme(plot.title = element_text(size = 10))

ggsave(file.path(output_dir, "PatientA_UMAP_KeyMarkers.pdf"),
       p_features, width = 16, height = 12)
cat("Saved: PatientA_UMAP_KeyMarkers.pdf\n")

# ------------------------------------------------------------------------------
# CLUSTER MARKER EXPRESSION HEATMAP
# ------------------------------------------------------------------------------

cat("\n8. Creating cluster heatmap...\n")

# Calculate mean expression per cluster
expr <- GetAssayData(obj, layer = "data")
cluster_means <- sapply(levels(obj$seurat_clusters), function(cl) {
  cells <- WhichCells(obj, idents = cl)
  if (length(cells) > 1) {
    rowMeans(as.matrix(expr[all_lineage_markers, cells]))
  } else {
    expr[all_lineage_markers, cells]
  }
})

cluster_means_scaled <- t(scale(t(cluster_means)))

# Cluster annotation bar
cluster_sizes <- table(obj$seurat_clusters)
col_labels <- paste0("C", names(cluster_sizes), "\n(n=", cluster_sizes, ")")

# Create heatmap
pdf(file.path(output_dir, "PatientA_Cluster_Heatmap.pdf"), width = 12, height = 14)

pheatmap(
  cluster_means_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  labels_col = col_labels,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  breaks = seq(-2, 2, length.out = 101),
  main = "Patient A: Cluster Marker Expression",
  fontsize_row = 9,
  fontsize_col = 10
)

dev.off()
cat("Saved: PatientA_Cluster_Heatmap.pdf\n")

# ------------------------------------------------------------------------------
# DOT PLOT
# ------------------------------------------------------------------------------

cat("\n9. Creating dot plot...\n")

p_dot <- DotPlot(obj, features = all_lineage_markers, group.by = "seurat_clusters") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Marker Expression by Cluster") +
  scale_color_gradient2(low = "navy", mid = "white", high = "firebrick3", midpoint = 0)

ggsave(file.path(output_dir, "PatientA_Cluster_DotPlot.pdf"),
       p_dot, width = 14, height = 16)
cat("Saved: PatientA_Cluster_DotPlot.pdf\n")

# ------------------------------------------------------------------------------
# CLUSTER COMPOSITION ANALYSIS
# ------------------------------------------------------------------------------

cat("\n10. Analyzing cluster composition...\n")

# Cluster proportions by timepoint
comp_timepoint <- obj@meta.data %>%
  group_by(Timepoint, seurat_clusters) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(Timepoint) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p_comp_tp <- ggplot(comp_timepoint, aes(x = Timepoint, y = pct, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "Cluster Composition: Bx2 vs Bx4", x = NULL, y = "% of Cells", fill = "Cluster") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Cluster proportions by Cancer FOV
comp_cancer <- obj@meta.data %>%
  group_by(Cancer_FOV, seurat_clusters) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(Cancer_FOV) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p_comp_cancer <- ggplot(comp_cancer, aes(x = Cancer_FOV, y = pct, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "Cluster Composition: Cancer vs Other FOVs", x = NULL, y = "% of Cells", fill = "Cluster") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Combined composition
p_comp_combined <- p_comp_tp + p_comp_cancer
ggsave(file.path(output_dir, "PatientA_Cluster_Composition.pdf"),
       p_comp_combined, width = 14, height = 6)
cat("Saved: PatientA_Cluster_Composition.pdf\n")

# ------------------------------------------------------------------------------
# CLUSTER SUMMARY TABLES
# ------------------------------------------------------------------------------

cat("\n11. Creating cluster summary tables...\n")

# Summary markers
summary_markers <- c("CD3", "CD4", "CD8", "FOXP3", "CD68", "CD163", "iNOS", "CD11c",
                     "EpCAM", "CD20", "CD56", "Vimentin", "CD31",
                     "PD-1", "PD-L1", "LAG3", "Tim-3", "Ki-67", "HLA-DR", "B7-H3")
summary_markers <- intersect(summary_markers, proteins)

# Mean expression per cluster
cluster_summary <- data.frame(
  Cluster = levels(obj$seurat_clusters),
  N_Cells = as.numeric(table(obj$seurat_clusters))
)

for (marker in summary_markers) {
  means <- sapply(levels(obj$seurat_clusters), function(cl) {
    cells <- WhichCells(obj, idents = cl)
    mean(expr[marker, cells])
  })
  cluster_summary[[marker]] <- round(means, 2)
}

# Add percentage in Cancer FOV vs Other
cluster_cancer_pct <- obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarize(Pct_CancerFOV = round(mean(Cancer_FOV == "Cancer_FOV") * 100, 1)) %>%
  pull(Pct_CancerFOV)

cluster_summary$Pct_in_CancerFOV <- cluster_cancer_pct

# Add percentage by timepoint
cluster_bx2_pct <- obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarize(Pct_Bx2 = round(mean(Timepoint == "Bx2") * 100, 1)) %>%
  pull(Pct_Bx2)

cluster_summary$Pct_Bx2 <- cluster_bx2_pct
cluster_summary$Pct_Bx4 <- 100 - cluster_bx2_pct

# Add empty annotation column
cluster_summary$Annotation <- ""

write.csv(cluster_summary, 
          file.path(output_dir, "PatientA_Cluster_Summary.csv"),
          row.names = FALSE)
cat("Saved: PatientA_Cluster_Summary.csv\n")

# Bx2 vs Bx4 comparison
bx2_bx4_comp <- comp_timepoint %>%
  select(Timepoint, seurat_clusters, pct) %>%
  pivot_wider(names_from = Timepoint, values_from = pct, values_fill = 0) %>%
  mutate(
    Change = Bx4 - Bx2,
    Pct_Change = round((Bx4 - Bx2) / pmax(Bx2, 0.1) * 100, 1),
    Direction = case_when(
      Pct_Change > 20 ~ "UP",
      Pct_Change < -20 ~ "DOWN",
      TRUE ~ "STABLE"
    )
  )

write.csv(bx2_bx4_comp, 
          file.path(output_dir, "PatientA_Bx2vsBx4_Clusters.csv"),
          row.names = FALSE)
cat("Saved: PatientA_Bx2vsBx4_Clusters.csv\n")

# Cancer FOV enrichment
cancer_enrichment <- comp_cancer %>%
  select(Cancer_FOV, seurat_clusters, pct) %>%
  pivot_wider(names_from = Cancer_FOV, values_from = pct, values_fill = 0) %>%
  mutate(
    Enrichment_CancerFOV = round(Cancer_FOV / pmax(Other_FOV, 0.1), 2),
    Direction = case_when(
      Enrichment_CancerFOV > 1.5 ~ "ENRICHED_in_Cancer",
      Enrichment_CancerFOV < 0.67 ~ "ENRICHED_in_Other",
      TRUE ~ "SIMILAR"
    )
  )

write.csv(cancer_enrichment,
          file.path(output_dir, "PatientA_CancerFOV_Enrichment.csv"),
          row.names = FALSE)
cat("Saved: PatientA_CancerFOV_Enrichment.csv\n")

# ------------------------------------------------------------------------------
# PRINT SUMMARIES
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("CLUSTER SUMMARY\n")
cat("============================================================\n\n")

print(cluster_summary[, c("Cluster", "N_Cells", "Pct_in_CancerFOV", "Pct_Bx2", "Pct_Bx4", 
                          summary_markers[1:min(6, length(summary_markers))])])

cat("\n\nBx2 vs Bx4 Cluster Changes:\n")
print(bx2_bx4_comp)

cat("\n\nCancer FOV Enrichment:\n")
print(cancer_enrichment)

# ------------------------------------------------------------------------------
# SAVE CLUSTERED OBJECT
# ------------------------------------------------------------------------------

cat("\n12. Saving clustered object...\n")
saveRDS(obj, file.path(output_dir, "PatientA_Clustered.rds"))
cat("Saved: PatientA_Clustered.rds\n")

# ------------------------------------------------------------------------------
# INSTRUCTIONS
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("NEXT STEPS: CLUSTER ANNOTATION\n")
cat("============================================================\n\n")

cat("1. Review the output files:\n")
cat("   - PatientA_UMAP_Overview.pdf: Cluster, timepoint, cancer FOV views\n")
cat("   - PatientA_UMAP_CancerFOVs.pdf: Large cancer vs other FOV plot\n")
cat("   - PatientA_UMAP_KeyMarkers.pdf: Key marker expression\n")
cat("   - PatientA_Cluster_Heatmap.pdf: Marker heatmap\n")
cat("   - PatientA_Cluster_DotPlot.pdf: Dot plot of all markers\n")
cat("   - PatientA_Cluster_Summary.csv: Mean expression + Cancer FOV %\n\n")

cat("2. Key annotation hints:\n")
cat("   - Clusters enriched in Cancer FOVs (high Pct_in_CancerFOV) = likely tumor/TME\n")
cat("   - High EpCAM + low CD45 = Epithelial/Cancer\n")
cat("   - High CD3 = T cells; check CD4/CD8 for subset\n")
cat("   - High CD68 = Macrophages; check CD163 (M2) vs iNOS (M1)\n")
cat("   - High Vimentin = Stromal\n\n")

cat("3. Fill in 'Annotation' column in PatientA_Cluster_Summary.csv\n")
cat("4. Upload annotated CSV for refined analysis!\n\n")

cat("============================================================\n")
cat("Done! Output directory:", output_dir, "\n")
cat("============================================================\n")
