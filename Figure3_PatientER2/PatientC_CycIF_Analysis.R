# ============================================================================
# Patient C CycIF Validation Analysis
# For refining CAF zones and validating CosMx findings
# ============================================================================

library(tidyverse)
library(data.table)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(patchwork)
library(scales)

set.seed(42)

# ============================================================================
# PATHS
# ============================================================================

base_path <- "supplementary_input_data"
data_path <- file.path(base_path, "CycIF_Data")
output_path <- file.path(base_path, "New Figures for paper/CycIF_Validation/PatientC")

if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

cat("Output directory:", output_path, "\n")

# ============================================================================
# SECTION 1: LOAD PATIENT C DATA
# ============================================================================

cat("\n========== LOADING PATIENT C CycIF DATA ==========\n\n")

# Load coordinates
coord_df <- fread(file.path(data_path, "features_SMT130-4_CentroidXY.csv"))
setnames(coord_df, "V1", "cell_id")
# Convert to microns (assuming 0.35 um/pixel)
coord_df[, X_um := DAPI_X * 0.35]
coord_df[, Y_um := DAPI_Y * 0.35]

cat("Coordinates loaded:", nrow(coord_df), "cells\n")
cat("Scenes:", paste(unique(coord_df$scene), collapse = ", "), "\n")

# Load protein intensities
protein_df <- fread(file.path(data_path, "features_SMT130-4_FilteredMeanIntensity_DAPI11_Nuclei1000.csv"))
setnames(protein_df, "UNIQID", "cell_id")

cat("Protein data loaded:", nrow(protein_df), "cells\n")

# Get marker names (exclude cell_id and QC columns)
protein_markers <- colnames(protein_df)[!colnames(protein_df) %in% c("cell_id")]
# Remove QC markers (R5Qc*)
protein_markers <- protein_markers[!grepl("^R5Qc", protein_markers)]

cat("Protein markers:", length(protein_markers), "\n")

# Merge coordinates and protein data
merged_df <- merge(protein_df, coord_df[, .(cell_id, DAPI_X, DAPI_Y, X_um, Y_um, scene)],
                   by = "cell_id", all.x = TRUE)

cat("Merged data:", nrow(merged_df), "cells with both coordinates and protein\n")

# ============================================================================
# SECTION 2: CELL TYPE CLASSIFICATION
# ============================================================================

cat("\n========== CELL TYPE CLASSIFICATION ==========\n\n")

# Define marker thresholds using quantiles
# Log transform for better distribution
for (marker in protein_markers) {
  merged_df[[paste0(marker, "_log")]] <- log1p(merged_df[[marker]])
}

# Calculate z-scores for key markers
key_markers <- c("CK8_Ring", "CK19_Ring", "CK7_Ring", "ER_Nuclei",  # Epithelial
                 "CD45_Ring", "CD3_Ring", "CD68_Ring", "CD20_Ring",  # Immune
                 "ColI_Ring", "ColIV_Ring", "Vim_Ring", "aSMA_Ring", "PDPN_Ring")  # Stroma

for (marker in key_markers) {
  if (marker %in% colnames(merged_df)) {
    merged_df[[paste0(marker, "_zscore")]] <- scale(log1p(merged_df[[marker]]))[,1]
  }
}

# Simple cell type classification
merged_df[, CellType := "Other"]

# Epithelial/Cancer: high CK8 or CK19
epi_threshold <- quantile(merged_df$CK8_Ring, 0.7, na.rm = TRUE)
merged_df[CK8_Ring > epi_threshold | CK19_Ring > epi_threshold, CellType := "Epithelial"]

# Immune: high CD45
cd45_threshold <- quantile(merged_df$CD45_Ring, 0.8, na.rm = TRUE)
merged_df[CD45_Ring > cd45_threshold, CellType := "Immune"]

# Fibroblast/Stroma: high Vim or ColI, low CK and CD45
vim_threshold <- quantile(merged_df$Vim_Ring, 0.7, na.rm = TRUE)
col_threshold <- quantile(merged_df$ColI_Ring, 0.7, na.rm = TRUE)
merged_df[CellType == "Other" & (Vim_Ring > vim_threshold | ColI_Ring > col_threshold),
          CellType := "Fibroblast"]

cat("\nCell type distribution:\n")
print(table(merged_df$CellType))
print(round(100 * prop.table(table(merged_df$CellType)), 1))

# ============================================================================
# SECTION 3: CAF SUBTYPE CLASSIFICATION
# ============================================================================

cat("\n========== CAF SUBTYPE CLASSIFICATION ==========\n\n")

# Focus on fibroblasts - make a copy to avoid reference issues
fib_df <- copy(merged_df[CellType == "Fibroblast"])
cat("Fibroblasts for CAF analysis:", nrow(fib_df), "\n")

# CAF subtype markers:
# - mCAF (myofibroblastic): aSMA+, ColI high
# - iCAF (inflammatory): lower aSMA, may have different collagen pattern
# - vCAF (vascular): CD31 associated
# - apCAF (antigen presenting): potentially MHC markers

# Calculate CAF scores
fib_df[, mCAF_score := scale(log1p(aSMA_Ring))[,1] + scale(log1p(ColI_Ring))[,1]]
fib_df[, ColI_score := scale(log1p(ColI_Ring))[,1]]
fib_df[, aSMA_score := scale(log1p(aSMA_Ring))[,1]]

# Simple mCAF vs non-mCAF classification
asma_med <- median(fib_df$aSMA_Ring, na.rm = TRUE)
fib_df[, CAF_subtype := ifelse(aSMA_Ring > asma_med, "mCAF-like", "iCAF-like")]

cat("\nCAF subtype distribution:\n")
print(table(fib_df$CAF_subtype))

# ============================================================================
# SECTION 4: SPATIAL ANALYSIS BY SCENE
# ============================================================================

cat("\n========== SPATIAL ANALYSIS BY SCENE ==========\n\n")

# Scene-level summary
scene_summary <- merged_df[, .(
  n_cells = .N,
  n_epithelial = sum(CellType == "Epithelial"),
  n_immune = sum(CellType == "Immune"),
  n_fibroblast = sum(CellType == "Fibroblast"),
  mean_CK8 = mean(CK8_Ring, na.rm = TRUE),
  mean_ColI = mean(ColI_Ring, na.rm = TRUE),
  mean_aSMA = mean(aSMA_Ring, na.rm = TRUE),
  mean_CD45 = mean(CD45_Ring, na.rm = TRUE)
), by = scene]

scene_summary[, `:=`(
  pct_epithelial = 100 * n_epithelial / n_cells,
  pct_immune = 100 * n_immune / n_cells,
  pct_fibroblast = 100 * n_fibroblast / n_cells
)]

cat("\nScene-level summary:\n")
print(scene_summary)

fwrite(scene_summary, file.path(output_path, "PatientC_Scene_Summary.csv"))

# ============================================================================
# SECTION 5: SPATIAL PLOTS
# ============================================================================

cat("\n========== CREATING SPATIAL PLOTS ==========\n\n")

# Plot each scene
for (sc in unique(merged_df$scene)) {

  scene_data <- merged_df[scene == sc]
  cat("Plotting", sc, "with", nrow(scene_data), "cells\n")

  # Cell type spatial plot
  p_celltype <- ggplot(scene_data, aes(x = X_um, y = Y_um, color = CellType)) +
    geom_point(size = 0.3, alpha = 0.6) +
    scale_color_manual(values = c("Epithelial" = "#E41A1C", "Immune" = "#377EB8",
                                   "Fibroblast" = "#4DAF4A", "Other" = "grey70")) +
    coord_fixed() +
    theme_void() +
    labs(title = paste0("Patient C CycIF - ", sc, ": Cell Types"),
         color = "Cell Type") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "right")

  ggsave(file.path(output_path, paste0("PatientC_", sc, "_CellTypes.pdf")),
         p_celltype, width = 10, height = 8)

  # Key markers spatial plot
  p_ck8 <- ggplot(scene_data, aes(x = X_um, y = Y_um, color = log1p(CK8_Ring))) +
    geom_point(size = 0.2, alpha = 0.7) +
    scale_color_viridis(option = "magma", name = "CK8\n(log)") +
    coord_fixed() + theme_void() +
    labs(title = "CK8 (Epithelial)") +
    theme(plot.title = element_text(hjust = 0.5, size = 10))

  p_col <- ggplot(scene_data, aes(x = X_um, y = Y_um, color = log1p(ColI_Ring))) +
    geom_point(size = 0.2, alpha = 0.7) +
    scale_color_viridis(option = "viridis", name = "ColI\n(log)") +
    coord_fixed() + theme_void() +
    labs(title = "Collagen I (Stroma)") +
    theme(plot.title = element_text(hjust = 0.5, size = 10))

  p_asma <- ggplot(scene_data, aes(x = X_um, y = Y_um, color = log1p(aSMA_Ring))) +
    geom_point(size = 0.2, alpha = 0.7) +
    scale_color_viridis(option = "inferno", name = "aSMA\n(log)") +
    coord_fixed() + theme_void() +
    labs(title = "aSMA (mCAF)") +
    theme(plot.title = element_text(hjust = 0.5, size = 10))

  p_cd45 <- ggplot(scene_data, aes(x = X_um, y = Y_um, color = log1p(CD45_Ring))) +
    geom_point(size = 0.2, alpha = 0.7) +
    scale_color_viridis(option = "cividis", name = "CD45\n(log)") +
    coord_fixed() + theme_void() +
    labs(title = "CD45 (Immune)") +
    theme(plot.title = element_text(hjust = 0.5, size = 10))

  p_combined <- (p_ck8 + p_col) / (p_asma + p_cd45) +
    plot_annotation(
      title = paste0("Patient C CycIF - ", sc, ": Key Markers"),
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))
    )

  ggsave(file.path(output_path, paste0("PatientC_", sc, "_KeyMarkers.pdf")),
         p_combined, width = 14, height = 12)
}

# ============================================================================
# SECTION 6: CAF ZONE ANALYSIS (Cancer-Adjacent vs Cancer-Distal)
# ============================================================================

cat("\n========== CAF ZONE ANALYSIS ==========\n\n")

# For each fibroblast, calculate distance to nearest epithelial cell
# This helps define cancer-adjacent vs cancer-distal

epi_cells <- merged_df[CellType == "Epithelial", .(X_um, Y_um)]
fib_cells <- copy(fib_df)  # Use the already processed fib_df with CAF_subtype

if (nrow(epi_cells) > 0 && nrow(fib_cells) > 0) {

  cat("Calculating distances to nearest cancer cell...\n")

  # Sample epithelial cells for efficiency
  set.seed(42)
  epi_sample <- epi_cells[sample(.N, min(5000, .N))]

  # Calculate minimum distance for each fibroblast
  fib_cells[, dist_to_cancer := {
    min_dists <- numeric(.N)
    for (i in 1:.N) {
      dists <- sqrt((X_um[i] - epi_sample$X_um)^2 + (Y_um[i] - epi_sample$Y_um)^2)
      min_dists[i] <- min(dists)
    }
    min_dists
  }]

  # Define zones based on distance
  dist_threshold <- quantile(fib_cells$dist_to_cancer, 0.5, na.rm = TRUE)
  fib_cells[, Zone := ifelse(dist_to_cancer < dist_threshold, "Cancer-Adjacent", "Cancer-Distal")]

  cat("\nDistance threshold (median):", round(dist_threshold, 1), "um\n")
  cat("\nFibroblast zone distribution:\n")
  print(table(fib_cells$Zone))

  # Compare CAF markers by zone
  zone_markers <- fib_cells[, .(
    n_cells = .N,
    mean_ColI = mean(ColI_Ring, na.rm = TRUE),
    mean_ColIV = mean(ColIV_Ring, na.rm = TRUE),
    mean_aSMA = mean(aSMA_Ring, na.rm = TRUE),
    mean_Vim = mean(Vim_Ring, na.rm = TRUE),
    mean_PDPN = mean(PDPN_Ring, na.rm = TRUE),
    pct_mCAF = 100 * mean(CAF_subtype == "mCAF-like", na.rm = TRUE)
  ), by = Zone]

  cat("\nCAF markers by zone:\n")
  print(zone_markers)

  fwrite(zone_markers, file.path(output_path, "PatientC_CAF_ZoneMarkers.csv"))

  # Statistical test
  cat("\nWilcoxon tests (Adjacent vs Distal):\n")
  for (marker in c("ColI_Ring", "ColIV_Ring", "aSMA_Ring", "Vim_Ring")) {
    test <- wilcox.test(fib_cells[[marker]][fib_cells$Zone == "Cancer-Adjacent"],
                        fib_cells[[marker]][fib_cells$Zone == "Cancer-Distal"])
    cat(sprintf("  %s: p = %.2e\n", marker, test$p.value))
  }

  # Plot CAF markers by zone
  zone_long <- fib_cells %>%
    select(Zone, ColI_Ring, ColIV_Ring, aSMA_Ring, Vim_Ring, PDPN_Ring) %>%
    pivot_longer(cols = -Zone, names_to = "Marker", values_to = "Intensity") %>%
    mutate(Marker = gsub("_Ring", "", Marker))

  p_zone_box <- ggplot(zone_long, aes(x = Marker, y = log1p(Intensity), fill = Zone)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
    scale_fill_manual(values = c("Cancer-Adjacent" = "#9E503C", "Cancer-Distal" = "#2D5A5A")) +
    theme_classic(base_size = 12) +
    labs(title = "Patient C CycIF: CAF Markers by Spatial Zone",
         subtitle = "Zone defined by distance to nearest epithelial cell",
         x = NULL, y = "Log Intensity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")

  ggsave(file.path(output_path, "PatientC_CAF_ZoneComparison_Boxplot.pdf"),
         p_zone_box, width = 10, height = 6)

  # Spatial plot of zones
  for (sc in unique(fib_cells$scene)) {
    scene_fib <- fib_cells[scene == sc]

    p_zone_spatial <- ggplot(scene_fib, aes(x = X_um, y = Y_um, color = Zone)) +
      geom_point(size = 0.5, alpha = 0.7) +
      scale_color_manual(values = c("Cancer-Adjacent" = "#9E503C", "Cancer-Distal" = "#2D5A5A")) +
      coord_fixed() +
      theme_void() +
      labs(title = paste0("Patient C CycIF - ", sc, ": Fibroblast Zones"),
           color = "Zone") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            legend.position = "right")

    ggsave(file.path(output_path, paste0("PatientC_", sc, "_FibroblastZones.pdf")),
           p_zone_spatial, width = 10, height = 8)
  }
}

# ============================================================================
# SECTION 7: MARKER EXPRESSION HEATMAP
# ============================================================================

cat("\n========== MARKER HEATMAP ==========\n\n")

# Calculate mean expression by cell type
celltype_expr <- merged_df[, lapply(.SD, function(x) mean(as.numeric(x), na.rm = TRUE)),
                           by = CellType,
                           .SDcols = protein_markers[1:min(30, length(protein_markers))]]

# Prepare matrix
expr_matrix <- as.matrix(celltype_expr[, -1])
rownames(expr_matrix) <- celltype_expr$CellType

# Log transform and scale
expr_matrix <- log1p(expr_matrix)
expr_matrix <- t(scale(t(expr_matrix)))

# Remove NA columns
expr_matrix <- expr_matrix[, !apply(is.na(expr_matrix), 2, all)]

pdf(file.path(output_path, "PatientC_Marker_Heatmap_CellType.pdf"), width = 14, height = 5)
pheatmap(expr_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 12,
         fontsize_col = 8,
         main = "Patient C CycIF: Marker Expression by Cell Type")
dev.off()

cat("Heatmap saved\n")

# ============================================================================
# SECTION 8: VALIDATION OF CosMx CAF MARKERS
# ============================================================================

cat("\n========== VALIDATION OF CosMx CAF MARKERS ==========\n\n")

# CosMx identified these as zone-specific in fibroblasts:
# Cancer Adjacent (Active CAF): COL1A1, COL1A2, COL3A1, FN1, THBS1
# Cancer Distal (Quiescent CAF): APOD, NDRG1, NPPC

# CycIF equivalents:
# ColI_Ring = Collagen I (COL1A1/COL1A2)
# ColIV_Ring = Collagen IV
# aSMA_Ring = alpha-SMA (myofibroblast marker)
# Vim_Ring = Vimentin
# PDPN_Ring = Podoplanin

cat("Testing CosMx-identified zone markers in CycIF data:\n\n")

if (exists("fib_cells") && nrow(fib_cells) > 0) {

  # ColI should be higher in Cancer-Adjacent (matches COL1A1/COL1A2)
  adj_col <- mean(fib_cells$ColI_Ring[fib_cells$Zone == "Cancer-Adjacent"], na.rm = TRUE)
  dist_col <- mean(fib_cells$ColI_Ring[fib_cells$Zone == "Cancer-Distal"], na.rm = TRUE)

  cat(sprintf("ColI (COL1A1/2 proxy):\n"))
  cat(sprintf("  Cancer-Adjacent: %.1f\n", adj_col))
  cat(sprintf("  Cancer-Distal:   %.1f\n", dist_col))
  cat(sprintf("  Fold-change:     %.2f (Expected: >1 if Adjacent enriched)\n\n", adj_col/dist_col))

  # aSMA should be higher in Cancer-Adjacent (myofibroblast/mCAF)
  adj_asma <- mean(fib_cells$aSMA_Ring[fib_cells$Zone == "Cancer-Adjacent"], na.rm = TRUE)
  dist_asma <- mean(fib_cells$aSMA_Ring[fib_cells$Zone == "Cancer-Distal"], na.rm = TRUE)

  cat(sprintf("aSMA (mCAF marker):\n"))
  cat(sprintf("  Cancer-Adjacent: %.1f\n", adj_asma))
  cat(sprintf("  Cancer-Distal:   %.1f\n", dist_asma))
  cat(sprintf("  Fold-change:     %.2f (Expected: >1 if Adjacent enriched)\n\n", adj_asma/dist_asma))

  # Create validation barplot
  validation_df <- data.frame(
    Marker = c("ColI", "ColI", "aSMA", "aSMA", "ColIV", "ColIV", "Vim", "Vim"),
    Zone = rep(c("Cancer-Adjacent", "Cancer-Distal"), 4),
    Mean_Intensity = c(adj_col, dist_col, adj_asma, dist_asma,
                       mean(fib_cells$ColIV_Ring[fib_cells$Zone == "Cancer-Adjacent"], na.rm = TRUE),
                       mean(fib_cells$ColIV_Ring[fib_cells$Zone == "Cancer-Distal"], na.rm = TRUE),
                       mean(fib_cells$Vim_Ring[fib_cells$Zone == "Cancer-Adjacent"], na.rm = TRUE),
                       mean(fib_cells$Vim_Ring[fib_cells$Zone == "Cancer-Distal"], na.rm = TRUE))
  )

  p_validation <- ggplot(validation_df, aes(x = Marker, y = Mean_Intensity, fill = Zone)) +
    geom_col(position = position_dodge(0.8), width = 0.7, color = "black") +
    scale_fill_manual(values = c("Cancer-Adjacent" = "#9E503C", "Cancer-Distal" = "#2D5A5A")) +
    theme_classic(base_size = 12) +
    labs(title = "Patient C CycIF: CAF Marker Validation",
         subtitle = "Comparing zones defined by distance to epithelial cells",
         x = NULL, y = "Mean Intensity") +
    theme(legend.position = "top")

  ggsave(file.path(output_path, "PatientC_CAF_Validation_Barplot.pdf"),
         p_validation, width = 8, height = 5)
}

# ============================================================================
# SECTION 9: SUMMARY
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("PATIENT C CycIF ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n")

cat("\nTotal cells analyzed:", nrow(merged_df), "\n")
cat("Cell types identified:\n")
print(table(merged_df$CellType))

if (exists("fib_cells") && nrow(fib_cells) > 0) {
  cat("\nFibroblast zones:\n")
  print(table(fib_cells$Zone))
}

cat("\nOutput files saved to:", output_path, "\n")
cat("\nKey outputs:\n")
cat("  - PatientC_Scene_Summary.csv\n")
cat("  - PatientC_[scene]_CellTypes.pdf\n")
cat("  - PatientC_[scene]_KeyMarkers.pdf\n")
cat("  - PatientC_[scene]_FibroblastZones.pdf\n")
cat("  - PatientC_CAF_ZoneMarkers.csv\n")
cat("  - PatientC_CAF_ZoneComparison_Boxplot.pdf\n")
cat("  - PatientC_CAF_Validation_Barplot.pdf\n")
cat("  - PatientC_Marker_Heatmap_CellType.pdf\n")
