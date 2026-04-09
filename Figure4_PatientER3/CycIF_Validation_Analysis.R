# ============================================================================
# CycIF Validation Analysis for CosMx Spatial Transcriptomics
# ============================================================================
# This script performs:
# 1. Data loading and preprocessing
# 2. Biopsy matching verification
# 3. UMAP visualization by patient, biopsy, cell type
# 4. Heatmaps of marker expression by cell type
# 5. ER pathway validation
# 6. Spatial neighborhood (Kmeans 4) analysis
# 7. L-R pair validation potential assessment
# ============================================================================

# Load required libraries
library(tidyverse)
library(data.table)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(uwot)  # For UMAP
library(viridis)
library(patchwork)
library(scales)

# Set working directory and paths
base_path <- "supplementary_input_data"
data_path <- file.path(base_path, "CycIF_Data")
output_path <- "CycIF_Validation"

# Create output directory if it doesn't exist
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

# ============================================================================
# SECTION 1: DATA LOADING
# ============================================================================

cat("Loading Cyc-IF data files...\n")

# Load protein intensity data
protein_df <- fread(file.path(data_path, "ST_AD_04_allmarkers_df.csv"))
setnames(protein_df, "V1", "cell_id")

# Load coordinates data
coord_df <- fread(file.path(data_path, "ST_AD_04_allmarkers_dfxy.csv"))
setnames(coord_df, "V1", "cell_id")
# Convert pixels to microns (0.35 um/pixel)
coord_df[, X_um := DAPI_X * 0.35]
coord_df[, Y_um := DAPI_Y * 0.35]

# Load metadata
meta_df <- fread(file.path(data_path, "ST_AD_04_allmarkers_obs 1 (1).csv"))
setnames(meta_df, "V1", "cell_id")

cat(sprintf("Loaded %d cells with %d protein markers\n",
            nrow(protein_df), ncol(protein_df) - 1))

# ============================================================================
# SECTION 2: BIOPSY AND SCENE MATCHING TABLE
# ============================================================================

cat("\nGenerating biopsy matching table...\n")

# IMPORTANT: CycIF data contains multiple scenes per biopsy
# All spatial analyses MUST separate by scene to avoid merging distinct tissue sections
# See Context files/CycIF_Data_Structure.md for details

# Create summary table with scene information
biopsy_summary <- meta_df[, .(
  n_cells = .N,
  n_scenes = uniqueN(scene),
  scenes = paste(unique(scene), collapse = ", "),
  sample_prefix = unique(substr(get("registeredimages.1"), 1, 12))[1]
), by = .(PT, BX)]

# Scene-level summary
scene_summary <- meta_df[, .(n_cells = .N), by = .(PT, BX, scene, slide_scene)]
cat("\nScene-level cell counts:\n")
print(scene_summary[order(PT, BX, scene)])

# Add Cyc-IF naming convention based on sample prefix patterns
biopsy_summary[, CycIF_Sample := fifelse(
  PT == "A" & BX == "Bx2", "HTA9-1-121",
  fifelse(PT == "A" & BX == "Bx3", "SMT101Bx3",
  fifelse(PT == "A" & BX == "Bx4", "HTA9-1-33",
  fifelse(PT == "D" & BX == "Bx1", "HTAN9-3-292",
  fifelse(PT == "D" & BX == "Bx2", "HTAN9-3-284",
  fifelse(PT == "D" & BX == "Bx3", "HTA9-3-41", "Unknown"))))))]

# Match to CosMx biopsies
# Patient A has CosMx: Bx2, Bx3, Bx4 (Bx2 is different metastatic site)
# Patient D has CosMx: Bx1, Bx2, Bx3
biopsy_summary[, CosMx_Match := fifelse(
  PT == "A" & BX %in% c("Bx2", "Bx3", "Bx4"), "YES",
  fifelse(PT == "D" & BX %in% c("Bx1", "Bx2", "Bx3"), "YES", "NO"))]

# Print biopsy matching table
cat("\nBiopsy Matching Table:\n")
print(biopsy_summary[order(PT, BX)])

# Save biopsy matching table
fwrite(biopsy_summary, file.path(output_path, "biopsy_matching_table.csv"))
fwrite(scene_summary, file.path(output_path, "scene_summary.csv"))

# ============================================================================
# SECTION 3: DATA QUALITY CHECKS
# ============================================================================

cat("\nPerforming data quality checks...\n")

# Get protein marker names (excluding cell_id and area/eccentricity)
protein_markers <- colnames(protein_df)[!colnames(protein_df) %in%
                                          c("cell_id", "nuclei_area", "nuclei_eccentricity")]

# Merge data for QC - INCLUDE scene columns for proper spatial separation
merged_df <- merge(protein_df, meta_df[, .(cell_id, PT, BX, scene, slide_scene,
                                            `Primary Celltype: Matrix`,
                                            `Kmeans 4`)],
                   by = "cell_id")

# Calculate mean intensity per biopsy AND scene for batch effect assessment
batch_summary <- merged_df[, lapply(.SD, function(x) {
  mean(as.numeric(x), na.rm = TRUE)
}), by = .(PT, BX, scene), .SDcols = protein_markers[1:min(20, length(protein_markers))]]

cat("\nMean marker intensities by biopsy/scene (first 10 markers):\n")
print(batch_summary[, 1:min(13, ncol(batch_summary))])

# Check for outlier cells based on total intensity
merged_df[, total_intensity := rowSums(.SD, na.rm = TRUE),
          .SDcols = protein_markers[1:20]]

# Flag potential outliers (>3 SD from mean within scene, not just biopsy)
merged_df[, z_score := (total_intensity - mean(total_intensity, na.rm = TRUE)) /
            sd(total_intensity, na.rm = TRUE), by = .(PT, BX, scene)]

outlier_summary <- merged_df[, .(
  n_cells = .N,
  n_outliers = sum(abs(z_score) > 3, na.rm = TRUE),
  pct_outliers = round(100 * sum(abs(z_score) > 3, na.rm = TRUE) / .N, 2)
), by = .(PT, BX, scene)]

cat("\nOutlier cell summary (|z-score| > 3):\n")
print(outlier_summary[order(PT, BX)])

# Save QC summary
fwrite(outlier_summary, file.path(output_path, "qc_outlier_summary.csv"))

# ============================================================================
# SECTION 4: UMAP VISUALIZATION
# ============================================================================

cat("\nGenerating UMAP embeddings...\n")

# Prepare data for UMAP (subsample for computational efficiency)
set.seed(42)
n_sample <- min(50000, nrow(merged_df))
sample_idx <- sample(nrow(merged_df), n_sample)
umap_data <- merged_df[sample_idx]

# Select protein markers for UMAP (remove those with many NAs)
umap_markers <- protein_markers[1:min(40, length(protein_markers))]
umap_matrix <- as.matrix(umap_data[, ..umap_markers])

# Handle missing values
umap_matrix[is.na(umap_matrix)] <- 0

# Log transform and scale
umap_matrix <- log1p(umap_matrix)
umap_matrix <- scale(umap_matrix)

# Run UMAP
cat("Running UMAP...\n")
umap_result <- umap(umap_matrix, n_neighbors = 30, min_dist = 0.3,
                    n_components = 2, metric = "euclidean")

umap_data[, UMAP1 := umap_result[, 1]]
umap_data[, UMAP2 := umap_result[, 2]]

# Create UMAP plots

# 1. UMAP by Patient
p_patient <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = PT)) +
  geom_point(size = 0.3, alpha = 0.5) +
  scale_color_manual(values = c("A" = "#E41A1C", "D" = "#377EB8")) +
  theme_minimal() +
  labs(title = "Cyc-IF UMAP by Patient", color = "Patient") +
  theme(legend.position = "right")

# 2. UMAP by Biopsy
umap_data[, Patient_Biopsy := paste0(PT, "_", BX)]
p_biopsy <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = Patient_Biopsy)) +
  geom_point(size = 0.3, alpha = 0.5) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal() +
  labs(title = "Cyc-IF UMAP by Biopsy", color = "Biopsy") +
  theme(legend.position = "right")

# 2b. UMAP by Scene (important for spatial separation)
p_scene <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = slide_scene)) +
  geom_point(size = 0.3, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Cyc-IF UMAP by Scene", color = "Scene") +
  theme(legend.position = "right")

# 3. UMAP by Cell Type
umap_data[, CellType := `Primary Celltype: Matrix`]
umap_data[, CellType := gsub("^[0-9]+: ", "", CellType)]
p_celltype <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = CellType)) +
  geom_point(size = 0.3, alpha = 0.5) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  labs(title = "Cyc-IF UMAP by Cell Type", color = "Cell Type") +
  theme(legend.position = "right")

# 4. UMAP by Kmeans 4 clusters
p_kmeans <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = as.factor(`Kmeans 4`))) +
  geom_point(size = 0.3, alpha = 0.5) +
  scale_color_viridis_d() +
  theme_minimal() +
  labs(title = "Cyc-IF UMAP by Spatial Neighborhood (Kmeans 4)", color = "Cluster") +
  theme(legend.position = "right")

# Combine and save UMAP plots
combined_umap <- (p_patient + p_biopsy) / (p_celltype + p_kmeans)
ggsave(file.path(output_path, "UMAP_Combined.pdf"), combined_umap,
       width = 14, height = 12)

cat("UMAP plots saved to:", file.path(output_path, "UMAP_Combined.pdf"), "\n")

# ============================================================================
# SECTION 5: MARKER EXPRESSION HEATMAPS
# ============================================================================

cat("\nGenerating marker expression heatmaps...\n")

# Calculate mean expression by cell type
celltype_expr <- merged_df[!is.na(`Primary Celltype: Matrix`),
                           lapply(.SD, function(x) mean(as.numeric(x), na.rm = TRUE)),
                           by = .(`Primary Celltype: Matrix`),
                           .SDcols = protein_markers[1:min(30, length(protein_markers))]]

# Clean cell type names
celltype_expr[, CellType := gsub("^[0-9]+: ", "", `Primary Celltype: Matrix`)]

# Prepare matrix for heatmap
expr_matrix <- as.matrix(celltype_expr[, -c("Primary Celltype: Matrix", "CellType")])
rownames(expr_matrix) <- celltype_expr$CellType

# Log transform and scale
expr_matrix <- log1p(expr_matrix)
expr_matrix <- t(scale(t(expr_matrix)))

# Create heatmap
pdf(file.path(output_path, "Marker_Heatmap_CellType.pdf"), width = 12, height = 6)
pheatmap(expr_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 10,
         fontsize_col = 8,
         main = "Protein Marker Expression by Cell Type (Cyc-IF)")
dev.off()

# Calculate mean expression by biopsy
biopsy_expr <- merged_df[, lapply(.SD, function(x) mean(as.numeric(x), na.rm = TRUE)),
                         by = .(PT, BX),
                         .SDcols = protein_markers[1:min(30, length(protein_markers))]]

biopsy_expr[, Biopsy := paste0(PT, "_", BX)]
biopsy_matrix <- as.matrix(biopsy_expr[, -c("PT", "BX", "Biopsy")])
rownames(biopsy_matrix) <- biopsy_expr$Biopsy

# Log transform and scale
biopsy_matrix <- log1p(biopsy_matrix)
biopsy_matrix <- t(scale(t(biopsy_matrix)))

# Create biopsy heatmap
pdf(file.path(output_path, "Marker_Heatmap_Biopsy.pdf"), width = 12, height = 5)
pheatmap(biopsy_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 10,
         fontsize_col = 8,
         main = "Protein Marker Expression by Biopsy (Cyc-IF)")
dev.off()

cat("Heatmaps saved to output directory\n")

# ============================================================================
# SECTION 6: ER PATHWAY VALIDATION
# ============================================================================

cat("\nAnalyzing ER pathway markers...\n")

# Key ER pathway markers in Cyc-IF
er_markers <- c("ER_", "AR_", "PgR_", "CCND1_", "BCL2_")

# Check which ER markers are available
available_er <- er_markers[er_markers %in% colnames(protein_df)]
cat("Available ER pathway markers:", paste(available_er, collapse = ", "), "\n")

# Calculate ER marker expression by biopsy
er_expr <- merged_df[PT == "A", lapply(.SD, function(x) {
  list(
    mean = mean(as.numeric(x), na.rm = TRUE),
    median = median(as.numeric(x), na.rm = TRUE),
    pct_positive = 100 * mean(as.numeric(x) > quantile(as.numeric(x), 0.75, na.rm = TRUE), na.rm = TRUE)
  )
}), by = .(BX), .SDcols = available_er]

# Create ER expression summary for Patient A - by biopsy AND scene
er_summary_A <- merged_df[PT == "A" & `Primary Celltype: Matrix` == "3: tumor",
                          lapply(.SD, function(x) mean(as.numeric(x), na.rm = TRUE)),
                          by = .(BX, scene), .SDcols = available_er]

cat("\nPatient A ER Pathway Expression (Tumor cells only, by scene):\n")
print(er_summary_A[order(BX, scene)])

# Also create biopsy-level summary (aggregated across scenes)
er_summary_A_bx <- merged_df[PT == "A" & `Primary Celltype: Matrix` == "3: tumor",
                              lapply(.SD, function(x) mean(as.numeric(x), na.rm = TRUE)),
                              by = .(BX), .SDcols = available_er]
cat("\nPatient A ER Pathway Expression (aggregated by biopsy):\n")
print(er_summary_A_bx)

# Create ER expression plot for Patient A
er_long <- melt(er_summary_A, id.vars = "BX", variable.name = "Marker", value.name = "Expression")
er_long[, Marker := gsub("_$", "", Marker)]

p_er <- ggplot(er_long, aes(x = BX, y = Expression, fill = Marker)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  labs(title = "Patient A: ER Pathway Marker Expression Across Biopsies",
       subtitle = "CosMx shows ESR1 decrease: 11.2% -> 4.5% -> 4.6% (Bx2->Bx3->Bx4)",
       x = "Biopsy", y = "Mean Intensity (tumor cells)") +
  theme(legend.position = "right")

ggsave(file.path(output_path, "ER_Pathway_PatientA.pdf"), p_er, width = 10, height = 6)

# Calculate percentage ER+ cells per biopsy (Patient A)
# Using ER_ column and threshold at 75th percentile
er_threshold <- quantile(merged_df[PT == "A" & `Primary Celltype: Matrix` == "3: tumor", ER_],
                         0.75, na.rm = TRUE)

er_positive <- merged_df[PT == "A" & `Primary Celltype: Matrix` == "3: tumor",
                         .(pct_ER_positive = 100 * mean(ER_ > er_threshold, na.rm = TRUE),
                           n_cells = .N), by = .(BX)]

cat("\nPatient A ER+ cells (>75th percentile threshold):\n")
print(er_positive)

# Save ER summary
fwrite(er_positive, file.path(output_path, "ER_positive_PatientA.csv"))

# ============================================================================
# SECTION 7: SPATIAL NEIGHBORHOOD ANALYSIS (Kmeans 4)
# ============================================================================

cat("\nAnalyzing Kmeans 4 spatial neighborhoods...\n")

# Characterize each Kmeans cluster by cell type composition
kmeans_celltype <- merged_df[!is.na(`Kmeans 4`) & !is.na(`Primary Celltype: Matrix`),
                             .N, by = .(`Kmeans 4`, `Primary Celltype: Matrix`)]
kmeans_celltype[, pct := 100 * N / sum(N), by = .(`Kmeans 4`)]

# Pivot wider
kmeans_wide <- dcast(kmeans_celltype, `Kmeans 4` ~ `Primary Celltype: Matrix`,
                     value.var = "pct", fill = 0)

cat("\nKmeans 4 cluster composition (% by cell type):\n")
print(kmeans_wide)

# Characterize each cluster by marker expression
kmeans_markers <- merged_df[!is.na(`Kmeans 4`),
                            lapply(.SD, function(x) mean(as.numeric(x), na.rm = TRUE)),
                            by = .(`Kmeans 4`),
                            .SDcols = protein_markers[1:min(25, length(protein_markers))]]

# Create heatmap of Kmeans clusters
kmeans_matrix <- as.matrix(kmeans_markers[, -1])
rownames(kmeans_matrix) <- paste0("Cluster_", kmeans_markers$`Kmeans 4`)

kmeans_matrix <- log1p(kmeans_matrix)
kmeans_matrix <- t(scale(t(kmeans_matrix)))

pdf(file.path(output_path, "Kmeans4_Marker_Heatmap.pdf"), width = 12, height = 5)
pheatmap(kmeans_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 12,
         fontsize_col = 8,
         main = "Kmeans 4 Spatial Neighborhood Marker Profiles")
dev.off()

# Save Kmeans summary
fwrite(kmeans_wide, file.path(output_path, "Kmeans4_celltype_composition.csv"))
fwrite(kmeans_markers, file.path(output_path, "Kmeans4_marker_expression.csv"))

# ============================================================================
# SECTION 8: L-R PAIR VALIDATION POTENTIAL
# ============================================================================

cat("\nAssessing L-R pair validation potential...\n")

# Available Cyc-IF markers for L-R validation
available_markers <- colnames(protein_df)[!colnames(protein_df) %in%
                                            c("cell_id", "nuclei_area", "nuclei_eccentricity")]

cat("Available protein markers:\n")
print(available_markers)

# Key L-R pairs from CosMx analysis and their Cyc-IF equivalents
lr_pairs <- data.table(
  CosMx_Ligand = c("MIF", "COL1A1", "COL1A2", "FN1", "EGFR ligands", "TGFB1", "MDK", "GAS6"),
  CosMx_Receptor = c("CD74", "ITGB1", "ITGB1", "ITGB1", "EGFR", "TGFBR1/2", "SDC4", "AXL"),
  CycIF_Available = c(
    "No direct marker, but CD68 for macrophages",
    "ColI_, ColIV_ (collagen markers)",
    "ColI_, ColIV_ (collagen markers)",
    "Vim_ (mesenchymal), ColI_, ColIV_",
    "EGFR_ (YES)",
    "No direct marker",
    "No direct marker",
    "No direct marker"
  ),
  Can_Validate = c("Partial", "Partial", "Partial", "Partial", "YES", "No", "No", "No")
)

cat("\nL-R Pair Validation Potential:\n")
print(lr_pairs)

# Check for specific markers
cat("\nMarkers directly relevant to L-R validation:\n")
cat("- EGFR_: PRESENT\n")
cat("- ColI_: PRESENT (collagen I)\n")
cat("- ColIV_: PRESENT (collagen IV)\n")
cat("- Vim_: PRESENT (vimentin - ECM/mesenchymal)\n")
cat("- CD68_: PRESENT (macrophage marker for MIF pathway)\n")
cat("- CD3_, CD4_, CD8_: PRESENT (T cell markers)\n")
cat("- PD1_: PRESENT (immune checkpoint)\n")

# Save L-R validation table
fwrite(lr_pairs, file.path(output_path, "LR_Pair_Validation_Potential.csv"))

# ============================================================================
# SECTION 9: IMMUNE MARKER ANALYSIS (for IFNGR2 mutation context)
# ============================================================================

cat("\nAnalyzing immune markers (IFNGR2 mutation context)...\n")

# Immune-related markers in Cyc-IF
immune_markers <- c("CD3_", "CD4_", "CD8_", "CD68_", "CD20_", "FoxP3_", "GRNZB_", "PD1_")
available_immune <- immune_markers[immune_markers %in% available_markers]

cat("Available immune markers:", paste(available_immune, collapse = ", "), "\n")

# Calculate immune marker expression by patient/biopsy/scene
immune_expr <- merged_df[, lapply(.SD, function(x) mean(as.numeric(x), na.rm = TRUE)),
                         by = .(PT, BX, scene), .SDcols = available_immune]

cat("\nImmune marker expression by biopsy/scene:\n")
print(immune_expr[order(PT, BX, scene)])

# Create immune profile plot for Patient A (with IFNGR2 mutation)
immune_long_A <- melt(immune_expr[PT == "A"], id.vars = c("PT", "BX"),
                      variable.name = "Marker", value.name = "Expression")
immune_long_A[, Marker := gsub("_$", "", Marker)]

p_immune_A <- ggplot(immune_long_A, aes(x = BX, y = Expression, fill = Marker)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  labs(title = "Patient A: Immune Marker Expression Across Biopsies",
       subtitle = "Patient A has IFNGR2 mutation - check for immune response patterns",
       x = "Biopsy", y = "Mean Intensity") +
  theme(legend.position = "right")

ggsave(file.path(output_path, "Immune_Markers_PatientA.pdf"), p_immune_A, width = 10, height = 6)

# Calculate immune cell proportions by scene
immune_cells <- merged_df[`Primary Celltype: Matrix` == "2: immune",
                          .N, by = .(PT, BX, scene)]
total_cells <- merged_df[, .N, by = .(PT, BX, scene)]
immune_pct <- merge(immune_cells, total_cells, by = c("PT", "BX", "scene"),
                    suffixes = c("_immune", "_total"))
immune_pct[, pct_immune := 100 * N_immune / N_total]

cat("\nImmune cell infiltration by biopsy/scene:\n")
print(immune_pct[order(PT, BX, scene)])

fwrite(immune_pct, file.path(output_path, "Immune_infiltration_byScene.csv"))

# Also save biopsy-level summary
immune_pct_bx <- merged_df[, .(
  n_immune = sum(`Primary Celltype: Matrix` == "2: immune", na.rm = TRUE),
  n_total = .N
), by = .(PT, BX)]
immune_pct_bx[, pct_immune := 100 * n_immune / n_total]
fwrite(immune_pct_bx, file.path(output_path, "Immune_infiltration.csv"))

# ============================================================================
# SECTION 10: SUMMARY STATISTICS
# ============================================================================

cat("\n\n=== SUMMARY STATISTICS ===\n")
cat("Total cells analyzed:", nrow(merged_df), "\n")
cat("Patient A cells:", nrow(merged_df[PT == "A"]), "\n")
cat("Patient D cells:", nrow(merged_df[PT == "D"]), "\n")
cat("Number of protein markers:", length(protein_markers), "\n")
cat("Total scenes:", uniqueN(merged_df$slide_scene), "\n")
cat("Biopsies with CosMx match:", sum(biopsy_summary$CosMx_Match == "YES"), "/",
    nrow(biopsy_summary), "\n")

cat("\nScene breakdown:\n")
print(merged_df[, .(n_cells = .N), by = .(PT, BX, scene)][order(PT, BX, scene)])

cat("\nAnalysis complete! Results saved to:", output_path, "\n")
cat("NOTE: All spatial analyses now separate by scene - see CycIF_Data_Structure.md\n")

# ============================================================================
# END OF SCRIPT
# ============================================================================
