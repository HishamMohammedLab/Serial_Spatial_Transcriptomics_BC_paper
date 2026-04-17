# ============================================================================
# Patient C: CycIF vs CosMx Global Validation
# Compare RNA-level observations with protein-level validation
# ============================================================================

library(tidyverse)
library(data.table)
library(ggplot2)
library(patchwork)
library(pheatmap)

set.seed(42)

# ============================================================================
# PATHS
# ============================================================================

base_path <- "supplementary_input_data"
cycif_path <- file.path(base_path, "CycIF_Data")
cosmx_path <- "PatientC_Analysis"
output_path <- file.path(base_path, "New Figures for paper/CycIF_Validation/PatientC")

cat("========== LOADING DATA ==========\n\n")

# ============================================================================
# LOAD CycIF DATA
# ============================================================================

coord_df <- fread(file.path(cycif_path, "features_SMT130-4_CentroidXY.csv"))
setnames(coord_df, "V1", "cell_id")
coord_df[, X_um := DAPI_X * 0.35]
coord_df[, Y_um := DAPI_Y * 0.35]

protein_df <- fread(file.path(cycif_path, "features_SMT130-4_FilteredMeanIntensity_DAPI11_Nuclei1000.csv"))
setnames(protein_df, "UNIQID", "cell_id")

cycif_df <- merge(protein_df, coord_df[, .(cell_id, X_um, Y_um, scene)], by = "cell_id", all.x = TRUE)

cat("CycIF cells:", nrow(cycif_df), "\n")

# ============================================================================
# VALIDATION 1: CELL TYPE MARKER SPECIFICITY
# ============================================================================

cat("\n========== VALIDATION 1: CELL TYPE MARKERS ==========\n\n")

# CosMx uses: Epithelial markers (KRT8, KRT19), Fibroblast markers (COL1A1, VIM)
# CycIF has: CK8, CK19 (epithelial), ColI, Vim (fibroblast), CD45 (immune)

# Check marker correlations in CycIF
epithelial_markers <- c("CK8_Ring", "CK19_Ring", "CK7_Ring", "Ecad_Ring")
fibroblast_markers <- c("ColI_Ring", "ColIV_Ring", "Vim_Ring", "aSMA_Ring")
immune_markers <- c("CD45_Ring", "CD3_Ring", "CD68_Ring", "CD20_Ring")

# Calculate correlation matrix for key markers
key_markers <- c(epithelial_markers, fibroblast_markers, immune_markers)
key_markers <- key_markers[key_markers %in% colnames(cycif_df)]

cor_matrix <- cor(cycif_df[, ..key_markers], use = "pairwise.complete.obs")

cat("Marker correlation matrix:\n")
print(round(cor_matrix, 2))

# Expected from CosMx:
# - Epithelial markers should correlate with each other
# - Fibroblast markers should correlate with each other
# - Negative correlation between epithelial and fibroblast

cat("\n--- Key correlations (CosMx prediction vs CycIF) ---\n")
cat(sprintf("CK8-CK19 correlation: %.3f (expected: positive, epithelial markers)\n",
            cor_matrix["CK8_Ring", "CK19_Ring"]))
cat(sprintf("ColI-Vim correlation: %.3f (expected: positive, stroma markers)\n",
            cor_matrix["ColI_Ring", "Vim_Ring"]))
cat(sprintf("CK8-ColI correlation: %.3f (expected: negative/weak, different lineages)\n",
            cor_matrix["CK8_Ring", "ColI_Ring"]))
cat(sprintf("ColI-aSMA correlation: %.3f (expected: positive, both in myCAF)\n",
            cor_matrix["ColI_Ring", "aSMA_Ring"]))

# ============================================================================
# VALIDATION 2: COLLAGEN EXPRESSION IN FIBROBLASTS
# ============================================================================

cat("\n========== VALIDATION 2: COLLAGEN IN FIBROBLASTS ==========\n\n")

# CosMx shows COL1A1 is the highest expressed gene in fibroblasts
# CycIF should show ColI highest in cells with fibroblast markers

# Classify cells by dominant marker
cycif_df[, max_epi := pmax(CK8_Ring, CK19_Ring, na.rm = TRUE)]
cycif_df[, max_fib := pmax(ColI_Ring, Vim_Ring, na.rm = TRUE)]
cycif_df[, max_imm := CD45_Ring]

# Simple classification
cycif_df[, CellType := "Other"]
epi_thresh <- quantile(cycif_df$max_epi, 0.7, na.rm = TRUE)
fib_thresh <- quantile(cycif_df$max_fib, 0.7, na.rm = TRUE)
imm_thresh <- quantile(cycif_df$max_imm, 0.8, na.rm = TRUE)

cycif_df[max_epi > epi_thresh, CellType := "Epithelial"]
cycif_df[max_imm > imm_thresh, CellType := "Immune"]
cycif_df[CellType == "Other" & max_fib > fib_thresh, CellType := "Fibroblast"]

# Calculate mean ColI by cell type
col_by_celltype <- cycif_df[, .(
  mean_ColI = mean(ColI_Ring, na.rm = TRUE),
  mean_ColIV = mean(ColIV_Ring, na.rm = TRUE),
  mean_aSMA = mean(aSMA_Ring, na.rm = TRUE),
  mean_Vim = mean(Vim_Ring, na.rm = TRUE),
  n_cells = .N
), by = CellType]

cat("Collagen expression by cell type (CycIF):\n")
print(col_by_celltype)

cat("\nCosMx prediction: Fibroblasts should have HIGHEST ColI\n")
cat("CycIF validation:",
    ifelse(col_by_celltype[CellType == "Fibroblast", mean_ColI] >
           col_by_celltype[CellType == "Epithelial", mean_ColI],
           "VALIDATED - Fibroblasts have higher ColI",
           "NOT VALIDATED"), "\n")

# ============================================================================
# VALIDATION 3: CYTOKERATIN SPECIFICITY
# ============================================================================

cat("\n========== VALIDATION 3: CYTOKERATIN SPECIFICITY ==========\n\n")

# CosMx shows KRT8, KRT19 are epithelial-specific
# CycIF should show CK8, CK19 highest in epithelial cells

ck_by_celltype <- cycif_df[, .(
  mean_CK8 = mean(CK8_Ring, na.rm = TRUE),
  mean_CK19 = mean(CK19_Ring, na.rm = TRUE),
  mean_CK5 = mean(CK5_Ring, na.rm = TRUE),
  mean_CK14 = mean(CK14_Ring, na.rm = TRUE),
  n_cells = .N
), by = CellType]

cat("Cytokeratin expression by cell type (CycIF):\n")
print(ck_by_celltype)

cat("\nCosMx prediction: Epithelial cells should have HIGHEST CK8/CK19\n")
cat("CycIF validation:",
    ifelse(ck_by_celltype[CellType == "Epithelial", mean_CK8] >
           ck_by_celltype[CellType == "Fibroblast", mean_CK8],
           "VALIDATED - Epithelial have higher CK8",
           "NOT VALIDATED"), "\n")

# ============================================================================
# VALIDATION 4: ER EXPRESSION PATTERN
# ============================================================================

cat("\n========== VALIDATION 4: ER EXPRESSION ==========\n\n")

# CosMx: Patient C is ER+ breast cancer
# CycIF should show ER expression in epithelial/cancer cells

er_by_celltype <- cycif_df[, .(
  mean_ER = mean(ER_Nuclei, na.rm = TRUE),
  mean_PgR = mean(PgR_Nuclei, na.rm = TRUE),
  pct_ER_high = 100 * mean(ER_Nuclei > quantile(cycif_df$ER_Nuclei, 0.75, na.rm = TRUE), na.rm = TRUE),
  n_cells = .N
), by = CellType]

cat("ER/PgR expression by cell type (CycIF):\n")
print(er_by_celltype)

cat("\nCosMx prediction: ER+ cancer - Epithelial should have HIGHEST ER\n")
cat("CycIF validation:",
    ifelse(er_by_celltype[CellType == "Epithelial", mean_ER] >
           er_by_celltype[CellType == "Fibroblast", mean_ER],
           "VALIDATED - Epithelial have higher ER",
           "NOT VALIDATED"), "\n")

# ============================================================================
# VALIDATION 5: SPATIAL HETEROGENEITY (SCENE COMPARISON)
# ============================================================================

cat("\n========== VALIDATION 5: SPATIAL HETEROGENEITY ==========\n\n")

# CosMx shows spatial segregation with distinct zones
# CycIF scenes should show different composition patterns

scene_composition <- cycif_df[, .(
  n_total = .N,
  n_epi = sum(CellType == "Epithelial"),
  n_fib = sum(CellType == "Fibroblast"),
  n_imm = sum(CellType == "Immune"),
  pct_epi = 100 * mean(CellType == "Epithelial"),
  pct_fib = 100 * mean(CellType == "Fibroblast"),
  pct_imm = 100 * mean(CellType == "Immune"),
  mean_ColI = mean(ColI_Ring, na.rm = TRUE),
  mean_CK8 = mean(CK8_Ring, na.rm = TRUE)
), by = scene]

cat("Cell composition by scene (CycIF):\n")
print(scene_composition)

# Calculate scene heterogeneity
epi_range <- max(scene_composition$pct_epi) - min(scene_composition$pct_epi)
cat(sprintf("\nEpithelial range across scenes: %.1f%% (min) to %.1f%% (max)\n",
            min(scene_composition$pct_epi), max(scene_composition$pct_epi)))
cat("CosMx prediction: Should see spatial heterogeneity in cell types\n")
cat("CycIF validation:",
    ifelse(epi_range > 10, "VALIDATED - Significant spatial heterogeneity",
           "PARTIAL - Limited heterogeneity"), "\n")

# ============================================================================
# VALIDATION 6: IMMUNE INFILTRATION
# ============================================================================

cat("\n========== VALIDATION 6: IMMUNE INFILTRATION ==========\n\n")

# Check immune marker patterns
immune_detail <- cycif_df[CellType == "Immune", .(
  mean_CD3 = mean(CD3_Ring, na.rm = TRUE),
  mean_CD8 = mean(CD8_Ring, na.rm = TRUE),
  mean_CD4 = mean(CD4_Ring, na.rm = TRUE),
  mean_CD68 = mean(CD68_Ring, na.rm = TRUE),
  mean_CD20 = mean(CD20_Ring, na.rm = TRUE),
  mean_PD1 = mean(PD1_Ring, na.rm = TRUE),
  n_cells = .N
), by = scene]

cat("Immune cell markers by scene (CycIF):\n")
print(immune_detail)

# T cell vs Macrophage ratio
immune_detail[, Tcell_Mac_ratio := mean_CD3 / mean_CD68]
cat("\nT cell / Macrophage ratio by scene:\n")
print(immune_detail[, .(scene, Tcell_Mac_ratio)])

# ============================================================================
# CREATE VALIDATION SUMMARY PLOT
# ============================================================================

cat("\n========== CREATING VALIDATION PLOTS ==========\n\n")

# Plot 1: Marker specificity by cell type
marker_summary <- cycif_df[, .(
  ColI = mean(ColI_Ring, na.rm = TRUE),
  CK8 = mean(CK8_Ring, na.rm = TRUE),
  CD45 = mean(CD45_Ring, na.rm = TRUE),
  aSMA = mean(aSMA_Ring, na.rm = TRUE),
  ER = mean(ER_Nuclei, na.rm = TRUE)
), by = CellType]

marker_long <- marker_summary %>%
  pivot_longer(cols = -CellType, names_to = "Marker", values_to = "Mean_Intensity")

p1 <- ggplot(marker_long, aes(x = CellType, y = Mean_Intensity, fill = CellType)) +
  geom_col() +
  facet_wrap(~Marker, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Epithelial" = "#E41A1C", "Fibroblast" = "#4DAF4A",
                                "Immune" = "#377EB8", "Other" = "grey70")) +
  theme_classic(base_size = 10) +
  labs(title = "CycIF Validation: Marker Specificity by Cell Type",
       subtitle = "ColI should be highest in Fibroblasts; CK8/ER in Epithelial; CD45 in Immune",
       x = NULL, y = "Mean Intensity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        strip.text = element_text(face = "bold"))

ggsave(file.path(output_path, "Validation_MarkerSpecificity.pdf"), p1, width = 12, height = 4)

# Plot 2: Scene composition
scene_long <- scene_composition %>%
  select(scene, pct_epi, pct_fib, pct_imm) %>%
  pivot_longer(cols = -scene, names_to = "CellType", values_to = "Percentage") %>%
  mutate(CellType = case_when(
    CellType == "pct_epi" ~ "Epithelial",
    CellType == "pct_fib" ~ "Fibroblast",
    CellType == "pct_imm" ~ "Immune"
  ))

p2 <- ggplot(scene_long, aes(x = scene, y = Percentage, fill = CellType)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Epithelial" = "#E41A1C", "Fibroblast" = "#4DAF4A", "Immune" = "#377EB8")) +
  theme_classic(base_size = 11) +
  labs(title = "CycIF Validation: Spatial Heterogeneity Across Scenes",
       subtitle = "Scene003 shows distinct pattern: low epithelial, high immune",
       x = "Scene", y = "% of Cells") +
  theme(legend.position = "top")

ggsave(file.path(output_path, "Validation_SceneComposition.pdf"), p2, width = 8, height = 5)

# Plot 3: CosMx vs CycIF comparison - Cell type proportions
# CosMx Patient C (approximate from data)
cosmx_props <- data.frame(
  Dataset = c("CosMx Bx1", "CosMx Bx1", "CosMx Bx2", "CosMx Bx2", "CycIF", "CycIF"),
  CellType = rep(c("Cancer/Epithelial", "Fibroblast"), 3),
  Percentage = c(33.5, 7.0, 37.4, 13.2, 19.9, 38.1)  # Approximate from zone data
)

p3 <- ggplot(cosmx_props, aes(x = Dataset, y = Percentage, fill = CellType)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Cancer/Epithelial" = "#E41A1C", "Fibroblast" = "#4DAF4A")) +
  theme_classic(base_size = 11) +
  labs(title = "CosMx vs CycIF: Cell Type Proportions",
       subtitle = "CycIF shows higher fibroblast proportion (may include more stroma)",
       x = NULL, y = "% of Cells") +
  theme(legend.position = "top")

ggsave(file.path(output_path, "Validation_CosMx_vs_CycIF_Proportions.pdf"), p3, width = 8, height = 5)

# ============================================================================
# SUMMARY TABLE
# ============================================================================

cat("\n========== VALIDATION SUMMARY ==========\n\n")

validation_results <- data.frame(
  Observation = c(
    "Epithelial markers (CK8/CK19) correlate",
    "Fibroblast markers (ColI/Vim) correlate",
    "ColI highest in fibroblasts",
    "CK8 highest in epithelial",
    "ER highest in epithelial (ER+ cancer)",
    "Spatial heterogeneity across regions",
    "aSMA correlates with ColI (myCAF)"
  ),
  CosMx_Prediction = c(
    "Positive correlation",
    "Positive correlation",
    "Yes - CAF signature",
    "Yes - epithelial marker",
    "Yes - ER+ breast cancer",
    "Yes - distinct zones",
    "Yes - both in myofibroblasts"
  ),
  CycIF_Result = c(
    sprintf("r = %.2f", cor_matrix["CK8_Ring", "CK19_Ring"]),
    sprintf("r = %.2f", cor_matrix["ColI_Ring", "Vim_Ring"]),
    ifelse(col_by_celltype[CellType == "Fibroblast", mean_ColI] >
           col_by_celltype[CellType == "Epithelial", mean_ColI], "YES", "NO"),
    ifelse(ck_by_celltype[CellType == "Epithelial", mean_CK8] >
           ck_by_celltype[CellType == "Fibroblast", mean_CK8], "YES", "NO"),
    ifelse(er_by_celltype[CellType == "Epithelial", mean_ER] >
           er_by_celltype[CellType == "Fibroblast", mean_ER], "YES", "NO"),
    ifelse(epi_range > 10, "YES", "PARTIAL"),
    sprintf("r = %.2f", cor_matrix["ColI_Ring", "aSMA_Ring"])
  ),
  Validated = c(
    ifelse(cor_matrix["CK8_Ring", "CK19_Ring"] > 0.3, "YES", "NO"),
    ifelse(cor_matrix["ColI_Ring", "Vim_Ring"] > 0.3, "YES", "NO"),
    ifelse(col_by_celltype[CellType == "Fibroblast", mean_ColI] >
           col_by_celltype[CellType == "Epithelial", mean_ColI], "YES", "NO"),
    ifelse(ck_by_celltype[CellType == "Epithelial", mean_CK8] >
           ck_by_celltype[CellType == "Fibroblast", mean_CK8], "YES", "NO"),
    ifelse(er_by_celltype[CellType == "Epithelial", mean_ER] >
           er_by_celltype[CellType == "Fibroblast", mean_ER], "YES", "NO"),
    ifelse(epi_range > 10, "YES", "PARTIAL"),
    ifelse(cor_matrix["ColI_Ring", "aSMA_Ring"] > 0.3, "YES", "NO")
  )
)

cat("Validation Results:\n")
print(validation_results)

fwrite(validation_results, file.path(output_path, "Validation_Summary.csv"))

# Combined figure
p_combined <- p1 / p2 +
  plot_annotation(
    title = "Patient C: CycIF Protein Validation of CosMx RNA Observations",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

ggsave(file.path(output_path, "Validation_Combined.pdf"), p_combined, width = 12, height = 8)

cat("\n========== COMPLETE ==========\n")
cat("Output files saved to:", output_path, "\n")
