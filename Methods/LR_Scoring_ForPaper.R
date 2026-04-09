#!/usr/bin/env Rscript
# =============================================================================
# Ligand-Receptor Interaction Scoring Across All Patients
# =============================================================================
#
# Computes L-R interaction scores for all patient-timepoint combinations
# using cell type-specific mean expression from CosMx spatial transcriptomics.
#
# Score = Mean_Ligand_Expr(sender) x Mean_Receptor_Expr(receiver)
#
# Input:
#   - Seurat object with Lineage annotations (SMMART_with_Lineage.rds)
#   - Curated L-R pair database (CosMx_Panel_LR_Pairs.csv)
#
# Output:
#   - LR_Scores_AllPatients_AllTimepoints.csv (per-pair, per-timepoint scores)
#   - LR_FoldChanges_PerPatient.csv (first vs last biopsy fold changes)
#   - LR_ConsistentChanges_Summary.csv (cross-patient consistency)
#
# Author: Hisham Mohammed
# =============================================================================

library(Seurat)
library(tidyverse)

# =============================================================================
# 1. Load Data
# =============================================================================

cat("Loading Seurat object...\n")
seurat_obj <- readRDS("data/CosMx_SMMART_345k_clean.rds")

lr_pairs <- read.csv("supplementary_input_data/LR_Scores/CosMx_Panel_LR_Pairs.csv",
                      stringsAsFactors = FALSE)

# Filter to pairs with both genes present
available_genes <- rownames(seurat_obj)
lr_pairs <- lr_pairs %>%
  filter(Ligand %in% available_genes & Receptor %in% available_genes) %>%
  distinct(Ligand, Receptor, .keep_all = TRUE)

cat(sprintf("  %d L-R pairs with both genes in 960-gene panel\n", nrow(lr_pairs)))

# =============================================================================
# 2. Pre-compute Expression Means by Group
# =============================================================================

cat("\nPre-computing expression means by Patient-Timepoint-Lineage...\n")

meta <- seurat_obj@meta.data %>%
  mutate(cell_id = rownames(seurat_obj@meta.data))
expr_matrix <- GetAssayData(seurat_obj, layer = "data")

genes_needed <- unique(c(lr_pairs$Ligand, lr_pairs$Receptor))
genes_needed <- genes_needed[genes_needed %in% rownames(expr_matrix)]

meta$col_idx <- match(meta$cell_id, colnames(expr_matrix))

groups <- meta %>%
  group_by(Patient, Timepoint, Lineage) %>%
  summarise(cell_indices = list(col_idx), n_cells = n(), .groups = "drop") %>%
  filter(n_cells >= 10)

# Compute mean expression and % expressing for all genes in each group
expr_means_list <- list()
pb <- txtProgressBar(min = 0, max = nrow(groups), style = 3)
for (i in 1:nrow(groups)) {
  idx <- groups$cell_indices[[i]]
  means <- rowMeans(expr_matrix[genes_needed, idx, drop = FALSE])
  pcts <- rowMeans(expr_matrix[genes_needed, idx, drop = FALSE] > 0) * 100

  expr_means_list[[i]] <- data.frame(
    Patient = groups$Patient[i],
    Timepoint = groups$Timepoint[i],
    Lineage = groups$Lineage[i],
    N_cells = groups$n_cells[i],
    Gene = genes_needed,
    Mean_Expr = as.numeric(means),
    Pct_Expr = as.numeric(pcts),
    stringsAsFactors = FALSE
  )
  setTxtProgressBar(pb, i)
}
close(pb)

expr_means <- bind_rows(expr_means_list)
rm(expr_matrix, expr_means_list); gc()

# =============================================================================
# 3. Define Cell Type Pairings
# =============================================================================

cell_type_pairings <- list(
  Fibroblast_to_Cancer = list(
    source = "Fibroblast", target = "Cancer",
    ligands = c("COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1",
                "COL5A2", "COL6A1", "COL6A2", "COL6A3", "COL11A1", "COL12A1",
                "COL14A1", "COL15A1", "FN1", "TIMP1", "THBS1", "THBS2", "LAMA4",
                "VTN", "CXCL12", "TGFB1", "TGFB2", "TGFB3", "BMP2", "BMP4",
                "HGF", "FGF2", "FGF7", "PDGFA", "PDGFB", "PDGFC", "PDGFD",
                "WNT5A", "WNT5B", "IL6", "IL11", "LIF", "OSM", "GDF15",
                "INHBA", "INHBB", "DCN", "LUM", "VCAN")
  ),
  Myeloid_to_Cancer = list(
    source = "Myeloid", target = "Cancer",
    ligands = c("SPP1", "IL1B", "IL1A", "TNF", "IL10", "IL6", "OSM", "TGFB1",
                "VEGFA", "CCL2", "CCL5", "CCL18", "CXCL8", "CXCL9", "CXCL10",
                "MMP9", "MMP2", "AREG", "EGF", "HGF", "GAS6", "APOE", "FN1",
                "THBS1", "S100A8", "S100A9", "LGALS9")
  ),
  Cancer_to_Myeloid = list(
    source = "Cancer", target = "Myeloid",
    ligands = c("MIF", "CSF1", "CSF2", "CCL2", "CCL5", "CXCL12", "IL6", "IL10",
                "TGFB1", "VEGFA", "CD274", "LGALS9", "FN1", "SPP1", "GAS6",
                "IL34", "CXCL8", "CCL18", "APOE")
  ),
  TCell_to_Cancer = list(
    source = "T Lymphocyte", target = "Cancer",
    ligands = c("IFNG", "TNF", "GZMB", "FASLG", "LTB", "IL2", "IL17A",
                "TNFSF10", "CD40LG")
  ),
  Cancer_to_TCell = list(
    source = "Cancer", target = "T Lymphocyte",
    ligands = c("CD274", "PDCD1LG2", "TGFB1", "IL10", "LGALS9", "MIF", "FAS",
                "TNFSF10", "CD80", "CD86", "IL6", "IL15", "CCL5")
  ),
  Cancer_to_Endothelial = list(
    source = "Cancer", target = c("Endothelial - Vascular", "Endothelial - Lymphatic"),
    ligands = c("VEGFA", "VEGFB", "VEGFC", "VEGFD", "PGF", "FGF2", "FGF1",
                "PDGFB", "CXCL8", "CXCL12", "HGF", "EFNB2", "DLL4", "JAG1",
                "TGFB1", "NRG1")
  ),
  Endothelial_to_Cancer = list(
    source = c("Endothelial - Vascular", "Endothelial - Lymphatic"), target = "Cancer",
    ligands = c("DLL4", "JAG1", "EFNB2", "EFNA1", "IL6", "IL11", "CXCL12",
                "TGFB1", "BMP2", "BMP4", "FGF2", "NRG1")
  ),
  Cancer_Autocrine = list(
    source = "Cancer", target = "Cancer",
    ligands = c("EGF", "AREG", "NRG1", "TGFB1", "TGFB2", "WNT3", "WNT7A",
                "WNT7B", "WNT5A", "WNT5B", "FGF1", "FGF2", "IGF1", "IGF2",
                "HGF", "MIF", "CDH1", "CXCL12", "IL6", "IL11", "LIF", "BMP2",
                "BMP4", "GDF15", "JAG1", "DLL1")
  ),
  BCell_to_Cancer = list(
    source = "B Lymphocyte", target = "Cancer",
    ligands = c("TNFSF13B", "CD40LG", "LTB", "TNF", "IL10")
  )
)

# =============================================================================
# 4. Calculate L-R Scores
# =============================================================================

cat("\nCalculating L-R scores...\n")

all_results <- list()
result_idx <- 1

for (pairing_name in names(cell_type_pairings)) {
  pairing <- cell_type_pairings[[pairing_name]]

  relevant_pairs <- lr_pairs %>% filter(Ligand %in% pairing$ligands)
  if (nrow(relevant_pairs) == 0) next

  # Weighted mean across multi-lineage sources/targets
  source_expr <- expr_means %>%
    filter(Lineage %in% pairing$source) %>%
    group_by(Patient, Timepoint, Gene) %>%
    summarise(Ligand_Mean = weighted.mean(Mean_Expr, N_cells),
              Ligand_Pct = weighted.mean(Pct_Expr, N_cells),
              N_Source = sum(N_cells), .groups = "drop")

  target_expr <- expr_means %>%
    filter(Lineage %in% pairing$target) %>%
    group_by(Patient, Timepoint, Gene) %>%
    summarise(Receptor_Mean = weighted.mean(Mean_Expr, N_cells),
              Receptor_Pct = weighted.mean(Pct_Expr, N_cells),
              N_Target = sum(N_cells), .groups = "drop")

  for (i in 1:nrow(relevant_pairs)) {
    lig_data <- source_expr %>% filter(Gene == relevant_pairs$Ligand[i]) %>%
      select(Patient, Timepoint, Ligand_Mean, Ligand_Pct, N_Source)
    rec_data <- target_expr %>% filter(Gene == relevant_pairs$Receptor[i]) %>%
      select(Patient, Timepoint, Receptor_Mean, Receptor_Pct, N_Target)

    lr_scores <- inner_join(lig_data, rec_data, by = c("Patient", "Timepoint")) %>%
      mutate(Score = Ligand_Mean * Receptor_Mean,
             Pairing = pairing_name,
             Ligand = relevant_pairs$Ligand[i],
             Receptor = relevant_pairs$Receptor[i],
             LR_Pair = paste0(relevant_pairs$Ligand[i], "\u2192", relevant_pairs$Receptor[i]),
             Source = paste(pairing$source, collapse = "/"),
             Target = paste(pairing$target, collapse = "/")) %>%
      filter(N_Source >= 10, N_Target >= 10)

    if (nrow(lr_scores) > 0) {
      all_results[[result_idx]] <- lr_scores
      result_idx <- result_idx + 1
    }
  }
}

results_df <- bind_rows(all_results)
cat(sprintf("  %d scores across %d unique L-R pairs\n",
            nrow(results_df), length(unique(results_df$LR_Pair))))

# =============================================================================
# 5. Calculate Temporal Fold Changes (First vs Last Biopsy)
# =============================================================================

cat("\nCalculating fold changes...\n")

timepoint_order <- c("Bx1", "Bx2", "Bx3", "Bx4")

fc_results <- results_df %>%
  group_by(Pairing, Ligand, Receptor, LR_Pair, Source, Target, Patient) %>%
  filter(n_distinct(Timepoint) >= 2) %>%
  arrange(factor(Timepoint, levels = timepoint_order)) %>%
  summarise(Timepoint_Early = first(Timepoint),
            Timepoint_Late = last(Timepoint),
            Score_Early = first(Score),
            Score_Late = last(Score),
            .groups = "drop") %>%
  filter(Score_Early > 0) %>%
  mutate(FC = Score_Late / Score_Early,
         Log2FC = log2(FC))

# =============================================================================
# 6. Consistency Classification Across Patients
# =============================================================================

cat("Classifying consistency across patients...\n")

consistent_analysis <- fc_results %>%
  group_by(Pairing, Ligand, Receptor, LR_Pair, Source, Target) %>%
  summarise(
    N_Patients = n(),
    N_Increased = sum(Log2FC > 0.5, na.rm = TRUE),
    N_Decreased = sum(Log2FC < -0.5, na.rm = TRUE),
    N_Stable = sum(abs(Log2FC) <= 0.5, na.rm = TRUE),
    Mean_Log2FC = mean(Log2FC, na.rm = TRUE),
    SD_Log2FC = sd(Log2FC, na.rm = TRUE),
    Median_Log2FC = median(Log2FC, na.rm = TRUE),
    Patients = paste(unique(Patient), collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(
    Direction = case_when(
      N_Increased >= 3 ~ "Consistent_INCREASE",
      N_Decreased >= 3 ~ "Consistent_DECREASE",
      N_Increased >= 2 & N_Decreased == 0 ~ "Trend_INCREASE",
      N_Decreased >= 2 & N_Increased == 0 ~ "Trend_DECREASE",
      TRUE ~ "Mixed/Variable"
    ),
    Consistency_Score = case_when(
      N_Increased >= 3 ~ N_Increased / N_Patients,
      N_Decreased >= 3 ~ N_Decreased / N_Patients,
      TRUE ~ 0
    )
  ) %>%
  arrange(desc(Consistency_Score), desc(abs(Mean_Log2FC)))

# =============================================================================
# 7. Save Results
# =============================================================================

output_dir <- "~/Desktop/Spatial CosmX project/Additional existing results/Temporal_AllPatients/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(results_df, file.path(output_dir, "LR_Scores_AllPatients_AllTimepoints.csv"),
          row.names = FALSE)
write.csv(fc_results, file.path(output_dir, "LR_FoldChanges_PerPatient.csv"),
          row.names = FALSE)
write.csv(consistent_analysis, file.path(output_dir, "LR_ConsistentChanges_Summary.csv"),
          row.names = FALSE)

cat("\n=== COMPLETE ===\n")
cat(sprintf("  Scores: %d entries\n", nrow(results_df)))
cat(sprintf("  Fold changes: %d entries\n", nrow(fc_results)))
cat(sprintf("  Consistent increase: %d pairs\n",
            sum(consistent_analysis$Direction == "Consistent_INCREASE")))
cat(sprintf("  Consistent decrease: %d pairs\n",
            sum(consistent_analysis$Direction == "Consistent_DECREASE")))
