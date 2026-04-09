# =============================================================================
# PATIENT D L-R HEATMAP - INTEGRATED INCREASING & DECREASING
# Purpose: Show L-R pairs changing with resistance, with selection criteria
# =============================================================================

library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. LOAD DATA
# -----------------------------------------------------------------------------
cat("Loading Seurat object...\n")
seurat_obj <- readRDS("data/CosMx_SMMART_345k_clean.rds")

lr_pairs <- read.csv("supplementary_input_data/LR_Scores/CosMx_Panel_LR_Pairs.csv", stringsAsFactors = FALSE)

meta <- seurat_obj@meta.data %>%
  mutate(cell_id = rownames(seurat_obj@meta.data))

expr_matrix <- GetAssayData(seurat_obj, layer = "data")

meta_D <- meta %>% filter(Patient == "Patient_D")

available_genes <- rownames(expr_matrix)
lr_pairs_filtered <- lr_pairs %>%
  filter(Ligand %in% available_genes & Receptor %in% available_genes) %>%
  distinct(Ligand, Receptor)

cat("L-R pairs available:", nrow(lr_pairs_filtered), "\n")

# -----------------------------------------------------------------------------
# 2. PRE-COMPUTE EXPRESSION MEANS
# -----------------------------------------------------------------------------
cat("\nPre-computing expression means...\n")

genes_needed <- unique(c(lr_pairs_filtered$Ligand, lr_pairs_filtered$Receptor))

expr_means <- meta_D %>%
  group_by(Timepoint, Lineage) %>%
  filter(n() >= 10) %>%
  group_modify(~{
    cells <- .x$cell_id
    means <- rowMeans(expr_matrix[genes_needed, cells, drop = FALSE])
    data.frame(Gene = genes_needed, Mean_Expr = as.numeric(means), N_cells = length(cells))
  }) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 3. DEFINE CELL TYPE PAIRINGS WITH RELEVANT LIGANDS
# -----------------------------------------------------------------------------
cell_type_pairings <- list(
  "Fibroblast → Cancer" = list(
    source = "Fibroblast", 
    target = "Cancer",
    ligands = c("COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL5A2",
                "COL6A1", "COL6A2", "COL6A3", "FN1", "TIMP1", "THBS1", "THBS2",
                "LAMA4", "VTN", "CXCL12", "TGFB1", "TGFB2", "TGFB3", "BMP2", "BMP4",
                "HGF", "FGF2", "FGF7", "PDGFA", "PDGFB", "WNT5A", "IL6", "IL11",
                "LIF", "OSM", "GDF15", "DCN", "LUM", "VCAN", "INHBA")
  ),
  
  "Myeloid → Cancer" = list(
    source = "Myeloid", 
    target = "Cancer",
    ligands = c("SPP1", "IL1B", "IL1A", "TNF", "IL10", "IL6", "OSM", "TGFB1",
                "VEGFA", "CCL2", "CCL5", "CCL18", "CXCL8", "CXCL9", "CXCL10",
                "MMP9", "MMP2", "AREG", "EGF", "HGF", "GAS6", "APOE", "THBS1",
                "S100A8", "S100A9", "LGALS9", "FN1", "TNFSF10")
  ),
  
  "Cancer → Myeloid" = list(
    source = "Cancer", 
    target = "Myeloid",
    ligands = c("MIF", "CSF1", "CSF2", "CCL2", "CCL5", "CXCL12", "IL6", "IL10",
                "TGFB1", "VEGFA", "CD274", "LGALS9", "FN1", "SPP1", "GAS6",
                "IL34", "CXCL8", "CCL18", "CX3CL1", "HMGB1", "THBS1")
  ),
  
  "Cancer → Endothelial" = list(
    source = "Cancer", 
    target = "Endothelial - Vascular",
    ligands = c("VEGFA", "VEGFB", "VEGFC", "VEGFD", "PGF", "FGF2", "FGF1",
                "PDGFB", "CXCL8", "CXCL12", "HGF", "EFNB2", "DLL4", "JAG1",
                "TGFB1", "NRG1", "ANGPT1", "ANGPT2")
  ),
  
  "T Cell → Cancer" = list(
    source = "T Lymphocyte", 
    target = "Cancer",
    ligands = c("IFNG", "TNF", "GZMB", "FASLG", "LTB", "IL2", "IL17A",
                "TNFSF10", "CD40LG", "TNFSF14", "IL21", "IL22", "CCL5")
  ),
  
  "Cancer → T Cell" = list(
    source = "Cancer", 
    target = "T Lymphocyte",
    ligands = c("CD274", "PDCD1LG2", "TGFB1", "IL10", "LGALS9", "MIF",
                "FAS", "TNFSF10", "CD80", "CD86", "IL6", "IL15", "CCL5",
                "CXCL9", "CXCL10", "CXCL12")
  ),
  
  "Cancer Autocrine" = list(
    source = "Cancer", 
    target = "Cancer",
    ligands = c("EGF", "AREG", "NRG1", "TGFB1", "TGFB2", "WNT3", "WNT7A",
                "WNT7B", "WNT5A", "WNT5B", "FGF1", "FGF2", "IGF1", "IGF2",
                "HGF", "MIF", "CDH1", "CXCL12", "IL6", "IL11", "LIF",
                "BMP2", "BMP4", "GDF15", "JAG1", "DLL1", "EFNB2", "EREG")
  ),
  
  "Endothelial → Cancer" = list(
    source = "Endothelial - Vascular", 
    target = "Cancer",
    ligands = c("DLL4", "JAG1", "EFNB2", "EFNA1", "IL6", "IL11", "CXCL12",
                "TGFB1", "BMP2", "BMP4", "FGF2", "NRG1", "WNT5A", "PDGFB",
                "ANGPT1", "ANGPT2", "WIF1", "SLIT2")
  )
)

category_colors <- c(
  "Fibroblast → Cancer" = "#D96941",
  "Myeloid → Cancer" = "#7570B3",
  "Cancer → Myeloid" = "#E7298A",
  "Cancer → Endothelial" = "#66A61E",
  "T Cell → Cancer" = "#E6AB02",
  "Cancer → T Cell" = "#A6761D",
  "Cancer Autocrine" = "#1B9E77",
  "Endothelial → Cancer" = "#666666"
)

# -----------------------------------------------------------------------------
# 4. CALCULATE L-R SCORES
# -----------------------------------------------------------------------------
cat("Calculating L-R scores...\n")

timepoints <- c("Bx1", "Bx2", "Bx3")
all_scores <- list()

for(pairing_name in names(cell_type_pairings)) {
  pairing <- cell_type_pairings[[pairing_name]]
  
  relevant_lr_pairs <- lr_pairs_filtered %>%
    filter(Ligand %in% pairing$ligands)
  
  if(nrow(relevant_lr_pairs) == 0) next
  
  cat("  ", pairing_name, ":", nrow(relevant_lr_pairs), "pairs\n")
  
  for(tp in timepoints) {
    source_expr <- expr_means %>%
      filter(Timepoint == tp, Lineage == pairing$source) %>%
      select(Gene, Ligand_Mean = Mean_Expr)
    
    target_expr <- expr_means %>%
      filter(Timepoint == tp, Lineage == pairing$target) %>%
      select(Gene, Receptor_Mean = Mean_Expr)
    
    if(nrow(source_expr) == 0 || nrow(target_expr) == 0) next
    
    scores <- relevant_lr_pairs %>%
      left_join(source_expr, by = c("Ligand" = "Gene")) %>%
      left_join(target_expr, by = c("Receptor" = "Gene")) %>%
      filter(!is.na(Ligand_Mean), !is.na(Receptor_Mean)) %>%
      mutate(
        Score = Ligand_Mean * Receptor_Mean,
        Pairing = pairing_name,
        Timepoint = tp,
        LR_Pair = paste0(Ligand, "→", Receptor)
      )
    
    all_scores[[length(all_scores) + 1]] <- scores
  }
}

scores_df <- bind_rows(all_scores)

# -----------------------------------------------------------------------------
# 5. CALCULATE TEMPORAL CHANGES WITH METRICS
# -----------------------------------------------------------------------------
cat("\nCalculating temporal changes...\n")

scores_wide <- scores_df %>%
  select(LR_Pair, Pairing, Ligand, Receptor, Timepoint, Score) %>%
  pivot_wider(names_from = Timepoint, values_from = Score, names_prefix = "Score_")

temporal_changes <- scores_wide %>%
  filter(!is.na(Score_Bx1), !is.na(Score_Bx3)) %>%
  mutate(
    # Fold changes
    FC_Bx1_Bx3 = ifelse(Score_Bx1 > 0, Score_Bx3 / Score_Bx1, NA),
    Log2FC_Bx1_Bx3 = ifelse(Score_Bx1 > 0, log2(Score_Bx3 / Score_Bx1), NA),
    FC_Bx2_Bx3 = ifelse(!is.na(Score_Bx2) & Score_Bx2 > 0, Score_Bx3 / Score_Bx2, NA),
    Log2FC_Bx2_Bx3 = ifelse(!is.na(Score_Bx2) & Score_Bx2 > 0, log2(Score_Bx3 / Score_Bx2), NA),
    
    # Expression metrics
    Max_Score = pmax(Score_Bx1, coalesce(Score_Bx2, 0), Score_Bx3),
    Mean_Score = (Score_Bx1 + coalesce(Score_Bx2, 0) + Score_Bx3) / (2 + as.numeric(!is.na(Score_Bx2))),
    
    # Absolute change
    Abs_Change = Score_Bx3 - Score_Bx1
  )

# -----------------------------------------------------------------------------
# 6. SELECTION CRITERIA
# -----------------------------------------------------------------------------
cat("\nApplying selection criteria...\n")

# PARAMETERS
N_PER_CATEGORY <- 4
MIN_EXPRESSION <- 0.1  # Minimum score in at least one timepoint
MIN_FC <- 0.5          # Minimum |Log2FC| to be considered changing

# Apply expression filter
expression_filtered <- temporal_changes %>%
  filter(Max_Score >= MIN_EXPRESSION)

cat("After expression filter:", nrow(expression_filtered), "of", nrow(temporal_changes), "\n")

# -----------------------------------------------------------------------------
# 7. SELECT TOP PAIRS WITH JUSTIFICATION
# -----------------------------------------------------------------------------
selected_pairs <- list()

for(pairing_name in names(cell_type_pairings)) {
  pairing_data <- expression_filtered %>% 
    filter(Pairing == pairing_name)
  
  if(nrow(pairing_data) == 0) next
  
  # TOP INCREASING: prioritize late-stage increases (Bx2→Bx3)
  top_increasing <- pairing_data %>%
    filter(Log2FC_Bx1_Bx3 > MIN_FC) %>%
    arrange(desc(coalesce(Log2FC_Bx2_Bx3, Log2FC_Bx1_Bx3))) %>%
    head(N_PER_CATEGORY) %>%
    mutate(
      Direction = "Increased",
      Selection_Reason = case_when(
        !is.na(Log2FC_Bx2_Bx3) & Log2FC_Bx2_Bx3 > 1 ~ 
          paste0("Late surge: ", round(2^Log2FC_Bx2_Bx3, 1), "x Bx2→Bx3"),
        Log2FC_Bx1_Bx3 > 1 ~ 
          paste0("Strong increase: ", round(2^Log2FC_Bx1_Bx3, 1), "x overall"),
        TRUE ~ paste0("Moderate increase: ", round(2^Log2FC_Bx1_Bx3, 1), "x")
      )
    )
  
  # TOP DECREASING: prioritize pairs that were high early but lost
  top_decreasing <- pairing_data %>%
    filter(Log2FC_Bx1_Bx3 < -MIN_FC) %>%
    filter(Score_Bx1 >= MIN_EXPRESSION | coalesce(Score_Bx2, 0) >= MIN_EXPRESSION) %>%  # Was expressed early
    arrange(Log2FC_Bx1_Bx3) %>%
    head(N_PER_CATEGORY) %>%
    mutate(
      Direction = "Decreased",
      Selection_Reason = case_when(
        !is.na(Log2FC_Bx2_Bx3) & Log2FC_Bx2_Bx3 < -1 ~ 
          paste0("Late loss: ", round(2^abs(Log2FC_Bx2_Bx3), 1), "x drop Bx2→Bx3"),
        Log2FC_Bx1_Bx3 < -1 ~ 
          paste0("Strong decrease: ", round(2^abs(Log2FC_Bx1_Bx3), 1), "x overall"),
        TRUE ~ paste0("Moderate decrease: ", round(2^abs(Log2FC_Bx1_Bx3), 1), "x")
      )
    )
  
  selected_pairs[[pairing_name]] <- bind_rows(top_increasing, top_decreasing)
}

top_pairs <- bind_rows(selected_pairs)
cat("Total pairs selected:", nrow(top_pairs), "\n")
cat("  Increasing:", sum(top_pairs$Direction == "Increased"), "\n")
cat("  Decreasing:", sum(top_pairs$Direction == "Decreased"), "\n")

# -----------------------------------------------------------------------------
# 8. PREPARE HEATMAP DATA
# -----------------------------------------------------------------------------
heatmap_data <- scores_df %>%
  inner_join(top_pairs %>% select(LR_Pair, Pairing, Direction, Selection_Reason, Log2FC_Bx1_Bx3), 
             by = c("LR_Pair", "Pairing")) %>%
  select(LR_Pair, Pairing, Direction, Selection_Reason, Log2FC_Bx1_Bx3, Timepoint, Score) %>%
  group_by(LR_Pair, Pairing) %>%
  mutate(
    Zscore = (Score - mean(Score)) / sd(Score),
    Zscore = ifelse(is.na(Zscore), 0, Zscore)
  ) %>%
  ungroup()

# Order: within each category, increasing first (sorted by FC), then decreasing (sorted by FC)
heatmap_data <- heatmap_data %>%
  arrange(Pairing, Direction, desc(Log2FC_Bx1_Bx3)) %>%
  mutate(
    LR_Pair = factor(LR_Pair, levels = unique(LR_Pair)),
    Pairing = factor(Pairing, levels = names(cell_type_pairings)),
    Timepoint = factor(Timepoint, levels = c("Bx1", "Bx2", "Bx3"))
  )

# -----------------------------------------------------------------------------
# 9. CREATE INTEGRATED HEATMAP
# -----------------------------------------------------------------------------
cat("\nCreating heatmap...\n")

output_dir <- "~/Topic_UMAP_Publication/LR_Analysis/PatientD_Analysis/LR_Temporal_Bx1_Bx3/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Direction annotation bar
direction_data <- heatmap_data %>%
  distinct(LR_Pair, Pairing, Direction) %>%
  mutate(x = 0.5)

p_direction <- ggplot(direction_data, aes(x = x, y = LR_Pair, fill = Direction)) +
  geom_tile(width = 1) +
  scale_fill_manual(values = c("Increased" = "#D73027", "Decreased" = "#4575B4"),
                    name = "Direction") +
  facet_grid(Pairing ~ ., scales = "free_y", space = "free_y") +
  theme_void() +
  theme(
    legend.position = "none",
    strip.text = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

# Main heatmap
p_heatmap <- ggplot(heatmap_data, aes(x = Timepoint, y = LR_Pair, fill = Zscore)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#3575B5", 
    mid = "white", 
    high = "#D73027",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish,
    name = "Z-score"
  ) +
  facet_grid(Pairing ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(
    title = "Patient D: L-R Signaling Changes During Treatment Resistance",
    subtitle = paste0("Top ", N_PER_CATEGORY, " increasing and decreasing pairs per lineage interaction",
                      " | Min expression: ", MIN_EXPRESSION, " | Min |Log2FC|: ", MIN_FC),
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 10),
    panel.grid = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 8),
    strip.background = element_rect(fill = "grey95", color = NA),
    legend.position = "right"
  )

# Combine
n_pairs <- length(unique(heatmap_data$LR_Pair))
plot_height <- max(10, n_pairs * 0.18 + 3)

p_combined <- p_direction + p_heatmap + 
  plot_layout(widths = c(0.05, 1), guides = "collect")

ggsave(file.path(output_dir, "LR_Heatmap_Integrated.pdf"), 
       p_combined, width = 8, height = plot_height)
cat("Saved: LR_Heatmap_Integrated.pdf\n")

# -----------------------------------------------------------------------------
# 10. SAVE DETAILED SUMMARY WITH JUSTIFICATION
# -----------------------------------------------------------------------------
summary_table <- top_pairs %>%
  select(Pairing, LR_Pair, Ligand, Receptor, Direction, Selection_Reason,
         Score_Bx1, Score_Bx2, Score_Bx3, Log2FC_Bx1_Bx3, Log2FC_Bx2_Bx3) %>%
  mutate(
    FC_Bx1_Bx3 = round(2^Log2FC_Bx1_Bx3, 2),
    FC_Bx2_Bx3 = round(2^coalesce(Log2FC_Bx2_Bx3, 0), 2)
  ) %>%
  arrange(Pairing, Direction, desc(abs(Log2FC_Bx1_Bx3)))

write.csv(summary_table, file.path(output_dir, "LR_Heatmap_Integrated_Summary.csv"), row.names = FALSE)
cat("Saved: LR_Heatmap_Integrated_Summary.csv\n")

# -----------------------------------------------------------------------------
# 11. CREATE SELECTION CRITERIA TABLE
# -----------------------------------------------------------------------------
criteria_summary <- data.frame(
  Parameter = c("Pairs per category", "Min expression (any Bx)", "Min |Log2FC|", 
                "Increasing selection", "Decreasing selection"),
  Value = c(N_PER_CATEGORY, MIN_EXPRESSION, MIN_FC,
            "Prioritize late-stage surge (Bx2→Bx3)",
            "Must have early expression (Bx1 or Bx2)")
)

write.csv(criteria_summary, file.path(output_dir, "LR_Selection_Criteria.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 12. PRINT SUMMARY
# -----------------------------------------------------------------------------
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("INTEGRATED L-R HEATMAP SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

cat("\n--- SELECTION CRITERIA ---\n")
cat("• Minimum expression (any timepoint):", MIN_EXPRESSION, "\n")
cat("• Minimum |Log2FC| for selection:", MIN_FC, "\n")
cat("• Pairs per category:", N_PER_CATEGORY, "increasing +", N_PER_CATEGORY, "decreasing\n")
cat("• Increasing: prioritized by Bx2→Bx3 change (late-stage surge)\n")
cat("• Decreasing: required early expression, prioritized by magnitude of loss\n")

cat("\n--- PAIRS PER CATEGORY ---\n")
summary_by_cat <- top_pairs %>%
  group_by(Pairing, Direction) %>%
  summarise(N = n(), .groups = "drop") %>%
  pivot_wider(names_from = Direction, values_from = N, values_fill = 0)
print(summary_by_cat)

cat("\n--- TOP INCREASING (by fold change) ---\n")
top_pairs %>%
  filter(Direction == "Increased") %>%
  arrange(desc(Log2FC_Bx1_Bx3)) %>%
  head(10) %>%
  select(Pairing, LR_Pair, Selection_Reason) %>%
  print()

cat("\n--- TOP DECREASING (by fold change) ---\n")
top_pairs %>%
  filter(Direction == "Decreased") %>%
  arrange(Log2FC_Bx1_Bx3) %>%
  head(10) %>%
  select(Pairing, LR_Pair, Selection_Reason) %>%
  print()

cat("\n=== DONE ===\n")
cat("Output directory:", output_dir, "\n")
