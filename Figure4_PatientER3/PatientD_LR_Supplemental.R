# =============================================================================
# PATIENT D L-R SUPPLEMENTAL FIGURE - ALL SIGNIFICANT CHANGES
# Purpose: Comprehensive view of all L-R changes with expression/lineage filters
# =============================================================================

library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(circlize)

# -----------------------------------------------------------------------------
# 1. LOAD DATA
# -----------------------------------------------------------------------------
cat("Loading data...\n")
seurat_obj <- readRDS("~/Desktop/Spatial CosmX project/code_of_publication/data/CosMx_SMMART_345k_clean.rds")
lr_pairs <- read.csv("~/Topic_UMAP_Publication/LR_Analysis/CosMx_Panel_LR_Pairs.csv", stringsAsFactors = FALSE)

meta <- seurat_obj@meta.data %>% mutate(cell_id = rownames(seurat_obj@meta.data))
expr_matrix <- GetAssayData(seurat_obj, layer = "data")
meta_D <- meta %>% filter(Patient == "Patient_D")

available_genes <- rownames(expr_matrix)
lr_pairs_filtered <- lr_pairs %>%
  filter(Ligand %in% available_genes & Receptor %in% available_genes) %>%
  distinct(Ligand, Receptor)

# -----------------------------------------------------------------------------
# 2. PRE-COMPUTE EXPRESSION MEANS
# -----------------------------------------------------------------------------
cat("Pre-computing expression...\n")
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
    source = "Fibroblast", target = "Cancer",
    ligands = c("COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL5A2",
                "COL6A1", "COL6A2", "COL6A3", "FN1", "TIMP1", "THBS1", "THBS2",
                "LAMA4", "VTN", "CXCL12", "TGFB1", "TGFB2", "TGFB3", "BMP2", "BMP4",
                "HGF", "FGF2", "FGF7", "PDGFA", "PDGFB", "WNT5A", "IL6", "IL11",
                "LIF", "OSM", "GDF15", "DCN", "LUM", "VCAN", "INHBA")
  ),
  "Myeloid → Cancer" = list(
    source = "Myeloid", target = "Cancer",
    ligands = c("SPP1", "IL1B", "IL1A", "TNF", "IL10", "IL6", "OSM", "TGFB1",
                "VEGFA", "CCL2", "CCL5", "CCL18", "CXCL8", "CXCL9", "CXCL10",
                "MMP9", "MMP2", "AREG", "EGF", "HGF", "GAS6", "APOE", "THBS1",
                "S100A8", "S100A9", "LGALS9", "FN1", "TNFSF10")
  ),
  "Cancer → Myeloid" = list(
    source = "Cancer", target = "Myeloid",
    ligands = c("MIF", "CSF1", "CSF2", "CCL2", "CCL5", "CXCL12", "IL6", "IL10",
                "TGFB1", "VEGFA", "CD274", "LGALS9", "FN1", "SPP1", "GAS6",
                "IL34", "CXCL8", "CCL18", "CX3CL1", "THBS1")
  ),
  "Cancer → Endothelial" = list(
    source = "Cancer", target = "Endothelial - Vascular",
    ligands = c("VEGFA", "VEGFB", "VEGFC", "VEGFD", "PGF", "FGF2", "FGF1",
                "PDGFB", "CXCL8", "CXCL12", "HGF", "EFNB2", "DLL4", "JAG1",
                "TGFB1", "NRG1", "ANGPT1", "ANGPT2")
  ),
  "T Cell → Cancer" = list(
    source = "T Lymphocyte", target = "Cancer",
    ligands = c("IFNG", "TNF", "GZMB", "FASLG", "LTB", "IL2", "IL17A",
                "TNFSF10", "CD40LG", "TNFSF14", "CCL5")
  ),
  "Cancer → T Cell" = list(
    source = "Cancer", target = "T Lymphocyte",
    ligands = c("CD274", "PDCD1LG2", "TGFB1", "IL10", "LGALS9", "MIF",
                "FAS", "TNFSF10", "CD80", "CD86", "IL6", "IL15", "CCL5",
                "CXCL9", "CXCL10", "CXCL12")
  ),
  "Cancer Autocrine" = list(
    source = "Cancer", target = "Cancer",
    ligands = c("EGF", "AREG", "NRG1", "TGFB1", "TGFB2", "WNT3", "WNT7A",
                "WNT7B", "WNT5A", "WNT5B", "FGF1", "FGF2", "IGF1", "IGF2",
                "HGF", "MIF", "CDH1", "CXCL12", "IL6", "IL11", "LIF",
                "BMP2", "BMP4", "GDF15", "JAG1", "DLL1", "EFNB2", "EREG")
  ),
  "Endothelial → Cancer" = list(
    source = "Endothelial - Vascular", target = "Cancer",
    ligands = c("DLL4", "JAG1", "EFNB2", "EFNA1", "IL6", "IL11", "CXCL12",
                "TGFB1", "BMP2", "BMP4", "FGF2", "NRG1", "WNT5A", "PDGFB",
                "ANGPT1", "ANGPT2", "WIF1")
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
# 4. CALCULATE ALL L-R SCORES
# -----------------------------------------------------------------------------
cat("Calculating L-R scores...\n")
timepoints <- c("Bx1", "Bx2", "Bx3")
all_scores <- list()

for(pairing_name in names(cell_type_pairings)) {
  pairing <- cell_type_pairings[[pairing_name]]
  relevant_lr_pairs <- lr_pairs_filtered %>% filter(Ligand %in% pairing$ligands)
  if(nrow(relevant_lr_pairs) == 0) next
  
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
      mutate(Score = Ligand_Mean * Receptor_Mean, Pairing = pairing_name,
             Timepoint = tp, LR_Pair = paste0(Ligand, "→", Receptor))
    all_scores[[length(all_scores) + 1]] <- scores
  }
}

scores_df <- bind_rows(all_scores)

# -----------------------------------------------------------------------------
# 5. CALCULATE TEMPORAL CHANGES & APPLY FILTERS
# -----------------------------------------------------------------------------
cat("Calculating changes and applying filters...\n")

# FILTER PARAMETERS - adjust these as needed
MIN_EXPRESSION <- 0.05   # Minimum score in at least one timepoint
MIN_FC <- 0.5            # Minimum |Log2FC| to include (~1.4x change)

scores_wide <- scores_df %>%
  select(LR_Pair, Pairing, Ligand, Receptor, Timepoint, Score) %>%
  pivot_wider(names_from = Timepoint, values_from = Score, names_prefix = "Score_")

temporal_changes <- scores_wide %>%
  filter(!is.na(Score_Bx1), !is.na(Score_Bx3)) %>%
  mutate(
    Log2FC = ifelse(Score_Bx1 > 0, log2(Score_Bx3 / Score_Bx1), NA),
    Max_Score = pmax(Score_Bx1, coalesce(Score_Bx2, 0), Score_Bx3),
    Mean_Score = (Score_Bx1 + coalesce(Score_Bx2, 0) + Score_Bx3) / (2 + as.numeric(!is.na(Score_Bx2)))
  ) %>%
  filter(!is.na(Log2FC))

# Apply filters
filtered_changes <- temporal_changes %>%
  filter(Max_Score >= MIN_EXPRESSION) %>%
  filter(abs(Log2FC) >= MIN_FC) %>%
  mutate(
    Direction = ifelse(Log2FC > 0, "Increased", "Decreased"),
    Pairing = factor(Pairing, levels = names(cell_type_pairings))
  )

cat("Total pairs after filtering:", nrow(filtered_changes), "\n")
cat("  Increased:", sum(filtered_changes$Direction == "Increased"), "\n")
cat("  Decreased:", sum(filtered_changes$Direction == "Decreased"), "\n")

# -----------------------------------------------------------------------------
# 6. OUTPUT DIRECTORY
# -----------------------------------------------------------------------------
output_dir <- "~/Topic_UMAP_Publication/LR_Analysis/PatientD_Analysis/LR_Supplemental/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# VISUALIZATION 1: CLUSTERED HEATMAP (ComplexHeatmap)
# =============================================================================
cat("\n=== Creating clustered heatmap ===\n")

# Prepare matrix
heatmap_scores <- scores_df %>%
  filter(paste(LR_Pair, Pairing) %in% paste(filtered_changes$LR_Pair, filtered_changes$Pairing)) %>%
  select(LR_Pair, Pairing, Timepoint, Score) %>%
  unite("Row_ID", LR_Pair, Pairing, sep = " | ") %>%
  pivot_wider(names_from = Timepoint, values_from = Score)

hm_matrix <- as.matrix(heatmap_scores[, c("Bx1", "Bx2", "Bx3")])
rownames(hm_matrix) <- heatmap_scores$Row_ID
hm_matrix_scaled <- t(scale(t(hm_matrix)))
hm_matrix_scaled[is.na(hm_matrix_scaled)] <- 0

# Row annotations
row_info <- filtered_changes %>%
  mutate(Row_ID = paste(LR_Pair, Pairing, sep = " | ")) %>%
  select(Row_ID, Pairing, Direction, Log2FC)

row_annot_df <- data.frame(
  Pairing = row_info$Pairing[match(rownames(hm_matrix_scaled), row_info$Row_ID)],
  Direction = row_info$Direction[match(rownames(hm_matrix_scaled), row_info$Row_ID)],
  Log2FC = row_info$Log2FC[match(rownames(hm_matrix_scaled), row_info$Row_ID)]
)

col_fun <- colorRamp2(c(-2, 0, 2), c("#3575B5", "white", "#D73027"))
fc_col_fun <- colorRamp2(c(-3, 0, 3), c("#3575B5", "white", "#D73027"))

ha_row <- rowAnnotation(
  Category = row_annot_df$Pairing,
  Direction = row_annot_df$Direction,
  Log2FC = row_annot_df$Log2FC,
  col = list(
    Category = category_colors,
    Direction = c("Increased" = "#D73027", "Decreased" = "#4575B4"),
    Log2FC = fc_col_fun
  ),
  show_annotation_name = TRUE,
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 8)
)

pdf(file.path(output_dir, "Supplemental_Clustered_Heatmap.pdf"), width = 10, height = 14)

ht <- Heatmap(
  hm_matrix_scaled,
  name = "Z-score",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  show_row_names = FALSE,
  column_names_rot = 0,
  column_names_centered = TRUE,
  left_annotation = ha_row,
  row_split = row_annot_df$Direction,
  row_gap = unit(2, "mm"),
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_title = paste0("Patient D: All Significant L-R Changes (n=", nrow(filtered_changes), 
                        ")\nExpression ≥", MIN_EXPRESSION, " | |Log2FC| ≥", MIN_FC),
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(title = "Z-score", legend_height = unit(4, "cm"))
)

draw(ht, annotation_legend_side = "right")
dev.off()
cat("Saved: Supplemental_Clustered_Heatmap.pdf\n")

# =============================================================================
# VISUALIZATION 2: DOT PLOT
# =============================================================================
cat("\n=== Creating dot plot ===\n")

dot_data <- filtered_changes %>%
  select(LR_Pair, Pairing, Direction, Log2FC, Score_Bx1, Score_Bx2, Score_Bx3) %>%
  pivot_longer(cols = starts_with("Score_"), names_to = "Timepoint", values_to = "Score") %>%
  mutate(Timepoint = gsub("Score_", "", Timepoint),
         Timepoint = factor(Timepoint, levels = c("Bx1", "Bx2", "Bx3"))) %>%
  filter(!is.na(Score))

pair_order <- filtered_changes %>% arrange(Pairing, desc(Log2FC)) %>% pull(LR_Pair) %>% unique()
dot_data$LR_Pair <- factor(dot_data$LR_Pair, levels = rev(pair_order))

p_dot <- ggplot(dot_data, aes(x = Timepoint, y = LR_Pair)) +
  geom_point(aes(size = Score, color = Log2FC)) +
  scale_size_continuous(range = c(0.5, 4), name = "Score") +
  scale_color_gradient2(low = "#3575B5", mid = "white", high = "#D73027", 
                        midpoint = 0, limits = c(-3, 3), oob = scales::squish, name = "Log2FC") +
  facet_grid(Pairing ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(title = "Patient D: L-R Pair Expression and Change",
       subtitle = paste0("n = ", nrow(filtered_changes), " pairs"), x = "", y = "") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 5),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7, face = "bold"),
        strip.placement = "outside", panel.grid.major.y = element_blank())

n_pairs <- nrow(filtered_changes)
ggsave(file.path(output_dir, "Supplemental_DotPlot.pdf"), 
       p_dot, width = 8, height = max(12, n_pairs * 0.12))
cat("Saved: Supplemental_DotPlot.pdf\n")

# =============================================================================
# VISUALIZATION 3: VOLCANO PLOT
# =============================================================================
cat("\n=== Creating volcano plot ===\n")

p_volcano <- ggplot(filtered_changes, aes(x = Log2FC, y = log10(Mean_Score + 0.01))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = Pairing), alpha = 0.7, size = 2) +
  geom_text(data = filtered_changes %>% group_by(Pairing) %>% 
              slice_max(abs(Log2FC), n = 2) %>% ungroup(),
            aes(label = LR_Pair), size = 2, hjust = -0.1, check_overlap = TRUE) +
  scale_color_manual(values = category_colors) +
  labs(title = "Patient D: L-R Changes Overview",
       x = "Log2 Fold Change (Bx3 / Bx1)", y = "Log10 Mean Score", color = "Interaction") +
  theme_classic() + theme(legend.position = "bottom") + guides(color = guide_legend(nrow = 2))

ggsave(file.path(output_dir, "Supplemental_Volcano.pdf"), p_volcano, width = 10, height = 8)
cat("Saved: Supplemental_Volcano.pdf\n")

# =============================================================================
# VISUALIZATION 4: SUMMARY BAR PLOT
# =============================================================================
cat("\n=== Creating summary bar plot ===\n")

summary_counts <- filtered_changes %>%
  group_by(Pairing, Direction) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(N_signed = ifelse(Direction == "Decreased", -N, N))

p_summary <- ggplot(summary_counts, aes(x = Pairing, y = N_signed, fill = Direction)) +
  geom_col(width = 0.7) + geom_hline(yintercept = 0) +
  geom_text(aes(label = abs(N_signed), y = N_signed + ifelse(Direction == "Increased", 1, -1)), size = 3) +
  scale_fill_manual(values = c("Increased" = "#D73027", "Decreased" = "#4575B4")) +
  coord_flip() +
  labs(title = "Number of Significant L-R Changes per Category", x = "", y = "Number of pairs") +
  theme_classic() + theme(legend.position = "bottom")

ggsave(file.path(output_dir, "Supplemental_Summary_Barplot.pdf"), p_summary, width = 10, height = 6)
cat("Saved: Supplemental_Summary_Barplot.pdf\n")

# =============================================================================
# VISUALIZATION 5: FACETED GGPLOT HEATMAP (readable row names)
# =============================================================================
cat("\n=== Creating faceted heatmap ===\n")

facet_data <- scores_df %>%
  filter(paste(LR_Pair, Pairing) %in% paste(filtered_changes$LR_Pair, filtered_changes$Pairing)) %>%
  group_by(LR_Pair, Pairing) %>%
  mutate(Zscore = (Score - mean(Score)) / sd(Score), Zscore = ifelse(is.na(Zscore), 0, Zscore)) %>%
  ungroup() %>%
  left_join(filtered_changes %>% select(LR_Pair, Pairing, Direction, Log2FC), by = c("LR_Pair", "Pairing")) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("Bx1", "Bx2", "Bx3")),
         Pairing = factor(Pairing, levels = names(cell_type_pairings))) %>%
  arrange(Pairing, Direction, desc(Log2FC)) %>%
  mutate(LR_Pair = factor(LR_Pair, levels = unique(LR_Pair)))

p_facet <- ggplot(facet_data, aes(x = Timepoint, y = LR_Pair, fill = Zscore)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#3575B5", mid = "white", high = "#D73027",
                       midpoint = 0, limits = c(-2, 2), oob = scales::squish) +
  facet_grid(Pairing ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(title = "Patient D: All Significant L-R Changes",
       subtitle = paste0("n = ", nrow(filtered_changes), " pairs | Ordered by direction then FC"),
       x = "", y = "", fill = "Z-score") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 4), axis.text.x = element_text(size = 10),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7, face = "bold"),
        strip.placement = "outside", strip.background = element_rect(fill = "grey95", color = NA),
        panel.grid = element_blank(), panel.spacing = unit(0.3, "lines"))

ggsave(file.path(output_dir, "Supplemental_Faceted_Heatmap.pdf"), 
       p_facet, width = 7, height = max(14, n_pairs * 0.12))
cat("Saved: Supplemental_Faceted_Heatmap.pdf\n")

# =============================================================================
# SAVE DATA TABLES
# =============================================================================
cat("\n=== Saving tables ===\n")

write.csv(filtered_changes %>% arrange(Pairing, Direction, desc(abs(Log2FC))) %>%
            select(Pairing, LR_Pair, Ligand, Receptor, Direction, Score_Bx1, Score_Bx2, Score_Bx3, Log2FC),
          file.path(output_dir, "Supplemental_AllChanges.csv"), row.names = FALSE)

write.csv(data.frame(
  Parameter = c("MIN_EXPRESSION", "MIN_FC", "Total_Pairs", "Increased", "Decreased"),
  Value = c(MIN_EXPRESSION, MIN_FC, nrow(filtered_changes), 
            sum(filtered_changes$Direction == "Increased"), sum(filtered_changes$Direction == "Decreased"))
), file.path(output_dir, "Supplemental_FilterCriteria.csv"), row.names = FALSE)

cat("\n=== SUMMARY ===\n")
cat("Total pairs:", nrow(filtered_changes), "\n")
cat("Output directory:", output_dir, "\n")
cat("=== DONE ===\n")
