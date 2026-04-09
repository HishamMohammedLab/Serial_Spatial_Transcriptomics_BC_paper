# ==============================================================================
# Patient B: Bar Plots with Error Bars and Significance Testing
# 1. Percent Positive Cells
# 2. Mean Expression
# Grouped by functional category, comparing Bx1 vs Bx4
# ==============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

output_dir <- "LR_Analysis/PatientB_Analysis/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load cancer cells with topics
cat("Loading data...\n")
obj <- subset(readRDS("data/CosMx_SMMART_345k_clean.rds"), broad_lineage == "Cancer")

# Subset to Patient B
patB <- subset(obj, subset = Patient == "Patient_B")
cat("Patient B cancer cells:", ncol(patB), "cells\n")
cat("By timepoint:\n")
print(table(patB$Timepoint))

# Verify CNV+ only
cat("\nClone distribution:\n")
print(table(patB$Clone))

# ==============================================================================
# Define Key Marker Genes by Functional Category
# ==============================================================================

gene_categories <- list(
  "Basal Keratins" = c("KRT5", "KRT15", "KRT16", "KRT17", "KRT23"),
  "EMT / Invasion" = c("VIM", "CD44", "MMP7", "TM4SF1", "SERPINA3"),
  "Stress / Survival" = c("CRYAB", "CD59", "HMGB2"),
  "Luminal" = c("GATA3", "ESR1", "KRT8", "XBP1", "ERBB3", "S100P"),
  "Proliferation" = c("CCND1", "MKI67", "PCNA", "HSPB1"),
  "Immune / HLA" = c("B2M", "HLA-A", "HLA-B", "HLA-C", "CD274")
)

# Flatten and check availability
all_genes <- unique(unlist(gene_categories))
available_genes <- all_genes[all_genes %in% rownames(patB)]
missing_genes <- all_genes[!all_genes %in% rownames(patB)]

cat("\nAvailable genes:", length(available_genes), "/", length(all_genes), "\n")
if(length(missing_genes) > 0) {
  cat("Missing genes:", paste(missing_genes, collapse = ", "), "\n")
}

# Update categories
gene_categories <- lapply(gene_categories, function(x) x[x %in% available_genes])
gene_categories <- gene_categories[sapply(gene_categories, length) > 0]

# ==============================================================================
# Extract Expression Data
# ==============================================================================

expr_matrix <- GetAssayData(patB, layer = "data")

# Create long-format data
expr_df <- as.data.frame(t(as.matrix(expr_matrix[available_genes, , drop = FALSE])))
expr_df$Timepoint <- patB$Timepoint
expr_df$cell_id <- rownames(expr_df)

# Get cell counts per timepoint
n_bx1 <- sum(expr_df$Timepoint == "Bx1")
n_bx4 <- sum(expr_df$Timepoint == "Bx4")
cat("\nCells: Bx1 =", n_bx1, ", Bx4 =", n_bx4, "\n")

# ==============================================================================
# Calculate Statistics with Error Bars
# ==============================================================================

cat("\nCalculating statistics with error bars...\n")

stats_list <- list()
pval_list <- list()

for(gene in available_genes) {
  
  bx1_vals <- expr_df[expr_df$Timepoint == "Bx1", gene]
  bx4_vals <- expr_df[expr_df$Timepoint == "Bx4", gene]
  
  # ----- Mean Expression Stats -----
  mean_bx1 <- mean(bx1_vals, na.rm = TRUE)
  mean_bx4 <- mean(bx4_vals, na.rm = TRUE)
  
  # Standard error of mean
  se_bx1 <- sd(bx1_vals, na.rm = TRUE) / sqrt(length(bx1_vals))
  se_bx4 <- sd(bx4_vals, na.rm = TRUE) / sqrt(length(bx4_vals))
  
  # Wilcoxon test for expression
  wilcox_expr <- wilcox.test(bx1_vals, bx4_vals)
  
  # ----- Percent Positive Stats -----
  pos_bx1 <- sum(bx1_vals > 0, na.rm = TRUE)
  pos_bx4 <- sum(bx4_vals > 0, na.rm = TRUE)
  
  pct_bx1 <- pos_bx1 / length(bx1_vals) * 100
  pct_bx4 <- pos_bx4 / length(bx4_vals) * 100
  
  # Standard error for proportion (binomial)
  se_pct_bx1 <- sqrt((pct_bx1 * (100 - pct_bx1)) / length(bx1_vals))
  se_pct_bx4 <- sqrt((pct_bx4 * (100 - pct_bx4)) / length(bx4_vals))
  
  # Chi-squared test for proportions
  contingency <- matrix(c(pos_bx1, length(bx1_vals) - pos_bx1,
                          pos_bx4, length(bx4_vals) - pos_bx4), 
                        nrow = 2, byrow = TRUE)
  chisq_test <- tryCatch(chisq.test(contingency), error = function(e) list(p.value = NA))
  
  # Store Bx1 stats
  stats_list[[paste0(gene, "_Bx1")]] <- data.frame(
    Gene = gene,
    Timepoint = "Bx1",
    Mean_Expression = mean_bx1,
    SE_Expression = se_bx1,
    Pct_Positive = pct_bx1,
    SE_Pct = se_pct_bx1,
    stringsAsFactors = FALSE
  )
  
  # Store Bx4 stats
  stats_list[[paste0(gene, "_Bx4")]] <- data.frame(
    Gene = gene,
    Timepoint = "Bx4",
    Mean_Expression = mean_bx4,
    SE_Expression = se_bx4,
    Pct_Positive = pct_bx4,
    SE_Pct = se_pct_bx4,
    stringsAsFactors = FALSE
  )
  
  # Store p-values
  pval_list[[gene]] <- data.frame(
    Gene = gene,
    p_expr = wilcox_expr$p.value,
    p_pct = chisq_test$p.value,
    stringsAsFactors = FALSE
  )
}

stats_df <- bind_rows(stats_list)
pval_df <- bind_rows(pval_list)

# Adjust p-values for multiple testing
pval_df$p_expr_adj <- p.adjust(pval_df$p_expr, method = "BH")
pval_df$p_pct_adj <- p.adjust(pval_df$p_pct, method = "BH")

# Add significance stars
pval_df <- pval_df %>%
  mutate(
    sig_expr = case_when(
      p_expr_adj < 0.001 ~ "***",
      p_expr_adj < 0.01 ~ "**",
      p_expr_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    sig_pct = case_when(
      p_pct_adj < 0.001 ~ "***",
      p_pct_adj < 0.01 ~ "**",
      p_pct_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# Add category
stats_df$Category <- NA
for(cat_name in names(gene_categories)) {
  stats_df$Category[stats_df$Gene %in% gene_categories[[cat_name]]] <- cat_name
}

# Merge p-values into stats
stats_df <- stats_df %>%
  left_join(pval_df, by = "Gene")

# ==============================================================================
# Prepare Data for Plotting
# ==============================================================================

# Order categories
category_order <- c("Basal Keratins", "EMT / Invasion", "Stress / Survival", 
                    "Luminal", "Proliferation", "Immune / HLA")
category_order <- category_order[category_order %in% unique(stats_df$Category)]

stats_df$Category <- factor(stats_df$Category, levels = category_order)

# Create gene order: by category, then by Bx4 expression within category
gene_order <- stats_df %>%
  filter(Timepoint == "Bx4") %>%
  arrange(Category, desc(Mean_Expression)) %>%
  pull(Gene) %>%
  unique()

stats_df$Gene <- factor(stats_df$Gene, levels = rev(gene_order))

# Colors (Patient B green shades)
timepoint_colors <- c("Bx1" = "#95D5B2", "Bx4" = "#1B4332")

# ==============================================================================
# Helper Function: Get significance annotation position
# ==============================================================================

get_sig_data <- function(data, value_col, se_col, sig_col) {
  data %>%
    group_by(Gene, Category) %>%
    summarise(
      y_pos = max(.data[[value_col]] + .data[[se_col]], na.rm = TRUE),
      sig = first(.data[[sig_col]]),
      .groups = "drop"
    ) %>%
    mutate(y_pos = y_pos * 1.05)  # Add 5% padding
}

# ==============================================================================
# PLOT 1: Percent Positive Cells with Error Bars
# ==============================================================================

cat("\nGenerating Percent Positive barplot with error bars...\n")

# Get significance positions for pct plot
sig_data_pct <- get_sig_data(stats_df, "Pct_Positive", "SE_Pct", "sig_pct")

p_pct <- ggplot(stats_df, aes(x = Gene, y = Pct_Positive, fill = Timepoint)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbar(
    aes(ymin = Pct_Positive - SE_Pct, ymax = Pct_Positive + SE_Pct),
    position = position_dodge(width = 0.8),
    width = 0.25,
    linewidth = 0.4
  ) +
  geom_text(
    data = sig_data_pct,
    aes(x = Gene, y = y_pos, label = sig),
    inherit.aes = FALSE,
    size = 3,
    fontface = "bold",
    hjust = 0.5
  ) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = timepoint_colors, name = "Biopsy") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  coord_flip() +
  labs(
    title = "Patient B: Percent Positive Cells by Gene",
    subtitle = "Mean ± SE | Chi-squared test: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant",
    x = NULL,
    y = "% Positive Cells"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(color = "grey40", size = 9, hjust = 0),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 9, hjust = 1),
    strip.placement = "outside",
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    panel.spacing = unit(0.4, "lines"),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.justification = "left"
  )

ggsave(file.path(output_dir, "PatientB_PctPositive_Barplot_withStats.pdf"),
       p_pct, width = 9, height = 12)

cat("Saved: PatientB_PctPositive_Barplot_withStats.pdf\n")

# ==============================================================================
# PLOT 2: Mean Expression with Error Bars
# ==============================================================================

cat("Generating Mean Expression barplot with error bars...\n")

# Get significance positions for expression plot
sig_data_expr <- get_sig_data(stats_df, "Mean_Expression", "SE_Expression", "sig_expr")

p_mean <- ggplot(stats_df, aes(x = Gene, y = Mean_Expression, fill = Timepoint)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbar(
    aes(ymin = Mean_Expression - SE_Expression, ymax = Mean_Expression + SE_Expression),
    position = position_dodge(width = 0.8),
    width = 0.25,
    linewidth = 0.4
  ) +
  geom_text(
    data = sig_data_expr,
    aes(x = Gene, y = y_pos, label = sig),
    inherit.aes = FALSE,
    size = 3,
    fontface = "bold",
    hjust = 0.5
  ) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = timepoint_colors, name = "Biopsy") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  coord_flip() +
  labs(
    title = "Patient B: Mean Expression by Gene",
    subtitle = "Mean ± SE | Wilcoxon test: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant",
    x = NULL,
    y = "Mean Normalized Expression"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(color = "grey40", size = 9, hjust = 0),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 9, hjust = 1),
    strip.placement = "outside",
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    panel.spacing = unit(0.4, "lines"),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.justification = "left"
  )

ggsave(file.path(output_dir, "PatientB_MeanExpression_Barplot_withStats.pdf"),
       p_mean, width = 9, height = 12)

cat("Saved: PatientB_MeanExpression_Barplot_withStats.pdf\n")

# ==============================================================================
# Save Full Statistics Table
# ==============================================================================

# Create comprehensive summary table
summary_table <- stats_df %>%
  select(Gene, Category, Timepoint, 
         Pct_Positive, SE_Pct,
         Mean_Expression, SE_Expression) %>%
  pivot_wider(
    names_from = Timepoint,
    values_from = c(Pct_Positive, SE_Pct, Mean_Expression, SE_Expression)
  ) %>%
  left_join(pval_df, by = "Gene") %>%
  mutate(
    Pct_Diff = Pct_Positive_Bx4 - Pct_Positive_Bx1,
    Mean_Log2FC = log2((Mean_Expression_Bx4 + 0.01) / (Mean_Expression_Bx1 + 0.01))
  ) %>%
  select(Gene, Category, 
         Pct_Positive_Bx1, Pct_Positive_Bx4, Pct_Diff, p_pct_adj, sig_pct,
         Mean_Expression_Bx1, Mean_Expression_Bx4, Mean_Log2FC, p_expr_adj, sig_expr) %>%
  arrange(Category, desc(Mean_Log2FC)) %>%
  mutate(
    across(starts_with("Pct"), ~round(., 1)),
    across(starts_with("Mean"), ~round(., 3)),
    Mean_Log2FC = round(Mean_Log2FC, 2),
    p_pct_adj = signif(p_pct_adj, 3),
    p_expr_adj = signif(p_expr_adj, 3)
  )

write.csv(summary_table, file.path(output_dir, "PatientB_Gene_Stats_Full.csv"),
          row.names = FALSE)

cat("\nSaved: PatientB_Gene_Stats_Full.csv\n")

# ==============================================================================
# Print Summary
# ==============================================================================

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("SUMMARY BY CATEGORY\n")
cat("=", rep("=", 70), "\n", sep = "")

for(cat_name in category_order) {
  cat_data <- summary_table %>% filter(Category == cat_name)
  
  cat("\n", cat_name, ":\n", sep = "")
  cat(rep("-", 50), "\n", sep = "")
  
  for(i in 1:nrow(cat_data)) {
    direction <- ifelse(cat_data$Mean_Log2FC[i] > 0, "↑", "↓")
    cat(sprintf("  %-10s %s %+.2f log2FC (%s) | %%.pos: %.1f → %.1f (%s)\n",
                cat_data$Gene[i],
                direction,
                cat_data$Mean_Log2FC[i],
                cat_data$sig_expr[i],
                cat_data$Pct_Positive_Bx1[i],
                cat_data$Pct_Positive_Bx4[i],
                cat_data$sig_pct[i]))
  }
}

cat("\n\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("OUTPUT FILES\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("\n")
cat("  1. PatientB_PctPositive_Barplot_withStats.pdf  - % positive with SE & significance\n")
cat("  2. PatientB_MeanExpression_Barplot_withStats.pdf - Mean expr with SE & significance\n")
cat("  3. PatientB_Gene_Stats_Full.csv - Complete statistics table\n")
cat("\n")
cat("Significance: BH-adjusted p-values\n")
cat("  *** p < 0.001\n")
cat("  **  p < 0.01\n")
cat("  *   p < 0.05\n")
cat("  ns  not significant\n")
cat("\n")
