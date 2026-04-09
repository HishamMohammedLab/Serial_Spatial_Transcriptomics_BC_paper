# ==============================================================================
# FINAL BARPLOT: CURATED FIBROBLAST DE GENES
# Cancer Adjacent vs Cancer Distal - Biologically plausible markers only
# ==============================================================================

library(Seurat)
library(tidyverse)
library(patchwork)

set.seed(42)

# ==============================================================================
# LOAD DATA
# ==============================================================================

message("Loading data...")

sobj <- readRDS("data/CosMx_SMMART_345k_clean.rds")

# Subset to Patient C Bx2 fibroblasts
patC_fib <- subset(sobj, Patient == "Patient_C" & Timepoint == "Bx2" & Lineage == "Fibroblast")
message("Patient C Bx2 Fibroblasts: ", ncol(patC_fib))

# Create TWO-GROUP zone classification
patC_fib$zone <- case_when(
  patC_fib$spatial_domain %in% c("2", "3") ~ "Cancer Distal",
  TRUE ~ "Cancer Adjacent"
)

message("\nCells per zone:")
print(table(patC_fib$zone))

rm(sobj)
gc()

# Output directory
out_dir <- "Final_Figures/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# DEFINE CURATED GENES (Biologically plausible for fibroblasts)
# ==============================================================================

# Cancer Adjacent - Active CAF signature
adjacent_genes <- c("COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", 
                    "FN1", "THBS1", "HSPB1")

# Cancer Distal - Quiescent CAF signature  
distal_genes <- c("APOD", "NDRG1", "NPPC", "MT1X", "CXCL2", "IL1RN", "IL18")

# Check availability
all_genes_available <- rownames(patC_fib)

adjacent_avail <- adjacent_genes[adjacent_genes %in% all_genes_available]
distal_avail <- distal_genes[distal_genes %in% all_genes_available]

message("\nAdjacent genes available (", length(adjacent_avail), "/", length(adjacent_genes), "): ")
message("  ", paste(adjacent_avail, collapse = ", "))

message("\nDistal genes available (", length(distal_avail), "/", length(distal_genes), "): ")
message("  ", paste(distal_avail, collapse = ", "))

all_markers <- c(adjacent_avail, distal_avail)

# ==============================================================================
# EXTRACT AND NORMALIZE EXPRESSION
# ==============================================================================

expr_data <- FetchData(patC_fib, vars = all_markers)
expr_data$zone <- patC_fib$zone
expr_data$nCount_RNA <- patC_fib$nCount_RNA

# Normalize by nCount
expr_norm <- expr_data %>%
  mutate(across(all_of(all_markers), ~ . / nCount_RNA * 1000))

# ==============================================================================
# CALCULATE STATS BY ZONE
# ==============================================================================

stats_by_zone <- expr_norm %>%
  pivot_longer(cols = all_of(all_markers), names_to = "gene", values_to = "expr") %>%
  group_by(zone, gene) %>%
  summarise(
    n_cells = n(),
    mean_expr = mean(expr, na.rm = TRUE),
    se_expr = sd(expr, na.rm = TRUE) / sqrt(n()),
    pct_expressing = mean(expr > 0) * 100,
    .groups = "drop"
  )

# Add category
stats_by_zone <- stats_by_zone %>%
  mutate(
    category = case_when(
      gene %in% adjacent_avail ~ "Cancer Adjacent\n(Active CAF)",
      gene %in% distal_avail ~ "Cancer Distal\n(Quiescent CAF)",
      TRUE ~ "Other"
    )
  )

stats_by_zone$zone <- factor(stats_by_zone$zone, 
                              levels = c("Cancer Adjacent", "Cancer Distal"))

# ==============================================================================
# COLOR PALETTE
# ==============================================================================

zone_colors <- c(
  "Cancer Adjacent" = "#9E503C",
  "Cancer Distal" = "#2D5A5A"
)

# ==============================================================================
# CREATE SEPARATE PANELS WITH DIFFERENT SCALES
# ==============================================================================

message("\nCreating final figure...")

# Order genes by expression difference
adjacent_order <- stats_by_zone %>%
  filter(gene %in% adjacent_avail, zone == "Cancer Adjacent") %>%
  arrange(desc(mean_expr)) %>%
  pull(gene)

distal_order <- stats_by_zone %>%
  filter(gene %in% distal_avail, zone == "Cancer Distal") %>%
  arrange(desc(mean_expr)) %>%
  pull(gene)

# Panel 1: Cancer Adjacent genes
p_adjacent <- stats_by_zone %>%
  filter(gene %in% adjacent_avail) %>%
  mutate(gene = factor(gene, levels = adjacent_order)) %>%
  ggplot(aes(x = gene, y = mean_expr, fill = zone)) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = pmax(0, mean_expr - se_expr), ymax = mean_expr + se_expr),
                position = position_dodge(0.8), width = 0.2, linewidth = 0.4) +
  scale_fill_manual(values = zone_colors, name = "Zone") +
  labs(
    title = "Cancer Adjacent Zone",
    subtitle = "Active CAF signature",
    x = NULL,
    y = "Normalized Expression"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic"),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    legend.position = "none"
  )

# Panel 2: Cancer Distal genes (separate scale)
p_distal <- stats_by_zone %>%
  filter(gene %in% distal_avail) %>%
  mutate(gene = factor(gene, levels = distal_order)) %>%
  ggplot(aes(x = gene, y = mean_expr, fill = zone)) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = pmax(0, mean_expr - se_expr), ymax = mean_expr + se_expr),
                position = position_dodge(0.8), width = 0.2, linewidth = 0.4) +
  scale_fill_manual(values = zone_colors, name = "Zone") +
  labs(
    title = "Cancer Distal Zone",
    subtitle = "Quiescent CAF signature",
    x = NULL,
    y = "Normalized Expression"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic"),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    legend.position = "none"
  )

# Combine with shared legend
p_final <- (p_adjacent | p_distal) +
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(
    title = "Fibroblast Zone-Specific Gene Expression",
    subtitle = "Patient C Bx2 | nCount-normalized (mean ± SE)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "grey30")
    )
  )

# Create legend separately
p_legend <- stats_by_zone %>%
  head(2) %>%
  ggplot(aes(x = gene, y = mean_expr, fill = zone)) +
  geom_col() +
  scale_fill_manual(values = zone_colors, name = "Zone") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"))

legend <- cowplot::get_legend(p_legend)

# Final assembly
p_with_legend <- cowplot::plot_grid(
  p_final,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.08)
)

pdf(paste0(out_dir, "Fig_Final_DE_Genes_Curated.pdf"), width = 14, height = 6)
print(p_with_legend)
dev.off()
message("Saved: Fig_Final_DE_Genes_Curated.pdf")

# ==============================================================================
# ALTERNATIVE: WITHOUT COWPLOT (if not installed)
# ==============================================================================

# Simpler version using patchwork
p_simple <- (p_adjacent + theme(legend.position = "bottom")) | 
            (p_distal + theme(legend.position = "bottom"))

p_simple <- p_simple +
  plot_layout(widths = c(1.2, 1), guides = "collect") +
  plot_annotation(
    title = "Fibroblast Zone-Specific Gene Expression",
    subtitle = "Patient C Bx2 | nCount-normalized (mean ± SE)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "grey30")
    )
  ) &
  theme(legend.position = "bottom")

pdf(paste0(out_dir, "Fig_Final_DE_Genes_Curated_v2.pdf"), width = 14, height = 6)
print(p_simple)
dev.off()
message("Saved: Fig_Final_DE_Genes_Curated_v2.pdf")

# ==============================================================================
# SAVE STATS
# ==============================================================================

stats_curated <- stats_by_zone %>%
  filter(gene %in% c(adjacent_avail, distal_avail)) %>%
  select(zone, gene, category, mean_expr, se_expr, pct_expressing) %>%
  arrange(category, zone, desc(mean_expr))

write_csv(stats_curated, paste0(out_dir, "Final_DE_Genes_Stats.csv"))

# ==============================================================================
# SUMMARY
# ==============================================================================
message("\n")
message(strrep("=", 70))
message("FINAL CURATED GENE LIST")
message(strrep("=", 70))

message("\nCancer Adjacent - Active CAF (", length(adjacent_avail), " genes):")
message("  ", paste(adjacent_avail, collapse = ", "))

message("\nCancer Distal - Quiescent CAF (", length(distal_avail), " genes):")
message("  ", paste(distal_avail, collapse = ", "))

message("\nOutputs in: ", out_dir)
message("  - Fig_Final_DE_Genes_Curated.pdf")
message("  - Fig_Final_DE_Genes_Curated_v2.pdf (simpler version)")
message("  - Final_DE_Genes_Stats.csv")
