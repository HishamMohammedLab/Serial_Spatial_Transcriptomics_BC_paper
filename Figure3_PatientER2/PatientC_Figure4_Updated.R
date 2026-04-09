# ============================================================================
# Patient C Figure 4: Updated Figures with Clinical Context
# ============================================================================
# Clinical context:
# - ER+ breast cancer with different ESR1 mutations in Bx1 vs Bx2
# - Germline BRCA2 and MSH6 mutant
# - Bx2 has loss of: CDKN2A, CDKN2B, MTAP, BRCA2
# - Treatment: Tamoxifen + Everolimus + Paclitaxel (pre-Bx1)
#              Fulvestrant + Abemaciclib (between Bx1-Bx2)
# ============================================================================

library(Seurat)
library(tidyverse)
library(patchwork)

set.seed(42)

output_dir <- "PatientC_Analysis/Figure4_Updated/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading data...\n")
obj <- readRDS("data/CosMx_SMMART_345k_clean.rds")
patC <- subset(obj, Patient == "Patient_C")
expr_mat <- GetAssayData(patC, layer = "data")

# Get cancer cells
cancer_bx1 <- which(patC$Lineage == "Cancer" & patC$Timepoint == "Bx1")
cancer_bx2 <- which(patC$Lineage == "Cancer" & patC$Timepoint == "Bx2")

cat("Cancer cells: Bx1 =", length(cancer_bx1), ", Bx2 =", length(cancer_bx2), "\n")

# ============================================================================
# HELPER FUNCTION: Calculate expression stats
# ============================================================================

calc_gene_stats <- function(genes, cells_bx1, cells_bx2, expr_mat) {
  genes_avail <- genes[genes %in% rownames(expr_mat)]

  results <- data.frame()
  for (gene in genes_avail) {
    bx1_vals <- expr_mat[gene, cells_bx1]
    bx2_vals <- expr_mat[gene, cells_bx2]

    bx1_mean <- mean(bx1_vals, na.rm = TRUE)
    bx2_mean <- mean(bx2_vals, na.rm = TRUE)
    bx1_se <- sd(bx1_vals, na.rm = TRUE) / sqrt(length(bx1_vals))
    bx2_se <- sd(bx2_vals, na.rm = TRUE) / sqrt(length(bx2_vals))

    results <- rbind(results, data.frame(
      Gene = gene,
      Bx1_Mean = bx1_mean, Bx1_SE = bx1_se,
      Bx2_Mean = bx2_mean, Bx2_SE = bx2_se,
      FC = bx2_mean / bx1_mean,
      Log2FC = log2(bx2_mean / bx1_mean)
    ))
  }
  return(results)
}

# ============================================================================
# COLOR SCHEME
# ============================================================================

bx_colors <- c("Bx1" = "#7BCCC4", "Bx2" = "#2B8CBE")  # Teal gradient

# ============================================================================
# PANEL 1: LUMINAL IDENTITY (ESR1 mutation escape)
# ============================================================================

cat("\n=== Panel 1: Luminal Identity ===\n")

luminal_genes <- c("ESR1", "GATA3", "AGR2", "KRT8", "KRT18", "KRT19")
luminal_stats <- calc_gene_stats(luminal_genes, cancer_bx1, cancer_bx2, expr_mat)

luminal_stats$Gene <- factor(luminal_stats$Gene, levels = luminal_genes)

luminal_long <- luminal_stats %>%
  select(Gene, Bx1_Mean, Bx2_Mean, Bx1_SE, Bx2_SE) %>%
  pivot_longer(cols = c(Bx1_Mean, Bx2_Mean), names_to = "Timepoint", values_to = "Expression") %>%
  mutate(
    SE = ifelse(grepl("Bx1", Timepoint), Bx1_SE, Bx2_SE),
    Timepoint = gsub("_Mean", "", Timepoint)
  ) %>%
  select(Gene, Timepoint, Expression, SE)

p_luminal <- ggplot(luminal_long, aes(x = Gene, y = Expression, fill = Timepoint)) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = pmax(0, Expression - SE), ymax = Expression + SE),
                position = position_dodge(0.8), width = 0.2) +
  scale_fill_manual(values = bx_colors, name = "Biopsy") +
  labs(
    title = "Luminal Identity Strengthening",
    subtitle = "Despite Fulvestrant (ER degrader) | Different ESR1 mutations Bx1 vs Bx2",
    x = NULL, y = "Mean Expression (Cancer Cells)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    legend.position = "top"
  )

# Add FC annotations
fc_labels <- luminal_stats %>%
  mutate(label = sprintf("%.1fx", FC), y_pos = pmax(Bx1_Mean, Bx2_Mean) * 1.15)

p_luminal <- p_luminal +
  geom_text(data = fc_labels, aes(x = Gene, y = y_pos, label = label),
            inherit.aes = FALSE, size = 3, fontface = "bold", color = "darkred")

ggsave(file.path(output_dir, "Panel1_Luminal_Identity.pdf"), p_luminal, width = 8, height = 5)

# ============================================================================
# PANEL 2: HEAT SHOCK RESPONSE (Treatment stress adaptation)
# ============================================================================

cat("\n=== Panel 2: Heat Shock Response ===\n")

hsp_genes <- c("HSPB1", "HSP90AA1", "HSP90AB1", "HSP90B1", "HSPA1A", "HSPA1B")
hsp_stats <- calc_gene_stats(hsp_genes, cancer_bx1, cancer_bx2, expr_mat)

hsp_stats$Gene <- factor(hsp_stats$Gene, levels = hsp_genes)

hsp_long <- hsp_stats %>%
  select(Gene, Bx1_Mean, Bx2_Mean, Bx1_SE, Bx2_SE) %>%
  pivot_longer(cols = c(Bx1_Mean, Bx2_Mean), names_to = "Timepoint", values_to = "Expression") %>%
  mutate(
    SE = ifelse(grepl("Bx1", Timepoint), Bx1_SE, Bx2_SE),
    Timepoint = gsub("_Mean", "", Timepoint)
  ) %>%
  select(Gene, Timepoint, Expression, SE)

p_hsp <- ggplot(hsp_long, aes(x = Gene, y = Expression, fill = Timepoint)) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = pmax(0, Expression - SE), ymax = Expression + SE),
                position = position_dodge(0.8), width = 0.2) +
  scale_fill_manual(values = bx_colors, name = "Biopsy") +
  labs(
    title = "Heat Shock Protein Response",
    subtitle = "Treatment stress adaptation | HSP90 as potential therapeutic target",
    x = NULL, y = "Mean Expression (Cancer Cells)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    legend.position = "top"
  )

# Add FC annotations
fc_hsp <- hsp_stats %>%
  mutate(label = sprintf("%.1fx", FC), y_pos = pmax(Bx1_Mean, Bx2_Mean) * 1.15)

p_hsp <- p_hsp +
  geom_text(data = fc_hsp, aes(x = Gene, y = y_pos, label = label),
            inherit.aes = FALSE, size = 3, fontface = "bold",
            color = ifelse(fc_hsp$FC > 1.5, "darkred", "grey40"))

ggsave(file.path(output_dir, "Panel2_HeatShock_Response.pdf"), p_hsp, width = 8, height = 5)

# ============================================================================
# PANEL 3: PROLIFERATION (Despite CDK4/6 inhibitor)
# ============================================================================

cat("\n=== Panel 3: Proliferation ===\n")

prolif_genes <- c("MKI67", "TOP2A", "PCNA", "CCND1", "CCNE1", "RB1")
prolif_stats <- calc_gene_stats(prolif_genes, cancer_bx1, cancer_bx2, expr_mat)

prolif_stats$Gene <- factor(prolif_stats$Gene, levels = prolif_genes)

prolif_long <- prolif_stats %>%
  select(Gene, Bx1_Mean, Bx2_Mean, Bx1_SE, Bx2_SE) %>%
  pivot_longer(cols = c(Bx1_Mean, Bx2_Mean), names_to = "Timepoint", values_to = "Expression") %>%
  mutate(
    SE = ifelse(grepl("Bx1", Timepoint), Bx1_SE, Bx2_SE),
    Timepoint = gsub("_Mean", "", Timepoint)
  ) %>%
  select(Gene, Timepoint, Expression, SE)

p_prolif <- ggplot(prolif_long, aes(x = Gene, y = Expression, fill = Timepoint)) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = pmax(0, Expression - SE), ymax = Expression + SE),
                position = position_dodge(0.8), width = 0.2) +
  scale_fill_manual(values = bx_colors, name = "Biopsy") +
  labs(
    title = "Proliferation Maintained Despite CDK4/6 Inhibitor",
    subtitle = "Abemaciclib between biopsies | CDKN2A loss may enable bypass",
    x = NULL, y = "Mean Expression (Cancer Cells)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    legend.position = "top"
  )

ggsave(file.path(output_dir, "Panel3_Proliferation.pdf"), p_prolif, width = 8, height = 5)

# ============================================================================
# PANEL 4: INTERFERON SUPPRESSION (Immune evasion)
# ============================================================================

cat("\n=== Panel 4: Interferon Response ===\n")

ifn_genes <- c("IFITM1", "MX1", "STAT1", "HLA-A", "HLA-B", "B2M", "OAS1", "OAS2")
ifn_stats <- calc_gene_stats(ifn_genes, cancer_bx1, cancer_bx2, expr_mat)

ifn_stats$Gene <- factor(ifn_stats$Gene, levels = ifn_genes)

ifn_long <- ifn_stats %>%
  select(Gene, Bx1_Mean, Bx2_Mean, Bx1_SE, Bx2_SE) %>%
  pivot_longer(cols = c(Bx1_Mean, Bx2_Mean), names_to = "Timepoint", values_to = "Expression") %>%
  mutate(
    SE = ifelse(grepl("Bx1", Timepoint), Bx1_SE, Bx2_SE),
    Timepoint = gsub("_Mean", "", Timepoint)
  ) %>%
  select(Gene, Timepoint, Expression, SE)

p_ifn <- ggplot(ifn_long, aes(x = Gene, y = Expression, fill = Timepoint)) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = pmax(0, Expression - SE), ymax = Expression + SE),
                position = position_dodge(0.8), width = 0.2) +
  scale_fill_manual(values = bx_colors, name = "Biopsy") +
  labs(
    title = "Interferon Response Suppression",
    subtitle = "Reduced antigen presentation and IFN signaling | Immune evasion",
    x = NULL, y = "Mean Expression (Cancer Cells)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    legend.position = "top"
  )

ggsave(file.path(output_dir, "Panel4_Interferon_Down.pdf"), p_ifn, width = 8, height = 5)

# ============================================================================
# PANEL 5: DNA DAMAGE RESPONSE (BRCA2 context)
# ============================================================================

cat("\n=== Panel 5: DNA Damage Response ===\n")

ddr_genes <- c("BRCA1", "ATM", "ATR", "RAD51", "CHEK1", "CHEK2", "PARP1")
ddr_stats <- calc_gene_stats(ddr_genes, cancer_bx1, cancer_bx2, expr_mat)

ddr_stats$Gene <- factor(ddr_stats$Gene, levels = ddr_genes)

ddr_long <- ddr_stats %>%
  select(Gene, Bx1_Mean, Bx2_Mean, Bx1_SE, Bx2_SE) %>%
  pivot_longer(cols = c(Bx1_Mean, Bx2_Mean), names_to = "Timepoint", values_to = "Expression") %>%
  mutate(
    SE = ifelse(grepl("Bx1", Timepoint), Bx1_SE, Bx2_SE),
    Timepoint = gsub("_Mean", "", Timepoint)
  ) %>%
  select(Gene, Timepoint, Expression, SE)

p_ddr <- ggplot(ddr_long, aes(x = Gene, y = Expression, fill = Timepoint)) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = pmax(0, Expression - SE), ymax = Expression + SE),
                position = position_dodge(0.8), width = 0.2) +
  scale_fill_manual(values = bx_colors, name = "Biopsy") +
  labs(
    title = "DNA Damage Response Compensation",
    subtitle = "Germline BRCA2 + somatic loss | PARP1 down suggests PARP inhibitor sensitivity",
    x = NULL, y = "Mean Expression (Cancer Cells)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    legend.position = "top"
  )

# Highlight PARP1 (down) differently
fc_ddr <- ddr_stats %>%
  mutate(label = sprintf("%.2fx", FC), y_pos = pmax(Bx1_Mean, Bx2_Mean) * 1.15)

p_ddr <- p_ddr +
  geom_text(data = fc_ddr, aes(x = Gene, y = y_pos, label = label),
            inherit.aes = FALSE, size = 3, fontface = "bold",
            color = ifelse(fc_ddr$FC < 1, "darkblue", "darkred"))

ggsave(file.path(output_dir, "Panel5_DDR_Response.pdf"), p_ddr, width = 8, height = 5)

# ============================================================================
# PANEL 6: ECM - CLARIFIED (Cancer vs Fibroblast)
# ============================================================================

cat("\n=== Panel 6: ECM by Cell Type ===\n")

# Get fibroblast cells
fibro_bx1 <- which(patC$Lineage == "Fibroblast" & patC$Timepoint == "Bx1")
fibro_bx2 <- which(patC$Lineage == "Fibroblast" & patC$Timepoint == "Bx2")

ecm_genes <- c("COL1A1", "COL1A2", "COL3A1", "FN1")

# Calculate for both cell types
ecm_cancer <- calc_gene_stats(ecm_genes, cancer_bx1, cancer_bx2, expr_mat) %>%
  mutate(CellType = "Cancer")
ecm_fibro <- calc_gene_stats(ecm_genes, fibro_bx1, fibro_bx2, expr_mat) %>%
  mutate(CellType = "Fibroblast")

ecm_combined <- bind_rows(ecm_cancer, ecm_fibro)
ecm_combined$Gene <- factor(ecm_combined$Gene, levels = ecm_genes)

ecm_long <- ecm_combined %>%
  select(Gene, CellType, Bx1_Mean, Bx2_Mean) %>%
  pivot_longer(cols = c(Bx1_Mean, Bx2_Mean), names_to = "Timepoint", values_to = "Expression") %>%
  mutate(Timepoint = gsub("_Mean", "", Timepoint))

# Create grouped barplot
ecm_long$Group <- paste(ecm_long$CellType, ecm_long$Timepoint, sep = "_")
ecm_long$Group <- factor(ecm_long$Group,
                          levels = c("Cancer_Bx1", "Cancer_Bx2", "Fibroblast_Bx1", "Fibroblast_Bx2"))

group_colors <- c("Cancer_Bx1" = "#FEE5D9", "Cancer_Bx2" = "#FB6A4A",
                  "Fibroblast_Bx1" = "#DEEBF7", "Fibroblast_Bx2" = "#3182BD")

p_ecm <- ggplot(ecm_long, aes(x = Gene, y = Expression, fill = Group)) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = group_colors,
                    name = "Cell Type / Biopsy",
                    labels = c("Cancer Bx1", "Cancer Bx2", "Fibroblast Bx1", "Fibroblast Bx2")) +
  labs(
    title = "ECM Gene Expression by Cell Type",
    subtitle = "Fibroblasts have 10-20x higher ECM | Fibroblast COL1A1 decreases Bx1->Bx2",
    x = NULL, y = "Mean Expression"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    legend.position = "top"
  ) +
  guides(fill = guide_legend(nrow = 1))

ggsave(file.path(output_dir, "Panel6_ECM_ByCellType.pdf"), p_ecm, width = 9, height = 5)

# ============================================================================
# COMBINED SUMMARY FIGURE
# ============================================================================

cat("\n=== Creating Combined Figure ===\n")

# Summary heatmap of fold changes
summary_data <- bind_rows(
  luminal_stats %>% mutate(Category = "Luminal"),
  hsp_stats %>% mutate(Category = "Heat Shock"),
  prolif_stats %>% mutate(Category = "Proliferation"),
  ifn_stats %>% mutate(Category = "Interferon"),
  ddr_stats %>% mutate(Category = "DDR")
) %>%
  select(Gene, Log2FC, Category) %>%
  filter(!is.na(Log2FC) & is.finite(Log2FC))

# Order by category and Log2FC
summary_data <- summary_data %>%
  arrange(Category, desc(Log2FC)) %>%
  mutate(Gene = factor(Gene, levels = unique(Gene)))

p_summary <- ggplot(summary_data, aes(x = Gene, y = Log2FC, fill = Category)) +
  geom_col(color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_brewer(palette = "Set2") +
  coord_flip() +
  labs(
    title = "Patient C: Gene Expression Changes (Bx1 to Bx2)",
    subtitle = "Cancer cells only | Log2 fold change",
    x = NULL, y = "Log2 Fold Change (Bx2/Bx1)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "Summary_FoldChanges.pdf"), p_summary, width = 10, height = 8)

# ============================================================================
# MAIN COMBINED PANEL
# ============================================================================

# Top row: Luminal + HSP
# Middle row: Proliferation + IFN
# Bottom row: DDR + ECM

p_combined <- (p_luminal + p_hsp) /
              (p_prolif + p_ifn) /
              (p_ddr + p_ecm) +
  plot_annotation(
    title = "Patient C: Treatment Response and Adaptation",
    subtitle = "ER+ breast cancer | ESR1 mutations | Germline BRCA2/MSH6 | Fulvestrant + Abemaciclib",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 11, color = "grey40", hjust = 0.5)
    )
  )

ggsave(file.path(output_dir, "Figure4_Combined_AllPanels.pdf"), p_combined,
       width = 16, height = 18)

# ============================================================================
# SAVE STATISTICS
# ============================================================================

all_stats <- bind_rows(
  luminal_stats %>% mutate(Category = "Luminal"),
  hsp_stats %>% mutate(Category = "Heat_Shock"),
  prolif_stats %>% mutate(Category = "Proliferation"),
  ifn_stats %>% mutate(Category = "Interferon"),
  ddr_stats %>% mutate(Category = "DDR")
)

write_csv(all_stats, file.path(output_dir, "All_Gene_Statistics.csv"))

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("PATIENT C FIGURE 4 - UPDATED FIGURES COMPLETE\n")
cat(strrep("=", 70), "\n")

cat("\nOutput directory:", output_dir, "\n")
cat("\nFigures created:\n")
cat("  1. Panel1_Luminal_Identity.pdf - ESR1 mutation escape\n")
cat("  2. Panel2_HeatShock_Response.pdf - Treatment adaptation\n")
cat("  3. Panel3_Proliferation.pdf - Despite CDK4/6i\n")
cat("  4. Panel4_Interferon_Down.pdf - Immune evasion\n")
cat("  5. Panel5_DDR_Response.pdf - BRCA2 compensation\n")
cat("  6. Panel6_ECM_ByCellType.pdf - Clarified ECM\n")
cat("  7. Summary_FoldChanges.pdf - All changes summary\n")
cat("  8. Figure4_Combined_AllPanels.pdf - Full figure\n")

cat("\nKey narratives:\n")
cat("  - Luminal identity UP despite Fulvestrant (ESR1 mutation escape)\n")
cat("  - Heat shock response (HSPB1 4x) - treatment stress adaptation\n")
cat("  - Proliferation maintained despite Abemaciclib (CDKN2A loss?)\n")
cat("  - Interferon/HLA down - immune evasion\n")
cat("  - DDR compensation for BRCA2 loss, PARP1 down\n")
