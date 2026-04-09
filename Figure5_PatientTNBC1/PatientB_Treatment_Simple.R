# ==============================================================================
# Patient B: Simplified Treatment Response Bar Plots
# Horizontal bars - negative = down in Bx4, positive = up in Bx4
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

cat("Loading data...\n")
obj <- subset(readRDS("data/CosMx_SMMART_345k_clean.rds"), broad_lineage == "Cancer")
patB <- subset(obj, subset = Patient == "Patient_B")

# ==============================================================================
# Define Genes by Treatment
# ==============================================================================

treatment_genes <- tribble(
  ~Gene,       ~Treatment,     ~Expected,   ~CNV,
  # OLAPARIB
  "PARP1",     "Olaparib",     "↓",         FALSE,
  "BRCA1",     "Olaparib",     "↑",         FALSE,
  "RAD51",     "Olaparib",     "↑",         FALSE,
  "ATR",       "Olaparib",     "↑",         TRUE,
  "ATM",       "Olaparib",     "↑",         FALSE,
  "CHEK1",     "Olaparib",     "↑",         FALSE,
  "CHEK2",     "Olaparib",     "↑",         FALSE,
  # DURVALUMAB
  "CD274",     "Durvalumab",   "~",         FALSE,
  "B2M",       "Durvalumab",   "↓",         FALSE,
  "HLA-A",     "Durvalumab",   "↓",         FALSE,
  "HLA-B",     "Durvalumab",   "↓",         FALSE,
  "TAP1",      "Durvalumab",   "↓",         FALSE,
  "CD47",      "Durvalumab",   "↑",         FALSE,
  # SACITUZUMAB
  "TACSTD2",   "Sacituzumab",  "↓",         FALSE,
  "TOP2A",     "Sacituzumab",  "~",         FALSE,
  # PI3K (CNV)
  "AKT1",      "PI3K/CNV",     "↑",         TRUE,
  "MTOR",      "PI3K/CNV",     "↑",         TRUE
)

# Filter to available
available <- treatment_genes$Gene[treatment_genes$Gene %in% rownames(patB)]
treatment_genes <- treatment_genes %>% filter(Gene %in% available)

# ==============================================================================
# Calculate Log2FC
# ==============================================================================

expr_mat <- GetAssayData(patB, layer = "data")
expr_df <- as.data.frame(t(as.matrix(expr_mat[available, , drop = FALSE])))
expr_df$Timepoint <- patB$Timepoint

stats_list <- lapply(available, function(g) {
  bx1 <- expr_df[[g]][expr_df$Timepoint == "Bx1"]
  bx4 <- expr_df[[g]][expr_df$Timepoint == "Bx4"]
  data.frame(
    Gene = g,
    Log2FC = log2((mean(bx4, na.rm = TRUE) + 0.01) / (mean(bx1, na.rm = TRUE) + 0.01))
  )
})

stats_df <- do.call(rbind, stats_list)
plot_data <- merge(stats_df, treatment_genes, by = "Gene")

# Create labels with expected direction
plot_data <- plot_data %>%
  mutate(
    Gene_Label = ifelse(CNV,
                        paste0(Gene, " [CNV] (exp: ", Expected, ")"),
                        paste0(Gene, " (exp: ", Expected, ")")),
    Treatment = factor(Treatment, levels = c("Olaparib", "Durvalumab", "Sacituzumab", "PI3K/CNV"))
  )

# ==============================================================================
# Create Simple Horizontal Bar Plot
# ==============================================================================

# Color by direction
plot_data$Direction <- ifelse(plot_data$Log2FC > 0, "Up", "Down")

# Order genes within each treatment by Log2FC
plot_data <- plot_data %>%
  arrange(Treatment, Log2FC) %>%
  mutate(Gene_Label = factor(Gene_Label, levels = unique(Gene_Label)))

p_simple <- ggplot(plot_data, aes(x = Log2FC, y = Gene_Label, fill = Direction)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  facet_grid(Treatment ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("Up" = "#D73027", "Down" = "#4575B4"),
                    name = "Bx4 vs Bx1") +
  scale_x_continuous(limits = c(-1.5, 1.5), breaks = seq(-1.5, 1.5, 0.5)) +
  labs(
    title = "Patient B: Treatment Target Expression Changes",
    subtitle = "Bx1 (Lymph Node) → Bx4 (Soft Tissue) | exp: expected direction if resistant",
    x = "Log2 Fold Change (← Down in Bx4 | Up in Bx4 →)",
    y = NULL,
    caption = "CNV gains: ATR, PIK3CB detected in Bx4"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 10),
    plot.caption = element_text(color = "#D73027", size = 9),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey90", color = NA),
    axis.text.y = element_text(size = 10),
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "PatientB_Treatment_Simple.pdf"),
       p_simple, width = 9, height = 10)

# ==============================================================================
# Even simpler - one plot per treatment
# ==============================================================================

make_simple_plot <- function(tx_name, tx_color) {
  tx_data <- plot_data %>% filter(Treatment == tx_name) %>%
    arrange(Log2FC) %>%
    mutate(Gene_Label = factor(Gene_Label, levels = Gene_Label))

  ggplot(tx_data, aes(x = Log2FC, y = Gene_Label, fill = Direction)) +
    geom_col(width = 0.6) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    scale_fill_manual(values = c("Up" = "#D73027", "Down" = "#4575B4"), guide = "none") +
    scale_x_continuous(limits = c(-1.5, 1.5)) +
    labs(title = tx_name, x = "Log2FC", y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

p1 <- make_simple_plot("Olaparib", "#2E7D32")
p2 <- make_simple_plot("Durvalumab", "#7B1FA2")
p3 <- make_simple_plot("Sacituzumab", "#E65100")
p4 <- make_simple_plot("PI3K/CNV", "#C62828")

p_grid <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "Patient B: Treatment Response",
    subtitle = "Gene (exp: expected if resistant) | [CNV] = copy number gain\nBlue = down in Bx4 | Red = up in Bx4",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(color = "grey40", size = 10, hjust = 0.5)
    )
  )

ggsave(file.path(output_dir, "PatientB_Treatment_Grid.pdf"),
       p_grid, width = 12, height = 8)

cat("\nSaved:\n")
cat("  - PatientB_Treatment_Simple.pdf (single figure)\n")
cat("  - PatientB_Treatment_Grid.pdf (grid layout)\n")
