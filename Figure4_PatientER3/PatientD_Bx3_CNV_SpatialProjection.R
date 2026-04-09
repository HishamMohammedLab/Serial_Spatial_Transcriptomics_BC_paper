#!/usr/bin/env Rscript
# =============================================================================
# Patient D Bx3 - Spatial Projection of CNV Probes
# One plot per probe showing spatial distribution of scores
# =============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# =============================================================================
# CONFIGURATION
# =============================================================================

SEURAT_PATH <- "data/CosMx_SMMART_345k_clean.rds"
GENE_ORDER_PATH <- "supplementary_input_data/CNV/gene_order.txt"
OUTPUT_DIR <- "CNV_Analysis_PatientD/Spatial_Probes"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=============================================================\n")
cat("PATIENT D BX3 - SPATIAL CNV PROBE PROJECTION\n")
cat("=============================================================\n\n")

# =============================================================================
# DEFINE PROBES
# =============================================================================

PROBES <- list(
  # GAIN probes
  "11q13_amp" = list(chr = "chr11", start = 69000000, end = 77000000, type = "GAIN",
                     description = "CCND1 amplicon (Report: ~7 copies)"),
  "17q22_23_amp" = list(chr = "chr17", start = 58000000, end = 75000000, type = "GAIN",
                        description = "RNF43/PPM1D/BRIP1 (Report: ~6-9 copies)"),
  "8q_gain" = list(chr = "chr8", start = 127000000, end = 145000000, type = "GAIN",
                   description = "MYC region (Literature)"),
  "1q_gain" = list(chr = "chr1", start = 142000000, end = 250000000, type = "GAIN",
                   description = "1q arm gain (Literature: >50% breast cancers)"),
  "7q34_BRAF" = list(chr = "chr7", start = 140000000, end = 145000000, type = "GAIN",
                     description = "BRAF amplification (Report: ~6-15 copies)"),

  # LOSS probes
  "22q12_CHEK2" = list(chr = "chr22", start = 28000000, end = 32000000, type = "LOSS",
                       description = "CHEK2 loss (Report: 0.4 copies, germline variant)"),
  "17p13_TP53" = list(chr = "chr17", start = 0, end = 15000000, type = "LOSS",
                      description = "TP53 loss (Report: 0.4 copies)"),
  "11q24_CHEK1" = list(chr = "chr11", start = 125000000, end = 135000000, type = "LOSS",
                       description = "CHEK1 loss (Report: 0.1 copies)"),
  "9p21_CDKN2A" = list(chr = "chr9", start = 20000000, end = 25000000, type = "LOSS",
                       description = "CDKN2A loss (Report: 0.1 copies)"),
  "16q_loss" = list(chr = "chr16", start = 46000000, end = 90000000, type = "LOSS",
                    description = "16q loss - CDH1 region (Literature)"),
  "8p_loss" = list(chr = "chr8", start = 0, end = 45000000, type = "LOSS",
                   description = "8p loss (Literature)"),
  "3p25_FANCD2" = list(chr = "chr3", start = 0, end = 15000000, type = "LOSS",
                       description = "FANCD2 loss (Report: 0.4 copies)")
)

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data...\n")

gene_pos <- read.table(GENE_ORDER_PATH, header = FALSE, sep = "\t",
                       col.names = c("Gene", "Chr", "Start", "End"))

seurat_obj <- readRDS(SEURAT_PATH)

# Filter to Patient D Bx3
cells_patD_bx3 <- colnames(seurat_obj)[seurat_obj$Patient == "Patient_D" &
                                        seurat_obj$Timepoint == "Bx3"]
seurat_sub <- subset(seurat_obj, cells = cells_patD_bx3)
cat("Patient D Bx3 cells:", ncol(seurat_sub), "\n")

# Get expression and coordinates
expr_mat <- GetAssayData(seurat_sub, assay = "RNA", layer = "data")
coords <- data.frame(
  Cell = colnames(seurat_sub),
  x = seurat_sub$x_global_px,
  y = seurat_sub$y_global_px,
  Lineage = seurat_sub$Lineage,
  CellType = ifelse(seurat_sub$Lineage == "Cancer", "Cancer", "Non-Cancer")
)

cat("Cells with coordinates:", sum(!is.na(coords$x)), "\n\n")

# =============================================================================
# MAP GENES TO PROBES
# =============================================================================

cat("Mapping genes to probes...\n")

probe_genes <- list()
for (probe_name in names(PROBES)) {
  probe <- PROBES[[probe_name]]
  genes <- gene_pos %>%
    filter(Chr == probe$chr, Start >= probe$start, Start <= probe$end) %>%
    pull(Gene)
  genes <- intersect(genes, rownames(expr_mat))
  probe_genes[[probe_name]] <- genes
  cat(sprintf("  %s: %d genes\n", probe_name, length(genes)))
}

# =============================================================================
# CALCULATE PROBE SCORES
# =============================================================================

cat("\nCalculating probe scores...\n")

# Reference: non-cancer cells
lineages <- seurat_sub$Lineage
non_cancer_idx <- which(lineages != "Cancer")

# Calculate scores for each probe
probe_scores <- data.frame(Cell = colnames(expr_mat))

for (probe_name in names(PROBES)) {
  genes <- probe_genes[[probe_name]]
  probe_type <- PROBES[[probe_name]]$type

  if (length(genes) == 0) {
    probe_scores[[probe_name]] <- NA
    next
  }

  # Calculate mean expression
  if (length(genes) == 1) {
    probe_expr <- as.numeric(expr_mat[genes, ])
  } else {
    probe_expr <- colMeans(as.matrix(expr_mat[genes, , drop = FALSE]))
  }

  # Z-score vs non-cancer reference
  ref_mean <- mean(probe_expr[non_cancer_idx], na.rm = TRUE)
  ref_sd <- sd(probe_expr[non_cancer_idx], na.rm = TRUE)
  if (ref_sd < 0.01) ref_sd <- 0.01

  z_scores <- (probe_expr - ref_mean) / ref_sd

  # For LOSS probes, invert so positive = more loss
  if (probe_type == "LOSS") {
    z_scores <- -z_scores
  }

  probe_scores[[probe_name]] <- z_scores
}

# Merge with coordinates
plot_data <- merge(coords, probe_scores, by = "Cell")
cat("Cells for plotting:", nrow(plot_data), "\n\n")

# =============================================================================
# CREATE INDIVIDUAL SPATIAL PLOTS
# =============================================================================

cat("Creating spatial plots...\n\n")

for (probe_name in names(PROBES)) {
  probe <- PROBES[[probe_name]]

  if (all(is.na(plot_data[[probe_name]]))) {
    cat(sprintf("  Skipping %s (no data)\n", probe_name))
    next
  }

  cat(sprintf("  Plotting %s...\n", probe_name))

  # Prepare data - sort so high values plot on top
  df <- plot_data[!is.na(plot_data[[probe_name]]) & !is.na(plot_data$x), ]
  df$Score <- df[[probe_name]]
  df <- df[order(df$Score), ]  # Low values first

  # Cap extreme values
  df$Score[df$Score > 3] <- 3
  df$Score[df$Score < -3] <- -3

  # Color scale based on probe type
  if (probe$type == "GAIN") {
    colors <- c("#2166AC", "#F7F7F7", "#B2182B")  # Blue-White-Red
    midpoint <- 0
  } else {
    colors <- c("#B2182B", "#F7F7F7", "#2166AC")  # Red-White-Blue (inverted for loss)
    midpoint <- 0
  }

  # Create plot
  p <- ggplot(df, aes(x = x, y = y, color = Score)) +
    geom_point(size = 0.3, alpha = 0.7) +
    scale_color_gradient2(
      low = colors[1], mid = colors[2], high = colors[3],
      midpoint = midpoint,
      limits = c(-3, 3),
      name = "Z-score",
      oob = scales::squish
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "right",
      plot.margin = margin(10, 10, 10, 10)
    ) +
    coord_fixed() +
    labs(
      title = sprintf("%s (%s)", probe_name, probe$type),
      subtitle = probe$description
    )

  ggsave(
    file.path(OUTPUT_DIR, sprintf("PatientD_Bx3_Spatial_%s.pdf", probe_name)),
    p, width = 10, height = 8
  )
}

# =============================================================================
# CREATE COMBINED FIGURE - ALL PROBES
# =============================================================================

cat("\nCreating combined figure...\n")

# Create list of plots
plot_list <- list()

for (probe_name in names(PROBES)) {
  probe <- PROBES[[probe_name]]

  if (all(is.na(plot_data[[probe_name]]))) next

  df <- plot_data[!is.na(plot_data[[probe_name]]) & !is.na(plot_data$x), ]
  df$Score <- df[[probe_name]]
  df <- df[order(df$Score), ]
  df$Score[df$Score > 3] <- 3
  df$Score[df$Score < -3] <- -3

  p <- ggplot(df, aes(x = x, y = y, color = Score)) +
    geom_point(size = 0.1, alpha = 0.6) +
    scale_color_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-3, 3),
      oob = scales::squish
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 9, face = "bold"),
      legend.position = "none"
    ) +
    coord_fixed() +
    labs(title = probe_name)

  plot_list[[probe_name]] <- p
}

# Arrange in grid
combined <- wrap_plots(plot_list, ncol = 4) +
  plot_annotation(
    title = "Patient D Bx3 - Spatial CNV Probe Scores",
    subtitle = "Blue = Low, Red = High (z-scores vs non-cancer reference)",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

ggsave(
  file.path(OUTPUT_DIR, "PatientD_Bx3_Spatial_AllProbes_Grid.pdf"),
  combined, width = 20, height = 15
)
cat("Saved: PatientD_Bx3_Spatial_AllProbes_Grid.pdf\n")

# =============================================================================
# CREATE SIDE-BY-SIDE WITH CELL TYPE OVERLAY
# =============================================================================

cat("\nCreating probe + cell type overlay figures...\n")

for (probe_name in names(PROBES)[1:5]) {  # Top 5 probes
  probe <- PROBES[[probe_name]]

  if (all(is.na(plot_data[[probe_name]]))) next

  df <- plot_data[!is.na(plot_data[[probe_name]]) & !is.na(plot_data$x), ]
  df$Score <- df[[probe_name]]
  df$Score[df$Score > 3] <- 3
  df$Score[df$Score < -3] <- -3

  # Plot 1: CNV score
  p1 <- ggplot(df[order(df$Score), ], aes(x = x, y = y, color = Score)) +
    geom_point(size = 0.3, alpha = 0.7) +
    scale_color_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-3, 3),
      name = "Z-score",
      oob = scales::squish
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "bottom"
    ) +
    coord_fixed() +
    labs(title = sprintf("%s Score", probe_name))

  # Plot 2: Cell type
  p2 <- ggplot(df, aes(x = x, y = y, color = CellType)) +
    geom_point(size = 0.3, alpha = 0.7) +
    scale_color_manual(values = c("Cancer" = "#E41A1C", "Non-Cancer" = "#377EB8")) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "bottom"
    ) +
    coord_fixed() +
    labs(title = "Cell Type")

  combined <- p1 + p2 +
    plot_annotation(
      title = sprintf("%s - %s", probe_name, probe$description),
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )

  ggsave(
    file.path(OUTPUT_DIR, sprintf("PatientD_Bx3_Spatial_%s_withCellType.pdf", probe_name)),
    combined, width = 16, height = 7
  )
}

# =============================================================================
# CREATE GAIN vs LOSS COMPARISON
# =============================================================================

cat("\nCreating GAIN vs LOSS comparison...\n")

# Best GAIN probes
gain_probes <- c("11q13_amp", "8q_gain", "17q22_23_amp")
loss_probes <- c("22q12_CHEK2", "16q_loss", "8p_loss")

gain_plots <- list()
for (probe_name in gain_probes) {
  df <- plot_data[!is.na(plot_data[[probe_name]]) & !is.na(plot_data$x), ]
  df$Score <- df[[probe_name]]
  df <- df[order(df$Score), ]
  df$Score[df$Score > 3] <- 3
  df$Score[df$Score < -3] <- -3

  p <- ggplot(df, aes(x = x, y = y, color = Score)) +
    geom_point(size = 0.2, alpha = 0.6) +
    scale_color_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-3, 3),
      oob = scales::squish
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      legend.position = "none"
    ) +
    coord_fixed() +
    labs(title = probe_name)

  gain_plots[[probe_name]] <- p
}

loss_plots <- list()
for (probe_name in loss_probes) {
  df <- plot_data[!is.na(plot_data[[probe_name]]) & !is.na(plot_data$x), ]
  df$Score <- df[[probe_name]]
  df <- df[order(df$Score), ]
  df$Score[df$Score > 3] <- 3
  df$Score[df$Score < -3] <- -3

  p <- ggplot(df, aes(x = x, y = y, color = Score)) +
    geom_point(size = 0.2, alpha = 0.6) +
    scale_color_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-3, 3),
      oob = scales::squish
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      legend.position = "none"
    ) +
    coord_fixed() +
    labs(title = probe_name)

  loss_plots[[probe_name]] <- p
}

# Combine
gain_combined <- wrap_plots(gain_plots, ncol = 3) +
  plot_annotation(title = "GAIN Probes (Red = Higher in Cancer)",
                  theme = theme(plot.title = element_text(face = "bold", hjust = 0.5)))

loss_combined <- wrap_plots(loss_plots, ncol = 3) +
  plot_annotation(title = "LOSS Probes (Blue = Lower in Cancer)",
                  theme = theme(plot.title = element_text(face = "bold", hjust = 0.5)))

final <- gain_combined / loss_combined +
  plot_annotation(
    title = "Patient D Bx3 - Top CNV Probes Spatial Distribution",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

ggsave(
  file.path(OUTPUT_DIR, "PatientD_Bx3_Spatial_TopProbes_GainVsLoss.pdf"),
  final, width = 18, height = 12
)
cat("Saved: PatientD_Bx3_Spatial_TopProbes_GainVsLoss.pdf\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=============================================================\n")
cat("OUTPUT FILES\n")
cat("=============================================================\n\n")

cat("Individual probe plots:\n")
for (probe_name in names(PROBES)) {
  cat(sprintf("  - PatientD_Bx3_Spatial_%s.pdf\n", probe_name))
}

cat("\nCombined figures:\n")
cat("  - PatientD_Bx3_Spatial_AllProbes_Grid.pdf\n")
cat("  - PatientD_Bx3_Spatial_TopProbes_GainVsLoss.pdf\n")

cat("\nProbe + Cell Type overlays:\n")
for (probe_name in names(PROBES)[1:5]) {
  cat(sprintf("  - PatientD_Bx3_Spatial_%s_withCellType.pdf\n", probe_name))
}

cat(sprintf("\nOutput directory: %s\n", OUTPUT_DIR))

cat("\n=============================================================\n")
cat("DONE\n")
cat("=============================================================\n")
