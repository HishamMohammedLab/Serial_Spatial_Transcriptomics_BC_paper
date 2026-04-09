# =============================================================================
# Patient D Bx3 - Spatial CNV Probe Projection with POLYGONS
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
# R1124_322091 is the RNA polygon file for Patient D Bx3 (verified by coordinate matching)
POLYGON_PATH <- "supplementary_input_data/Polygons/RNA/R1124_322091-polygons.csv"
OUTPUT_DIR <- "CNV_Analysis_PatientD/Spatial_Polygons"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=============================================================\n")
cat("PATIENT D BX3 - SPATIAL CNV POLYGONS\n")
cat("=============================================================\n\n")

# =============================================================================
# DEFINE PROBES (Top performing ones)
# =============================================================================

PROBES <- list(
  # GAIN probes
  "11q13_amp" = list(chr = "chr11", start = 69000000, end = 77000000, type = "GAIN",
                     description = "CCND1 amplicon (~7 copies)"),
  "8q_gain" = list(chr = "chr8", start = 127000000, end = 145000000, type = "GAIN",
                   description = "MYC region"),
  "17q22_23_amp" = list(chr = "chr17", start = 58000000, end = 75000000, type = "GAIN",
                        description = "RNF43/PPM1D/BRIP1 (~6-9 copies)"),
  "1q_gain" = list(chr = "chr1", start = 142000000, end = 250000000, type = "GAIN",
                   description = "1q arm gain (>50% breast cancers)"),
  "7q34_BRAF" = list(chr = "chr7", start = 140000000, end = 145000000, type = "GAIN",
                     description = "BRAF amplification (~6-15 copies)"),

  # LOSS probes
  "22q12_CHEK2" = list(chr = "chr22", start = 28000000, end = 32000000, type = "LOSS",
                       description = "CHEK2 loss (0.4 copies, germline variant)"),
  "16q_loss" = list(chr = "chr16", start = 46000000, end = 90000000, type = "LOSS",
                    description = "16q loss - CDH1 region"),
  "8p_loss" = list(chr = "chr8", start = 0, end = 45000000, type = "LOSS",
                   description = "8p loss"),
  "17p13_TP53" = list(chr = "chr17", start = 0, end = 15000000, type = "LOSS",
                      description = "TP53 loss (0.4 copies)"),
  "9p21_CDKN2A" = list(chr = "chr9", start = 20000000, end = 25000000, type = "LOSS",
                       description = "CDKN2A loss (0.1 copies)"),
  "11q24_CHEK1" = list(chr = "chr11", start = 125000000, end = 135000000, type = "LOSS",
                       description = "CHEK1 loss (0.1 copies)"),
  "3p25_FANCD2" = list(chr = "chr3", start = 0, end = 15000000, type = "LOSS",
                       description = "FANCD2 loss (0.4 copies)")
)

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading polygon data...\n")
polygons <- read.csv(POLYGON_PATH)
cat("Polygon vertices:", nrow(polygons), "\n")

# Create cell key for polygon data
polygons$cell_key <- paste(polygons$fov, polygons$cellID, sep = "_")
unique_poly_cells <- unique(polygons$cell_key)
cat("Unique cells in polygon file:", length(unique_poly_cells), "\n")

cat("\nLoading gene positions...\n")
gene_pos <- read.table(GENE_ORDER_PATH, header = FALSE, sep = "\t",
                       col.names = c("Gene", "Chr", "Start", "End"))

cat("Loading Seurat object...\n")
seurat_obj <- readRDS(SEURAT_PATH)

# Filter to Patient D Bx3
cells_patD_bx3 <- colnames(seurat_obj)[seurat_obj$Patient == "Patient_D" &
                                        seurat_obj$Timepoint == "Bx3"]
seurat_sub <- subset(seurat_obj, cells = cells_patD_bx3)
cat("Patient D Bx3 cells:", ncol(seurat_sub), "\n")

# Parse Seurat cell names: c_5_fov_cellID
cell_names <- colnames(seurat_sub)
cell_parts <- strsplit(cell_names, "_")
cell_mapping <- data.frame(
  Cell = cell_names,
  fov = as.integer(sapply(cell_parts, function(x) x[3])),
  cellID = as.integer(sapply(cell_parts, function(x) x[4])),
  stringsAsFactors = FALSE
)
cell_mapping$cell_key <- paste(cell_mapping$fov, cell_mapping$cellID, sep = "_")

# Find matching cells
matched_cells <- cell_mapping$cell_key %in% unique_poly_cells
cat("Cells matched to polygons:", sum(matched_cells), "\n")
cat("Cells without polygons:", sum(!matched_cells), "\n\n")

# =============================================================================
# GET EXPRESSION AND METADATA
# =============================================================================

expr_mat <- GetAssayData(seurat_sub, assay = "RNA", layer = "data")

cell_metadata <- data.frame(
  Cell = colnames(seurat_sub),
  Lineage = seurat_sub$Lineage,
  CellType = ifelse(seurat_sub$Lineage == "Cancer", "Cancer", "Non-Cancer"),
  stringsAsFactors = FALSE
)
cell_metadata <- merge(cell_metadata, cell_mapping, by = "Cell")

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

# Merge metadata with scores
plot_data <- merge(cell_metadata, probe_scores, by = "Cell")
cat("Cells with scores:", nrow(plot_data), "\n")

# Filter to cells with polygons
plot_data <- plot_data %>% filter(cell_key %in% unique_poly_cells)
cat("Cells with polygons:", nrow(plot_data), "\n\n")

# =============================================================================
# MERGE POLYGON DATA WITH SCORES
# =============================================================================

cat("Merging polygon coordinates with scores...\n")

# Get probe columns
probe_cols <- names(PROBES)
score_cols <- c("Cell", "cell_key", "Lineage", "CellType", probe_cols)

# Merge polygons with cell scores
poly_scores <- merge(polygons, plot_data[, score_cols], by = "cell_key")
cat("Total polygon vertices with scores:", nrow(poly_scores), "\n\n")

# =============================================================================
# CREATE SPATIAL POLYGON PLOTS
# =============================================================================

cat("Creating spatial polygon plots...\n\n")

# Function to create a single polygon plot
create_polygon_plot <- function(data, probe_name, probe_info, point_size = 0.1) {

  df <- data[!is.na(data[[probe_name]]), ]
  df$Score <- df[[probe_name]]

  # Cap extreme values
  df$Score[df$Score > 3] <- 3
  df$Score[df$Score < -3] <- -3

  # Color scale based on probe type
  if (probe_info$type == "GAIN") {
    colors <- c("#2166AC", "#F7F7F7", "#B2182B")  # Blue-White-Red
  } else {
    colors <- c("#2166AC", "#F7F7F7", "#B2182B")  # Same scale, interpretation differs
  }

  p <- ggplot(df, aes(x = x_global_px, y = y_global_px,
                       group = cell_key, fill = Score)) +
    geom_polygon(color = NA, linewidth = 0) +
    scale_fill_gradient2(
      low = colors[1], mid = colors[2], high = colors[3],
      midpoint = 0,
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
      title = sprintf("%s (%s)", probe_name, probe_info$type),
      subtitle = probe_info$description
    )

  return(p)
}

# Create individual plots for each probe
for (probe_name in names(PROBES)) {
  probe <- PROBES[[probe_name]]

  if (all(is.na(plot_data[[probe_name]]))) {
    cat(sprintf("  Skipping %s (no data)\n", probe_name))
    next
  }

  cat(sprintf("  Plotting %s...\n", probe_name))

  p <- create_polygon_plot(poly_scores, probe_name, probe)

  ggsave(
    file.path(OUTPUT_DIR, sprintf("PatientD_Bx3_Polygon_%s.pdf", probe_name)),
    p, width = 12, height = 10
  )
}

# =============================================================================
# CREATE CELL TYPE REFERENCE PLOT
# =============================================================================

cat("\nCreating cell type reference plot...\n")

p_celltype <- ggplot(poly_scores, aes(x = x_global_px, y = y_global_px,
                                       group = cell_key, fill = CellType)) +
  geom_polygon(color = NA, linewidth = 0) +
  scale_fill_manual(values = c("Cancer" = "#E41A1C", "Non-Cancer" = "#377EB8"),
                    name = "Cell Type") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  coord_fixed() +
  labs(title = "Patient D Bx3 - Cell Type (Polygons)")

ggsave(
  file.path(OUTPUT_DIR, "PatientD_Bx3_Polygon_CellType.pdf"),
  p_celltype, width = 12, height = 10
)

# =============================================================================
# CREATE TOP PROBES COMBINED FIGURE
# =============================================================================

cat("\nCreating combined top probes figure...\n")

# Best GAIN and LOSS probes
top_gain <- c("11q13_amp", "8q_gain", "17q22_23_amp")
top_loss <- c("22q12_CHEK2", "16q_loss", "8p_loss")

# Create mini plots without legends
create_mini_plot <- function(data, probe_name, probe_info) {
  df <- data[!is.na(data[[probe_name]]), ]
  df$Score <- df[[probe_name]]
  df$Score[df$Score > 3] <- 3
  df$Score[df$Score < -3] <- -3

  p <- ggplot(df, aes(x = x_global_px, y = y_global_px,
                       group = cell_key, fill = Score)) +
    geom_polygon(color = NA, linewidth = 0) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-3, 3),
      oob = scales::squish
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
      legend.position = "none"
    ) +
    coord_fixed() +
    labs(title = probe_name)

  return(p)
}

# Create gain plots
gain_plots <- lapply(top_gain, function(pn) {
  create_mini_plot(poly_scores, pn, PROBES[[pn]])
})
names(gain_plots) <- top_gain

# Create loss plots
loss_plots <- lapply(top_loss, function(pn) {
  create_mini_plot(poly_scores, pn, PROBES[[pn]])
})
names(loss_plots) <- top_loss

# Cell type mini plot
celltype_mini <- ggplot(poly_scores, aes(x = x_global_px, y = y_global_px,
                                          group = cell_key, fill = CellType)) +
  geom_polygon(color = NA, linewidth = 0) +
  scale_fill_manual(values = c("Cancer" = "#E41A1C", "Non-Cancer" = "#377EB8")) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
    legend.position = "none"
  ) +
  coord_fixed() +
  labs(title = "Cell Type Reference")

# Combine all plots
all_plots <- c(list(CellType = celltype_mini), gain_plots, loss_plots)

combined <- wrap_plots(all_plots, ncol = 4) +
  plot_annotation(
    title = "Patient D Bx3 - CNV Probe Scores (Cell Polygons)",
    subtitle = "Row 1: Cell Type + GAIN probes | Row 2: LOSS probes\nBlue = Lower Z-score, Red = Higher Z-score",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

ggsave(
  file.path(OUTPUT_DIR, "PatientD_Bx3_Polygon_TopProbes_Grid.pdf"),
  combined, width = 24, height = 14
)
cat("Saved: PatientD_Bx3_Polygon_TopProbes_Grid.pdf\n")

# =============================================================================
# CREATE SIDE-BY-SIDE: CNV PROBE + CELL TYPE
# =============================================================================

cat("\nCreating side-by-side comparisons...\n")

for (probe_name in c(top_gain, top_loss)) {
  probe <- PROBES[[probe_name]]

  # CNV score plot
  p1 <- create_polygon_plot(poly_scores, probe_name, probe)
  p1 <- p1 + theme(legend.position = "bottom") +
    labs(subtitle = NULL)

  # Cell type plot
  p2 <- ggplot(poly_scores, aes(x = x_global_px, y = y_global_px,
                                 group = cell_key, fill = CellType)) +
    geom_polygon(color = NA, linewidth = 0) +
    scale_fill_manual(values = c("Cancer" = "#E41A1C", "Non-Cancer" = "#377EB8"),
                      name = "Cell Type") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "bottom",
      plot.margin = margin(10, 10, 10, 10)
    ) +
    coord_fixed() +
    labs(title = "Cell Type")

  combined <- p1 + p2 +
    plot_annotation(
      title = sprintf("%s - %s", probe_name, probe$description),
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )

  ggsave(
    file.path(OUTPUT_DIR, sprintf("PatientD_Bx3_Polygon_%s_withCellType.pdf", probe_name)),
    combined, width = 20, height = 9
  )
}

# =============================================================================
# CREATE GAIN vs LOSS COMPARISON
# =============================================================================

cat("\nCreating GAIN vs LOSS comparison figure...\n")

gain_row <- wrap_plots(gain_plots, ncol = 3) +
  plot_annotation(
    title = "GAIN Probes (Red = Higher in Cancer)",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12))
  )

loss_row <- wrap_plots(loss_plots, ncol = 3) +
  plot_annotation(
    title = "LOSS Probes (Blue = Lower in Cancer)",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12))
  )

final_comparison <- gain_row / loss_row +
  plot_annotation(
    title = "Patient D Bx3 - Top CNV Probes (Cell Polygons)",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

ggsave(
  file.path(OUTPUT_DIR, "PatientD_Bx3_Polygon_GainVsLoss.pdf"),
  final_comparison, width = 22, height = 16
)
cat("Saved: PatientD_Bx3_Polygon_GainVsLoss.pdf\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=============================================================\n")
cat("OUTPUT FILES\n")
cat("=============================================================\n\n")

cat("Individual probe plots (with polygons):\n")
for (probe_name in names(PROBES)) {
  cat(sprintf("  - PatientD_Bx3_Polygon_%s.pdf\n", probe_name))
}

cat("\nReference:\n")
cat("  - PatientD_Bx3_Polygon_CellType.pdf\n")

cat("\nCombined figures:\n")
cat("  - PatientD_Bx3_Polygon_TopProbes_Grid.pdf\n")
cat("  - PatientD_Bx3_Polygon_GainVsLoss.pdf\n")

cat("\nSide-by-side (Probe + Cell Type):\n")
for (probe_name in c(top_gain, top_loss)) {
  cat(sprintf("  - PatientD_Bx3_Polygon_%s_withCellType.pdf\n", probe_name))
}

cat(sprintf("\nOutput directory: %s\n", OUTPUT_DIR))

cat("\n=============================================================\n")
cat("DONE\n")
cat("=============================================================\n")
