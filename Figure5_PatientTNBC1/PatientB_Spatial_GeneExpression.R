# ==============================================================================
# Patient B: Spatial Gene Expression Plots
# Cancer cells colored by expression (grey to orange to red)
# Non-cancer cells in light grey
# Same scale for each gene across Bx1 and Bx4
# ==============================================================================

library(Seurat)
library(tidyverse)
library(sf)
library(patchwork)

set.seed(42)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

output_dir <- "PatientB_Analysis/Spatial_Genes/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

POLYGON_DIR <- "supplementary_input_data/Polygons/RNA/"
SEURAT_PATH <- "data/CosMx_SMMART_345k_clean.rds"

# Polygon files for Patient B
POLYGON_FILES <- list(
  "Bx1" = "R1124_265607-polygons.csv",
  "Bx4" = "R1124_322118-polygons.csv"
)

# ==============================================================================
# GENES TO PLOT
# ==============================================================================

genes_to_plot <- c(
  # Basal identity (UP)
  "KRT17", "KRT5",
  # Antigen presentation (UP)
  "B2M", "HLA-A",
  # Stress adaptation
  "SQSTM1",
  # Down-regulated
  "KRT8", "NPPC",
  # Treatment relevant
  "TACSTD2"
)

cat("Genes to plot:", paste(genes_to_plot, collapse = ", "), "\n")

# ==============================================================================
# COLOR SCALE: Grey to Orange to Red
# ==============================================================================

expression_colors <- c("grey90", "#FED976", "#FD8D3C", "#E31A1C", "#800026")

# ==============================================================================
# COMMON THEME
# ==============================================================================

spatial_theme <- theme_void(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(face = "bold", size = 12)
  )

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("\n========== LOADING DATA ==========\n\n")

cat("Loading Seurat object...\n")
obj_all <- readRDS(SEURAT_PATH)
obj_all <- UpdateSeuratObject(obj_all)

cat("Total cells:", ncol(obj_all), "\n")

# Filter to Patient B
obj_patB <- subset(obj_all, subset = Patient == "Patient_B")
cat("Patient B cells:", ncol(obj_patB), "\n")

# Extract metadata
meta <- obj_patB@meta.data
meta$cell_id <- rownames(meta)

# Parse cell IDs
parts <- strsplit(meta$cell_id, "_")
meta$seurat_sample <- sapply(parts, function(x) as.integer(x[2]))
meta$seurat_fov <- sapply(parts, function(x) as.integer(x[3]))
meta$seurat_cellID <- sapply(parts, function(x) as.integer(x[4]))

# Get expression data
expr_mat <- GetAssayData(obj_patB, layer = "data")

# Check gene availability
genes_avail <- genes_to_plot[genes_to_plot %in% rownames(expr_mat)]
cat("Genes available:", paste(genes_avail, collapse = ", "), "\n")

# Add expression to metadata
for(gene in genes_avail) {
  meta[[gene]] <- as.numeric(expr_mat[gene, rownames(meta)])
}

# Check lineage column
cat("\nLineage distribution:\n")
print(table(meta$Lineage))

# Identify cancer cells
meta$is_cancer <- grepl("Cancer|cancer", meta$Lineage, ignore.case = TRUE)
cat("\nCancer cells:", sum(meta$is_cancer), "\n")

# ==============================================================================
# FUNCTION: Create SF polygons (GLOBAL coordinates)
# ==============================================================================

create_sf_global <- function(poly_df) {

  cell_polys <- poly_df %>%
    group_by(fov, cellID) %>%
    filter(n() >= 3) %>%
    summarise(
      geometry = list(tryCatch({
        coords <- cbind(x_global_px, y_global_px)
        if (nrow(coords) < 3) return(NULL)
        if (!all(coords[1,] == coords[nrow(coords),])) {
          coords <- rbind(coords, coords[1,])
        }
        st_polygon(list(coords))
      }, error = function(e) NULL)),
      .groups = "drop"
    ) %>%
    filter(!sapply(geometry, is.null))

  if (nrow(cell_polys) == 0) return(NULL)

  st_sf(fov = cell_polys$fov, cellID = cell_polys$cellID,
        geometry = st_sfc(cell_polys$geometry))
}

# ==============================================================================
# FUNCTION: Create Gene Expression Plot
# ==============================================================================

create_gene_plot <- function(sf_polys, biopsy_meta, gene, expr_limits, title_text) {

  # Create matching key
  sf_polys$match_key <- paste0(sf_polys$fov, "_", sf_polys$cellID)

  meta_lookup <- biopsy_meta %>%
    mutate(match_key = paste0(seurat_fov, "_", seurat_cellID)) %>%
    select(match_key, is_cancer, all_of(gene))

  sf_matched <- sf_polys %>%
    left_join(meta_lookup, by = "match_key")

  # Separate cancer and non-cancer
  sf_cancer <- sf_matched %>% filter(is_cancer == TRUE)
  sf_noncancer <- sf_matched %>% filter(is_cancer == FALSE | is.na(is_cancer))

  # Plot
  p <- ggplot() +
    # Non-cancer cells in light grey
    geom_sf(data = sf_noncancer, fill = "grey85", color = "grey70", linewidth = 0.01) +
    # Cancer cells colored by expression
    geom_sf(data = sf_cancer, aes(fill = .data[[gene]]), color = "grey30", linewidth = 0.01) +
    scale_fill_gradientn(
      colors = expression_colors,
      limits = expr_limits,
      name = "Expression",
      na.value = "grey85",
      oob = scales::squish
    ) +
    labs(title = title_text) +
    spatial_theme

  return(p)
}

# ==============================================================================
# CALCULATE GLOBAL EXPRESSION LIMITS PER GENE
# ==============================================================================

cat("\n========== CALCULATING EXPRESSION RANGES ==========\n\n")

# Calculate limits across both biopsies for each gene
expr_limits_list <- list()

for(gene in genes_avail) {
  cancer_expr <- meta[[gene]][meta$is_cancer]
  # Use 0 to 99th percentile to avoid outliers
  expr_limits_list[[gene]] <- c(0, quantile(cancer_expr, 0.99, na.rm = TRUE))
  cat(sprintf("%s: 0 to %.2f (99th percentile)\n", gene, expr_limits_list[[gene]][2]))
}

# ==============================================================================
# PROCESS EACH BIOPSY
# ==============================================================================

cat("\n========== CREATING SPATIAL PLOTS ==========\n\n")

# Store plots for combined figure
all_plots <- list()

for(biopsy in c("Bx1", "Bx4")) {

  cat("\n--- Processing", biopsy, "---\n")

  # Load polygon file
  poly_path <- file.path(POLYGON_DIR, POLYGON_FILES[[biopsy]])

  if(!file.exists(poly_path)) {
    cat("WARNING: Polygon file not found:", poly_path, "\n")
    next
  }

  poly_df <- read.csv(poly_path)
  cat("Polygon vertices:", nrow(poly_df), "\n")

  # Get metadata for this biopsy
  biopsy_meta <- meta %>% filter(Timepoint == biopsy)
  cat("Cells in biopsy:", nrow(biopsy_meta), "\n")
  cat("Cancer cells:", sum(biopsy_meta$is_cancer), "\n")

  # Create SF polygons
  sf_polys <- create_sf_global(poly_df)

  if(is.null(sf_polys)) {
    cat("WARNING: Could not create polygons\n")
    next
  }

  cat("Polygons created:", nrow(sf_polys), "\n")

  # Create plot for each gene
  for(gene in genes_avail) {

    site_label <- ifelse(biopsy == "Bx1", "Lymph Node", "Soft Tissue")
    title_text <- paste0(gene, " - ", biopsy, " (", site_label, ")")

    p <- create_gene_plot(
      sf_polys = sf_polys,
      biopsy_meta = biopsy_meta,
      gene = gene,
      expr_limits = expr_limits_list[[gene]],
      title_text = title_text
    )

    # Save individual plot
    ggsave(
      file.path(output_dir, paste0("PatientB_", biopsy, "_", gene, ".pdf")),
      p, width = 8, height = 8, bg = "white"
    )

    # Store for combined
    all_plots[[paste0(gene, "_", biopsy)]] <- p

    cat("  Saved:", gene, "\n")
  }

  rm(poly_df, sf_polys)
  gc()
}

# ==============================================================================
# CREATE COMBINED FIGURES (Bx1 vs Bx4 side by side)
# ==============================================================================

cat("\n========== CREATING COMBINED FIGURES ==========\n\n")

for(gene in genes_avail) {

  p_bx1 <- all_plots[[paste0(gene, "_Bx1")]]
  p_bx4 <- all_plots[[paste0(gene, "_Bx4")]]

  if(!is.null(p_bx1) && !is.null(p_bx4)) {

    # Get DE stats for subtitle
    de_results <- read.csv(file.path("LR_Analysis/PatientB_Analysis/",
                                      "PatientB_DE_results.csv"))
    gene_stats <- de_results %>% filter(Gene == gene)

    if(nrow(gene_stats) > 0) {
      subtitle_text <- sprintf("Log2FC: %+.2f | Cancer cells only (non-cancer in grey)",
                               gene_stats$avg_log2FC[1])
    } else {
      subtitle_text <- "Cancer cells only (non-cancer in grey)"
    }

    combined <- p_bx1 + p_bx4 +
      plot_layout(ncol = 2, guides = "collect") +
      plot_annotation(
        title = paste0("Patient B: ", gene, " Expression"),
        subtitle = subtitle_text,
        theme = theme(
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          plot.subtitle = element_text(color = "grey40", size = 11, hjust = 0.5)
        )
      )

    ggsave(
      file.path(output_dir, paste0("PatientB_Combined_", gene, ".pdf")),
      combined, width = 16, height = 8, bg = "white"
    )

    cat("Combined:", gene, "\n")
  }
}

# ==============================================================================
# CREATE MULTI-GENE PANEL
# ==============================================================================

cat("\n========== CREATING MULTI-GENE PANEL ==========\n\n")

# Select key genes for panel
panel_genes <- c("KRT17", "KRT5", "B2M", "KRT8")
panel_genes <- panel_genes[panel_genes %in% genes_avail]

if(length(panel_genes) >= 2) {

  # Create 2-row layout: Bx1 on top, Bx4 on bottom
  panel_plots <- list()

  for(gene in panel_genes) {
    panel_plots[[paste0(gene, "_Bx1")]] <- all_plots[[paste0(gene, "_Bx1")]] +
      labs(title = gene) +
      theme(legend.position = "none")
    panel_plots[[paste0(gene, "_Bx4")]] <- all_plots[[paste0(gene, "_Bx4")]] +
      labs(title = NULL) +
      theme(legend.position = "none")
  }

  # Arrange: genes as columns, biopsies as rows
  n_genes <- length(panel_genes)

  # Row 1: Bx1
  row1 <- wrap_plots(panel_plots[paste0(panel_genes, "_Bx1")], ncol = n_genes)
  # Row 2: Bx4
  row2 <- wrap_plots(panel_plots[paste0(panel_genes, "_Bx4")], ncol = n_genes)

  multi_panel <- row1 / row2 +
    plot_annotation(
      title = "Patient B: Key Gene Expression Changes",
      subtitle = "Top: Bx1 (Lymph Node) | Bottom: Bx4 (Soft Tissue) | Cancer cells only",
      caption = "Same color scale used for each gene across both biopsies",
      theme = theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(color = "grey40", size = 11, hjust = 0.5),
        plot.caption = element_text(color = "grey50", size = 9, hjust = 0.5)
      )
    )

  ggsave(
    file.path(output_dir, "PatientB_MultiGene_Panel.pdf"),
    multi_panel, width = 4 * n_genes, height = 10, bg = "white"
  )

  cat("Multi-gene panel saved\n")
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n========== COMPLETE ==========\n\n")
cat("Output directory:", output_dir, "\n\n")
cat("Files created:\n")
cat("  Individual plots: PatientB_[Bx]_[Gene].pdf\n")
cat("  Combined (side-by-side): PatientB_Combined_[Gene].pdf\n")
cat("  Multi-gene panel: PatientB_MultiGene_Panel.pdf\n")
cat("\nColor scale: grey90 → yellow → orange → red → dark red\n")
cat("Non-cancer cells shown in grey85\n")
cat("Same scale used for each gene across both biopsies\n")
