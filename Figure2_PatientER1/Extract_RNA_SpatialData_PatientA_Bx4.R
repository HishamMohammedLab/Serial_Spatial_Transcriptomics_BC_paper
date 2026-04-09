# Extract RNA expression data for spatial plotting - Patient A Bx4

library(Seurat)
library(tidyverse)

# =============================================================================
# CONFIGURATION
# =============================================================================

SEURAT_PATH <- "data/CosMx_SMMART_345k_clean.rds"
POLYGON_PATH <- "supplementary_input_data/Polygons/RNA/R1134_322078-polygons.csv"
OUTPUT_DIR <- "Spatial_Expression_Plots"

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Target genes
TARGET_GENES <- c("ESR1", "TGFBR2", "JAG1", "SNAI1", "MKI67", "MX1", "STAT1",
                  "HLA-A", "HLA-B", "SPP1", "WNT5A", "HGF", "MET", "MIF",
                  "CD74", "CDH1")

# =============================================================================
# LOAD DATA
# =============================================================================
cat("Loading Seurat object...\n")
seurat_obj <- readRDS(SEURAT_PATH)

cat("Total cells:", ncol(seurat_obj), "\n")

# Filter to Patient A Bx4
cat("Filtering to Patient A Bx4...\n")
cells_patA_bx4 <- colnames(seurat_obj)[seurat_obj$Patient == "Patient_A" &
                                        seurat_obj$Timepoint == "Bx4"]
cat("Patient A Bx4 cells:", length(cells_patA_bx4), "\n")

# Subset
seurat_bx4 <- subset(seurat_obj, cells = cells_patA_bx4)

# =============================================================================
# EXTRACT EXPRESSION DATA
# =============================================================================
cat("Extracting expression data...\n")

# Get expression matrix (use layer for Seurat v5)
expr_matrix <- GetAssayData(seurat_bx4, layer = "data")

# Check which genes are available
available_genes <- TARGET_GENES[TARGET_GENES %in% rownames(expr_matrix)]
missing_genes <- TARGET_GENES[!TARGET_GENES %in% rownames(expr_matrix)]

cat("\nAvailable genes:", paste(available_genes, collapse = ", "), "\n")
if (length(missing_genes) > 0) {
  cat("Missing genes:", paste(missing_genes, collapse = ", "), "\n")
}

# Extract expression for available genes
expr_df <- as.data.frame(t(as.matrix(expr_matrix[available_genes, , drop = FALSE])))
expr_df$cell_id <- rownames(expr_df)

# Add metadata
meta <- seurat_bx4@meta.data
meta$cell_id <- rownames(meta)

# Merge
expr_df <- merge(expr_df, meta[, c("cell_id", "x_global_px", "y_global_px", "fov",
                                    "Lineage", "predicted.id")],
                 by = "cell_id")

cat("\nExpression data shape:", nrow(expr_df), "cells x", ncol(expr_df), "columns\n")

# =============================================================================
# LOAD POLYGONS
# =============================================================================
cat("\nLoading polygon data...\n")
polygons <- read.csv(POLYGON_PATH)
cat("Polygon rows:", nrow(polygons), "\n")

# Create cell_id matching key
polygons$match_key <- paste0(polygons$fov, "_", polygons$cellID)

# =============================================================================
# SAVE DATA
# =============================================================================

# Save expression data
expr_output <- file.path("supplementary_input_data/PatientA_Bx4_RNA_Expression.csv")
write.csv(expr_df, expr_output, row.names = FALSE)
cat("\nSaved expression data to:", expr_output, "\n")

# Save polygon data (subset to unique cells for centroid plotting)
centroids <- polygons %>%
  group_by(fov, cellID) %>%
  summarize(
    x_centroid = mean(x_global_px),
    y_centroid = mean(y_global_px),
    .groups = "drop"
  )

centroid_output <- file.path("supplementary_input_data/PatientA_Bx4_RNA_Centroids.csv")
write.csv(centroids, centroid_output, row.names = FALSE)
cat("Saved centroids to:", centroid_output, "\n")

# Save full polygon data for polygon plotting
polygon_output <- file.path("supplementary_input_data/PatientA_Bx4_RNA_Polygons.csv")
write.csv(polygons, polygon_output, row.names = FALSE)
cat("Saved full polygons to:", polygon_output, "\n")

cat("\n=== DONE ===\n")
cat("Available genes for plotting:\n")
for (g in available_genes) {
  expr_vals <- expr_df[[g]]
  cat(sprintf("  %s: %.1f%% expressing, mean=%.3f\n",
              g, 100*mean(expr_vals > 0), mean(expr_vals)))
}

