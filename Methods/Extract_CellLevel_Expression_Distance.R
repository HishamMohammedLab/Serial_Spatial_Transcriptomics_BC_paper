#!/usr/bin/env Rscript
# =============================================================================
# EXTRACT CELL-LEVEL EXPRESSION LINKED TO DISTANCE DATA (v2 — All Pairs)
# =============================================================================
#
# PURPOSE: Create a dataset that links cell-pair distances with the actual
#          expression levels of ligands and receptors in those specific cells.
#
# v2 CHANGES:
#   - Points to AllPairs raw distance file (462 pairs, 3.2M rows)
#   - Handles simpler column format (adds _um conversion, direct_contact, etc.)
#   - Optional PAIR_FILTER: set to a CSV path to restrict to specific pairs
#   - Outputs to LR_Distance_Expression_AllPairs.csv
#
# OPTIMIZED: Uses data.table joins for maximum speed
#
# =============================================================================

library(tidyverse)
library(Seurat)
library(data.table)

# =============================================================================
# CONFIGURATION
# =============================================================================

COMBINED_SEURAT <- "data/CosMx_SMMART_345k_clean.rds"
DISTANCE_FILE <- "supplementary_input_data/LR_Scores/LR_Distances_Raw_AllPairs.csv"
OUTPUT_DIR <- "Spatial_Distance_Complete"
OUTPUT_FILE <- file.path(OUTPUT_DIR, "LR_Distance_Expression_AllPairs.csv")

# Optional: set to a CSV file path to filter to specific pairs only (for test runs)
# The CSV should have an "LR_Pair" column with pairs in "GENE1->GENE2" format
PAIR_FILTER <- NULL  # e.g., "/path/to/curated_pairs.csv"

# Pixel to micron conversion factor (CosMx: 0.18 µm/px = 5.5556 px/µm)
PX_TO_UM <- 0.18

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading distance data...\n")
distance_data <- fread(DISTANCE_FILE)
cat(sprintf("  Loaded %d cell pairs, %d unique LR pairs\n",
            nrow(distance_data), uniqueN(distance_data$LR_Pair)))

# Apply pair filter if specified
if (!is.null(PAIR_FILTER) && file.exists(PAIR_FILTER)) {
  filter_pairs <- fread(PAIR_FILTER)$LR_Pair
  cat(sprintf("  Filtering to %d pairs from %s\n", length(filter_pairs), PAIR_FILTER))
  distance_data <- distance_data[LR_Pair %in% filter_pairs]
  cat(sprintf("  After filtering: %d cell pairs, %d unique LR pairs\n",
              nrow(distance_data), uniqueN(distance_data$LR_Pair)))
}

# Handle AllPairs format: convert column names and add derived columns
if ("polygon_distance" %in% names(distance_data) && !"polygon_distance_um" %in% names(distance_data)) {
  cat("  Detected AllPairs format — adding unit conversions and derived columns...\n")

  # Rename raw distance columns to _px suffix
  setnames(distance_data, "centroid_distance", "centroid_distance_px", skip_absent = TRUE)
  setnames(distance_data, "polygon_distance", "polygon_distance_px", skip_absent = TRUE)

  # Convert to microns
  distance_data[, centroid_distance_um := centroid_distance_px * PX_TO_UM]
  distance_data[, polygon_distance_um := polygon_distance_px * PX_TO_UM]

  # Add direct_contact flag (polygon distance = 0)
  distance_data[, direct_contact := polygon_distance_px == 0]

  # Add distance_category
  distance_data[, distance_category := fcase(
    direct_contact == TRUE, "Direct_Contact",
    polygon_distance_um < 15, "Juxtacrine (<15\u00b5m)",
    polygon_distance_um < 30, "Close_Paracrine (15-30\u00b5m)",
    polygon_distance_um < 50, "Mid_Paracrine (30-50\u00b5m)",
    default = "Distal_Paracrine (>50\u00b5m)"
  )]

  # Rename Source_Type/Target_Type to source_lineage/target_lineage for consistency
  # (keep both — Source_Type as-is, add lineage columns)
  distance_data[, source_lineage := Source_Type]
  distance_data[, target_lineage := Target_Type]

  # Normalize arrow notation: -> to →
  distance_data[, LR_Pair := gsub("->", "\u2192", LR_Pair)]

  # Drop helper columns not needed
  distance_data[, c("source_match_key", "target_match_key") := NULL]

  cat(sprintf("  Converted distances: median %.1f um, range %.1f-%.1f um\n",
              median(distance_data$polygon_distance_um),
              min(distance_data$polygon_distance_um),
              max(distance_data$polygon_distance_um)))
}

# Get unique genes needed
unique_ligands <- unique(distance_data$Ligand)
unique_receptors <- unique(distance_data$Receptor)
genes_needed <- unique(c(unique_ligands, unique_receptors))
cat(sprintf("  Need expression for %d genes\n", length(genes_needed)))

# =============================================================================
# LOAD SEURAT AND EXTRACT EXPRESSION
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("LOADING SEURAT OBJECT\n")
cat(strrep("=", 60), "\n")

cat("Loading combined Seurat object...\n")
obj <- readRDS(COMBINED_SEURAT)
cat(sprintf("  Loaded %d cells, %d genes\n", ncol(obj), nrow(obj)))

# Get expression matrix
cat("Extracting expression matrix...\n")
expr_matrix <- tryCatch({
  GetAssayData(obj, layer = "data")
}, error = function(e) {
  GetAssayData(obj, slot = "data")
})

# Check genes
available_genes <- intersect(genes_needed, rownames(expr_matrix))
missing_genes <- setdiff(genes_needed, rownames(expr_matrix))
cat(sprintf("  Genes found: %d / %d\n", length(available_genes), length(genes_needed)))

# =============================================================================
# CREATE LONG-FORMAT EXPRESSION TABLE
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("CREATING EXPRESSION LOOKUP TABLE\n")
cat(strrep("=", 60), "\n")

# Get expression for genes we need
cat("Extracting relevant genes from expression matrix...\n")
expr_subset <- expr_matrix[available_genes, , drop=FALSE]

# Convert to long format data.table for efficient joining
cat("Converting to long format (this may take a moment)...\n")
expr_long <- as.data.table(as.matrix(expr_subset), keep.rownames = "Gene")
expr_long <- melt(expr_long, id.vars = "Gene", variable.name = "Cell", value.name = "Expression")
expr_long[, Cell := as.character(Cell)]

cat(sprintf("  Expression table: %d rows\n", nrow(expr_long)))

# Set keys for fast joining
setkey(expr_long, Cell, Gene)

# =============================================================================
# JOIN EXPRESSION TO DISTANCE DATA
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("JOINING EXPRESSION TO DISTANCE DATA\n")
cat(strrep("=", 60), "\n")

# Create lookup tables for ligand and receptor separately
ligand_expr <- copy(expr_long)
setnames(ligand_expr, c("Gene", "Cell", "Expression"), c("Ligand", "source_cell", "Ligand_Expr"))
setkey(ligand_expr, source_cell, Ligand)

receptor_expr <- copy(expr_long)
setnames(receptor_expr, c("Gene", "Cell", "Expression"), c("Receptor", "target_cell", "Receptor_Expr"))
setkey(receptor_expr, target_cell, Receptor)

# Set keys on distance data
setkey(distance_data, source_cell, Ligand)

# Join ligand expression
cat("Joining ligand expression...\n")
distance_data <- ligand_expr[distance_data, on = .(source_cell, Ligand)]

# Reset key and join receptor expression
setkey(distance_data, target_cell, Receptor)
cat("Joining receptor expression...\n")
distance_data <- receptor_expr[distance_data, on = .(target_cell, Receptor)]

# Calculate LR product
cat("Calculating LR product...\n")
distance_data[, LR_Product := Ligand_Expr * Receptor_Expr]

# =============================================================================
# SUMMARY AND SAVE
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("SUMMARY\n")
cat(strrep("=", 60), "\n")

n_complete <- sum(!is.na(distance_data$Ligand_Expr) & !is.na(distance_data$Receptor_Expr))
cat(sprintf("Total cell pairs: %d\n", nrow(distance_data)))
cat(sprintf("Pairs with ligand expression: %d\n", sum(!is.na(distance_data$Ligand_Expr))))
cat(sprintf("Pairs with receptor expression: %d\n", sum(!is.na(distance_data$Receptor_Expr))))
cat(sprintf("Pairs with complete expression: %d (%.1f%%)\n",
            n_complete, 100 * n_complete / nrow(distance_data)))

cat(sprintf("\nMean Ligand Expression: %.3f\n", mean(distance_data$Ligand_Expr, na.rm=TRUE)))
cat(sprintf("Mean Receptor Expression: %.3f\n", mean(distance_data$Receptor_Expr, na.rm=TRUE)))
cat(sprintf("Mean LR Product: %.3f\n", mean(distance_data$LR_Product, na.rm=TRUE)))

# Save main output
cat(sprintf("\nSaving to: %s\n", OUTPUT_FILE))
fwrite(distance_data, OUTPUT_FILE)

# Create summary by distance bin
cat("Creating distance-binned summary...\n")
distance_data[, Distance_Bin := cut(polygon_distance_um,
                                     breaks = c(-Inf, 0, 5, 10, 15, 20, 30, 50, 100, Inf),
                                     labels = c("Contact", "0-5", "5-10", "10-15", "15-20", "20-30", "30-50", "50-100", "100+"))]

summary_by_distance <- distance_data[!is.na(Ligand_Expr) & !is.na(Receptor_Expr),
                                      .(N_Pairs = .N,
                                        Mean_Ligand_Expr = mean(Ligand_Expr, na.rm = TRUE),
                                        Mean_Receptor_Expr = mean(Receptor_Expr, na.rm = TRUE),
                                        Mean_LR_Product = mean(LR_Product, na.rm = TRUE),
                                        SD_LR_Product = sd(LR_Product, na.rm = TRUE)),
                                      by = .(LR_Pair, Distance_Bin)]

summary_file <- file.path(OUTPUT_DIR, "LR_Distance_Expression_AllPairs_Summary.csv")
fwrite(summary_by_distance, summary_file)
cat(sprintf("Saved summary to: %s\n", summary_file))

# Also save a summary by LR_Pair + Source_Type + Target_Type + Distance_Bin
cat("Creating context-level distance-binned summary...\n")
summary_by_context <- distance_data[!is.na(Ligand_Expr) & !is.na(Receptor_Expr),
                                     .(N_Pairs = .N,
                                       Mean_Ligand_Expr = mean(Ligand_Expr, na.rm = TRUE),
                                       Mean_Receptor_Expr = mean(Receptor_Expr, na.rm = TRUE),
                                       Mean_LR_Product = mean(LR_Product, na.rm = TRUE),
                                       SD_LR_Product = sd(LR_Product, na.rm = TRUE)),
                                     by = .(LR_Pair, Source_Type, Target_Type, Distance_Bin)]
context_summary_file <- file.path(OUTPUT_DIR, "LR_Distance_Expression_AllPairs_ContextSummary.csv")
fwrite(summary_by_context, context_summary_file)
cat(sprintf("Saved context summary to: %s\n", context_summary_file))

# Quick preview of the relationship
cat("\n", strrep("=", 60), "\n")
cat("PREVIEW: EXPRESSION BY DISTANCE\n")
cat(strrep("=", 60), "\n")

overall_summary <- distance_data[!is.na(LR_Product),
                                  .(N = .N,
                                    Mean_LR_Product = mean(LR_Product)),
                                  by = Distance_Bin][order(Distance_Bin)]
print(overall_summary)

cat("\n", strrep("=", 60), "\n")
cat("DONE\n")
cat(strrep("=", 60), "\n")
