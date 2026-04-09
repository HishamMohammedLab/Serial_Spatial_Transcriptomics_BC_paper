#!/usr/bin/env Rscript
# =============================================================================
# PROMPT 1: SPATIAL GRADIENT ANALYSIS - Patient D Bx3
# Purpose: Characterize the spatial gradient from ER+ Cancer Core to ER-negative
#          invasive Tumor Nests in Patient D's final liver biopsy
# Output:  6 figure panels + CSV tables for downstream use
# =============================================================================

library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(msigdbr)
library(cowplot)
library(patchwork)
library(ggrepel)
library(viridis)

cat("=", rep("=", 69), "\n", sep = "")
cat("PROMPT 1: SPATIAL GRADIENT ANALYSIS - Patient D Bx3\n")
cat("=", rep("=", 69), "\n\n", sep = "")

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------

SEURAT_PATH <- "data/CosMx_SMMART_345k_clean.rds"
OUTPUT_DIR <- "Final_Figure/spatial_gradient_results/"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Domain mapping (verified in prompt0)
DOMAIN_TO_GROUP <- c(
  "0" = "Cancer_Core",
  "1" = "Diploid_Epithelial",
  "2" = "Cancer_Core",
  "3" = "Stromal_Interface",
  "4" = "Diploid_Epithelial",
  "5" = "Tumor_Nest",
  "6" = "Diploid_Epithelial",
  "7" = "Diploid_Epithelial",
  "8" = "Diploid_Epithelial",
  "9" = "Cancer_Core"
)

# Zone palette: diverging blue-to-red (core to nest)
ZONE_COLORS <- c(
  "Deep Core"          = "#2166AC",
  "Core Edge"          = "#67A9CF",
  "Stromal Interface"  = "#F7F7F7",
  "Proximal Nests"     = "#EF8A62",
  "Distal Nests"       = "#B2182B"
)

# Pixel to micron conversion
PX_TO_UM <- 0.18

# =============================================================================
# STEP 1: DEFINE SPATIAL ZONES
# =============================================================================
cat("\n--- STEP 1: Define Spatial Zones ---\n")

# Load Seurat object
cat("Loading Seurat object...\n")
seu <- readRDS(SEURAT_PATH)
cat(sprintf("  Total cells: %d\n", ncol(seu)))

# Subset to Patient D Bx3
d_bx3 <- subset(seu, Patient == "Patient_D" & Timepoint == "Bx3")
cat(sprintf("  Patient D Bx3 cells: %d\n", ncol(d_bx3)))

# Assign domain groups (unname to avoid Seurat metadata name mismatch)
d_bx3$Domain_Group <- unname(DOMAIN_TO_GROUP[as.character(d_bx3$spatial_domain)])
cat("  Domain group counts (all cell types):\n")
print(table(d_bx3$Domain_Group))

# Subset to CANCER cells only
d_bx3_cancer <- subset(d_bx3, Lineage == "Cancer")
cat(sprintf("\n  Cancer cells in D Bx3: %d\n", ncol(d_bx3_cancer)))
cat("  Cancer cells per domain group:\n")
print(table(d_bx3_cancer$Domain_Group))

# Get spatial coordinates in microns
meta_cancer <- d_bx3_cancer@meta.data
meta_cancer$x_um <- meta_cancer$x_global_px * PX_TO_UM
meta_cancer$y_um <- meta_cancer$y_global_px * PX_TO_UM

# Compute Cancer Core centroid
core_cells <- meta_cancer[meta_cancer$Domain_Group == "Cancer_Core", ]
centroid_x <- mean(core_cells$x_um, na.rm = TRUE)
centroid_y <- mean(core_cells$y_um, na.rm = TRUE)
cat(sprintf("\n  Cancer Core centroid: (%.1f, %.1f) um\n", centroid_x, centroid_y))

# Compute distance to centroid for every cancer cell
meta_cancer$dist_to_core <- sqrt(
  (meta_cancer$x_um - centroid_x)^2 +
  (meta_cancer$y_um - centroid_y)^2
)

# Assign 5 zones
assign_zone <- function(domain_group, dist, core_median, nest_median) {
  case_when(
    domain_group == "Cancer_Core" & dist <= core_median   ~ "Deep Core",
    domain_group == "Cancer_Core" & dist > core_median    ~ "Core Edge",
    domain_group == "Stromal_Interface"                   ~ "Stromal Interface",
    domain_group == "Tumor_Nest" & dist <= nest_median    ~ "Proximal Nests",
    domain_group == "Tumor_Nest" & dist > nest_median     ~ "Distal Nests",
    domain_group == "Diploid_Epithelial"                  ~ NA_character_,
    TRUE                                                  ~ NA_character_
  )
}

# Compute medians for splitting
core_dist_median <- median(
  meta_cancer$dist_to_core[meta_cancer$Domain_Group == "Cancer_Core"],
  na.rm = TRUE
)
nest_cells_mask <- meta_cancer$Domain_Group == "Tumor_Nest"
nest_dist_median <- if (sum(nest_cells_mask) > 0) {
  median(meta_cancer$dist_to_core[nest_cells_mask], na.rm = TRUE)
} else {
  Inf
}

cat(sprintf("  Core distance median: %.1f um\n", core_dist_median))
cat(sprintf("  Nest distance median: %.1f um\n", nest_dist_median))

meta_cancer$Zone <- assign_zone(
  meta_cancer$Domain_Group,
  meta_cancer$dist_to_core,
  core_dist_median,
  nest_dist_median
)

# Handle Stromal Interface: merge into Core Edge if <30 cells
si_count <- sum(meta_cancer$Zone == "Stromal Interface", na.rm = TRUE)
cat(sprintf("\n  Stromal Interface cancer cells: %d\n", si_count))
if (si_count < 30) {
  cat("  WARNING: <30 cells in Stromal Interface - merging into Core Edge\n")
  meta_cancer$Zone[meta_cancer$Zone == "Stromal Interface"] <- "Core Edge"
}

# Remove cells not in any zone (Diploid_Epithelial)
meta_cancer <- meta_cancer[!is.na(meta_cancer$Zone), ]

# Set zone order
zone_levels <- c("Deep Core", "Core Edge", "Stromal Interface", "Proximal Nests", "Distal Nests")
# Only keep levels that exist
zone_levels <- zone_levels[zone_levels %in% unique(meta_cancer$Zone)]
meta_cancer$Zone <- factor(meta_cancer$Zone, levels = zone_levels)

# Zone index for correlation analysis
meta_cancer$Zone_Index <- as.integer(meta_cancer$Zone)

# Print cell counts
cat("\n  Zone cell counts:\n")
zone_counts <- table(meta_cancer$Zone)
print(zone_counts)

# Flag zones with <50 cells
small_zones <- names(zone_counts[zone_counts < 50])
if (length(small_zones) > 0) {
  cat("  WARNING: Zones with <50 cells:", paste(small_zones, collapse = ", "), "\n")
}

# QC: mean nCount_RNA per zone
cat("\n  QC - Mean nCount_RNA per zone:\n")
qc_by_zone <- meta_cancer %>%
  group_by(Zone) %>%
  summarise(
    n_cells = n(),
    mean_nCount = mean(nCount_RNA, na.rm = TRUE),
    median_nCount = median(nCount_RNA, na.rm = TRUE),
    .groups = "drop"
  )
print(as.data.frame(qc_by_zone))

# Check for >2-fold depth bias
z1_depth <- qc_by_zone$mean_nCount[qc_by_zone$Zone == zone_levels[1]]
z5_depth <- qc_by_zone$mean_nCount[qc_by_zone$Zone == tail(zone_levels, 1)]
depth_ratio <- max(z1_depth, z5_depth) / min(z1_depth, z5_depth)
cat(sprintf("  Depth ratio (Zone1 vs Zone5): %.2f\n", depth_ratio))
if (depth_ratio > 2) {
  cat("  WARNING: >2-fold transcript depth difference between zones!\n")
}

# Save zone assignments
zone_df <- meta_cancer %>%
  select(Zone, Zone_Index, Domain_Group, dist_to_core, x_um, y_um,
         nCount_RNA, nFeature_RNA) %>%
  mutate(cell_id = rownames(meta_cancer))
write.csv(zone_df, file.path(OUTPUT_DIR, "zone_assignments.csv"), row.names = FALSE)
cat("  Saved: zone_assignments.csv\n")

# Add zone to Seurat metadata for the cancer subset
cells_with_zones <- rownames(meta_cancer)

# Free memory - keep only what we need
rm(d_bx3)
gc()

# =============================================================================
# STEP 2: GENE EXPRESSION PROFILES PER ZONE
# =============================================================================
cat("\n--- STEP 2: Gene Expression Profiles Per Zone ---\n")

# Get normalized expression for cancer cells with zones
expr_data <- GetAssayData(seu, layer = "data")[, cells_with_zones]
cat(sprintf("  Expression matrix: %d genes x %d cells\n", nrow(expr_data), ncol(expr_data)))

# Compute mean expression per gene per zone
cat("  Computing mean expression per zone...\n")
zone_labels <- meta_cancer$Zone
names(zone_labels) <- rownames(meta_cancer)

# Build genes x zones matrix
gene_by_zone <- matrix(NA, nrow = nrow(expr_data), ncol = length(zone_levels),
                       dimnames = list(rownames(expr_data), zone_levels))

for (z in zone_levels) {
  z_cells <- names(zone_labels[zone_labels == z])
  if (length(z_cells) > 0) {
    gene_by_zone[, z] <- rowMeans(expr_data[, z_cells, drop = FALSE])
  }
}

cat(sprintf("  Gene-by-zone matrix: %d x %d\n", nrow(gene_by_zone), ncol(gene_by_zone)))

# Compute mean topic scores per zone
topic_cols <- grep("^Topic_", colnames(meta_cancer), value = TRUE)
cat(sprintf("  Topic columns found: %d\n", length(topic_cols)))

topic_by_zone <- matrix(NA, nrow = length(topic_cols), ncol = length(zone_levels),
                        dimnames = list(topic_cols, zone_levels))

for (z in zone_levels) {
  z_idx <- meta_cancer$Zone == z
  for (tc in topic_cols) {
    topic_by_zone[tc, z] <- mean(meta_cancer[[tc]][z_idx], na.rm = TRUE)
  }
}

# Save both
write.csv(gene_by_zone, file.path(OUTPUT_DIR, "gene_expression_by_zone.csv"))
write.csv(topic_by_zone, file.path(OUTPUT_DIR, "topic_scores_by_zone.csv"))
saveRDS(list(gene_by_zone = gene_by_zone, topic_by_zone = topic_by_zone),
        file.path(OUTPUT_DIR, "expression_profiles.rds"))
cat("  Saved: gene_expression_by_zone.csv, topic_scores_by_zone.csv, expression_profiles.rds\n")

# =============================================================================
# STEP 3: IDENTIFY SPATIALLY GRADIENT GENES
# =============================================================================
cat("\n--- STEP 3: Identify Spatially Gradient Genes ---\n")

# For each gene: Spearman correlation of zone-level mean expression vs zone index
zone_indices <- 1:length(zone_levels)
n_genes <- nrow(gene_by_zone)

gradient_results <- data.frame(
  gene = rownames(gene_by_zone),
  rho = NA_real_,
  pvalue = NA_real_,
  stringsAsFactors = FALSE
)

for (i in 1:n_genes) {
  expr_vals <- gene_by_zone[i, ]
  # Only compute if there is variance
  if (sd(expr_vals, na.rm = TRUE) > 0) {
    ct <- cor.test(zone_indices, expr_vals, method = "spearman", exact = FALSE)
    gradient_results$rho[i] <- ct$estimate
    gradient_results$pvalue[i] <- ct$p.value
  } else {
    gradient_results$rho[i] <- 0
    gradient_results$pvalue[i] <- 1
  }
}

# Assign direction
gradient_results <- gradient_results %>%
  mutate(
    direction = case_when(
      rho < 0 ~ "core-enriched",
      rho > 0 ~ "nest-enriched",
      TRUE    ~ "flat"
    ),
    category = case_when(
      abs(rho) > 0.8   ~ "strongly_monotonic",
      abs(rho) > 0.5   ~ "moderate",
      TRUE             ~ "weak_or_nonmonotonic"
    ),
    # Add expression info
    mean_expr = rowMeans(gene_by_zone),
    max_zone_expr = apply(gene_by_zone, 1, max)
  ) %>%
  arrange(rho)

# Print top 20 core-enriched, top 20 nest-enriched
cat("\n  Top 20 CORE-enriched genes (rho < 0):\n")
core_top <- gradient_results %>% filter(rho < 0) %>% head(20)
for (i in 1:nrow(core_top)) {
  cat(sprintf("    %s: rho=%.3f\n", core_top$gene[i], core_top$rho[i]))
}

cat("\n  Top 20 NEST-enriched genes (rho > 0):\n")
nest_top <- gradient_results %>% filter(rho > 0) %>% arrange(desc(rho)) %>% head(20)
for (i in 1:nrow(nest_top)) {
  cat(sprintf("    %s: rho=%.3f\n", nest_top$gene[i], nest_top$rho[i]))
}

# Summary counts
cat("\n  Gradient gene summary:\n")
cat(sprintf("    Strongly monotonic (|rho|>0.8): %d\n",
            sum(gradient_results$category == "strongly_monotonic")))
cat(sprintf("    Moderate (0.5<|rho|<=0.8): %d\n",
            sum(gradient_results$category == "moderate")))
cat(sprintf("    Core-enriched (rho<-0.5): %d\n",
            sum(gradient_results$rho < -0.5, na.rm = TRUE)))
cat(sprintf("    Nest-enriched (rho>0.5): %d\n",
            sum(gradient_results$rho > 0.5, na.rm = TRUE)))

# Save full table
write.csv(gradient_results, file.path(OUTPUT_DIR, "gradient_genes.csv"), row.names = FALSE)
cat("  Saved: gradient_genes.csv\n")

# =============================================================================
# STEP 4: PATHWAY ENRICHMENT
# =============================================================================
cat("\n--- STEP 4: Pathway Enrichment ---\n")

# Get Hallmark gene sets from msigdbr
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, gene_symbol) %>%
  distinct()

cat(sprintf("  Hallmark gene sets: %d\n", length(unique(hallmark_sets$gs_name))))

# CRITICAL: Background = CosMx panel genes (960), not whole genome
background_genes <- rownames(expr_data)
cat(sprintf("  Background universe: %d genes (CosMx panel)\n", length(background_genes)))

# Filter hallmark sets to only include genes in our panel
hallmark_filtered <- hallmark_sets %>%
  filter(gene_symbol %in% background_genes)
cat(sprintf("  Hallmark genes in panel: %d\n", length(unique(hallmark_filtered$gene_symbol))))

# Define core-enriched and nest-enriched gene lists
core_genes <- gradient_results %>% filter(rho < -0.5) %>% pull(gene)
nest_genes <- gradient_results %>% filter(rho > 0.5) %>% pull(gene)
cat(sprintf("  Core-enriched genes (rho < -0.5): %d\n", length(core_genes)))
cat(sprintf("  Nest-enriched genes (rho > 0.5): %d\n", length(nest_genes)))

# Fisher's exact test for enrichment
run_fisher_enrichment <- function(gene_list, hallmark_df, background) {
  pathways <- unique(hallmark_df$gs_name)
  results <- data.frame(
    pathway = character(),
    overlap = integer(),
    pathway_size_in_panel = integer(),
    gene_list_size = integer(),
    background_size = integer(),
    pvalue = numeric(),
    overlap_genes = character(),
    stringsAsFactors = FALSE
  )

  for (pw in pathways) {
    pw_genes <- hallmark_df %>% filter(gs_name == pw) %>% pull(gene_symbol)
    pw_in_panel <- pw_genes[pw_genes %in% background]

    if (length(pw_in_panel) < 2) next

    # 2x2 contingency table
    a <- length(intersect(gene_list, pw_in_panel))  # in list AND in pathway
    b <- length(setdiff(gene_list, pw_in_panel))     # in list NOT in pathway
    c <- length(setdiff(pw_in_panel, gene_list))     # in pathway NOT in list
    d <- length(background) - a - b - c              # neither

    if (a == 0) next

    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")

    overlap_genes_str <- paste(intersect(gene_list, pw_in_panel), collapse = "; ")

    results <- rbind(results, data.frame(
      pathway = pw,
      overlap = a,
      pathway_size_in_panel = length(pw_in_panel),
      gene_list_size = length(gene_list),
      background_size = length(background),
      pvalue = ft$p.value,
      overlap_genes = overlap_genes_str,
      stringsAsFactors = FALSE
    ))
  }

  if (nrow(results) > 0) {
    results$padj <- p.adjust(results$pvalue, method = "BH")
    results$neg_log10p <- -log10(results$pvalue)
    results <- results %>% arrange(pvalue)
  }

  return(results)
}

# Run enrichment for both directions
cat("\n  Running Fisher's exact tests...\n")
core_enrichment <- run_fisher_enrichment(core_genes, hallmark_filtered, background_genes)
nest_enrichment <- run_fisher_enrichment(nest_genes, hallmark_filtered, background_genes)

cat("\n  Top 10 pathways enriched in CORE-enriched genes:\n")
if (nrow(core_enrichment) > 0) {
  core_top10 <- head(core_enrichment, 10)
  for (i in 1:nrow(core_top10)) {
    cat(sprintf("    %s: overlap=%d, p=%.2e\n",
                gsub("HALLMARK_", "", core_top10$pathway[i]),
                core_top10$overlap[i],
                core_top10$pvalue[i]))
  }
} else {
  cat("    (none)\n")
}

cat("\n  Top 10 pathways enriched in NEST-enriched genes:\n")
if (nrow(nest_enrichment) > 0) {
  nest_top10 <- head(nest_enrichment, 10)
  for (i in 1:nrow(nest_top10)) {
    cat(sprintf("    %s: overlap=%d, p=%.2e\n",
                gsub("HALLMARK_", "", nest_top10$pathway[i]),
                nest_top10$overlap[i],
                nest_top10$pvalue[i]))
  }
} else {
  cat("    (none)\n")
}

# Save enrichment results
core_enrichment$direction <- "core-enriched"
nest_enrichment$direction <- "nest-enriched"
all_enrichment <- bind_rows(core_enrichment, nest_enrichment)
write.csv(all_enrichment, file.path(OUTPUT_DIR, "pathway_enrichment.csv"), row.names = FALSE)
cat("  Saved: pathway_enrichment.csv\n")

# =============================================================================
# STEP 5: CLUSTER GENES BY GRADIENT PROFILE
# =============================================================================
cat("\n--- STEP 5: Cluster Genes by Gradient Profile ---\n")

# Filter genes with |rho| > 0.4
gradient_genes <- gradient_results %>% filter(abs(rho) > 0.4) %>% pull(gene)
cat(sprintf("  Genes with |rho| > 0.4: %d\n", length(gradient_genes)))

if (length(gradient_genes) >= 10) {
  # Z-score normalize each gene's expression across zones
  gradient_expr <- gene_by_zone[gradient_genes, , drop = FALSE]
  gradient_z <- t(apply(gradient_expr, 1, function(x) {
    if (sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }))
  colnames(gradient_z) <- zone_levels

  # Hierarchical clustering
  dist_mat <- as.dist(1 - cor(t(gradient_z), method = "pearson"))
  # Handle NaN in distance matrix
  dist_mat[is.nan(dist_mat)] <- 2
  hc <- hclust(dist_mat, method = "ward.D2")

  # Cut into 5 clusters
  n_clusters <- min(5, length(gradient_genes))
  clusters <- cutree(hc, k = n_clusters)

  # Build cluster summary
  cluster_df <- data.frame(
    gene = names(clusters),
    cluster = clusters,
    stringsAsFactors = FALSE
  ) %>%
    left_join(gradient_results %>% select(gene, rho, direction, category), by = "gene")

  cat("\n  Cluster summary:\n")
  for (k in 1:n_clusters) {
    cl_genes <- cluster_df$gene[cluster_df$cluster == k]
    cl_z <- gradient_z[cl_genes, , drop = FALSE]
    mean_profile <- colMeans(cl_z, na.rm = TRUE)
    mean_rho <- mean(cluster_df$rho[cluster_df$cluster == k], na.rm = TRUE)

    # Assign descriptive label
    label <- case_when(
      mean_rho < -0.5 ~ "Core-high",
      mean_rho > 0.5  ~ "Nest-high",
      mean_profile[1] < -0.5 & mean_profile[length(zone_levels)] < -0.5 &
        max(mean_profile[2:(length(zone_levels)-1)]) > 0.3 ~ "Interface-peak",
      mean_rho < 0 ~ "Core-trending",
      TRUE ~ "Nest-trending"
    )

    cat(sprintf("    Cluster %d (%s): %d genes, mean rho=%.3f\n",
                k, label, length(cl_genes), mean_rho))
    cat(sprintf("      Mean z-score profile: %s\n",
                paste(sprintf("%.2f", mean_profile), collapse = " -> ")))
    if (length(cl_genes) <= 15) {
      cat(sprintf("      Genes: %s\n", paste(cl_genes, collapse = ", ")))
    } else {
      cat(sprintf("      First 15: %s ...\n", paste(head(cl_genes, 15), collapse = ", ")))
    }

    cluster_df$label[cluster_df$cluster == k] <- label
  }

  # Save cluster assignments
  write.csv(cluster_df, file.path(OUTPUT_DIR, "gene_clusters.csv"), row.names = FALSE)
  cat("  Saved: gene_clusters.csv\n")
} else {
  cat("  WARNING: <10 gradient genes found. Skipping clustering.\n")
  gradient_z <- gene_by_zone[gradient_genes, , drop = FALSE]
  cluster_df <- data.frame(gene = gradient_genes, cluster = 1, label = "all")
}

# =============================================================================
# STEP 6: FIGURES
# =============================================================================
cat("\n--- STEP 6: Generating Figures ---\n")

# Common theme for Nature Cell Biology style
theme_nature <- function(base_size = 8) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      text = element_text(family = "Helvetica"),
      axis.text = element_text(size = 7, color = "black"),
      axis.title = element_text(size = 8),
      plot.title = element_text(size = 9, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 7, hjust = 0),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7),
      legend.key.size = unit(0.3, "cm"),
      strip.text = element_text(size = 7, face = "bold"),
      strip.background = element_blank(),
      panel.grid = element_blank()
    )
}

# ---- Panel 1A: Spatial zone map ----
cat("  Panel 1A: Spatial zone map...\n")

plot_df_1a <- meta_cancer %>%
  select(x_um, y_um, Zone) %>%
  filter(!is.na(Zone))

# Use only colors for zones that exist
zone_colors_used <- ZONE_COLORS[zone_levels]

p1a <- ggplot(plot_df_1a, aes(x = x_um, y = y_um, color = Zone)) +
  geom_point(size = 0.3, alpha = 0.7) +
  scale_color_manual(values = zone_colors_used) +
  coord_fixed() +
  theme_nature() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), nrow = 1)) +
  labs(title = "Spatial zones - Patient D Bx3 cancer cells")

ggsave(file.path(OUTPUT_DIR, "Panel_1A_SpatialZones.pdf"),
       p1a, width = 5, height = 4.5)
cat("    Saved: Panel_1A_SpatialZones.pdf\n")

# ---- Panel 1B: 2x3 marker expression grid ----
cat("  Panel 1B: Marker expression spatial maps...\n")

marker_genes <- c("ESR1", "GATA3", "DUSP4", "DUSP6", "MKI67", "KRT8")
marker_genes_available <- marker_genes[marker_genes %in% rownames(expr_data)]
cat(sprintf("    Markers available: %s\n", paste(marker_genes_available, collapse = ", ")))

if (length(marker_genes_available) > 0) {
  # Get expression for each marker
  marker_plots <- list()
  for (gene in marker_genes_available) {
    expr_vec <- as.numeric(expr_data[gene, rownames(meta_cancer)])
    plot_df_1b <- data.frame(
      x = meta_cancer$x_um,
      y = meta_cancer$y_um,
      expr = expr_vec
    )

    # Cap at 99th percentile for visualization
    cap <- quantile(plot_df_1b$expr, 0.99, na.rm = TRUE)
    plot_df_1b$expr_capped <- pmin(plot_df_1b$expr, cap)

    marker_plots[[gene]] <- ggplot(plot_df_1b, aes(x = x, y = y, color = expr_capped)) +
      geom_point(size = 0.2, alpha = 0.7) +
      scale_color_viridis(option = "D", name = "Expr") +
      coord_fixed() +
      theme_nature() +
      theme(
        legend.position = "right",
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.15, "cm"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()
      ) +
      labs(title = gene)
  }

  p1b <- wrap_plots(marker_plots, ncol = 3) +
    plot_annotation(title = "Marker gene expression",
                    theme = theme(plot.title = element_text(size = 9, face = "bold")))

  ggsave(file.path(OUTPUT_DIR, "Panel_1B_MarkerExpression.pdf"),
         p1b, width = 8, height = 5)
  cat("    Saved: Panel_1B_MarkerExpression.pdf\n")
}

# ---- Panel 1C: Gradient heatmap (ComplexHeatmap) ----
cat("  Panel 1C: Gradient heatmap...\n")

if (length(gradient_genes) >= 10) {
  # Order genes by cluster, then by rho within cluster
  cluster_order <- cluster_df %>%
    arrange(cluster, rho) %>%
    pull(gene)

  heatmap_z <- gradient_z[cluster_order, , drop = FALSE]

  # Cap z-scores for visualization
  heatmap_z[heatmap_z > 2] <- 2
  heatmap_z[heatmap_z < -2] <- -2

  # Cluster annotation
  cluster_anno <- cluster_df %>%
    filter(gene %in% cluster_order) %>%
    arrange(match(gene, cluster_order))

  cluster_colors <- c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A",
                       "4" = "#984EA3", "5" = "#FF7F00")

  # Row annotation: cluster
  row_ha <- rowAnnotation(
    Cluster = factor(cluster_anno$cluster),
    col = list(Cluster = cluster_colors[1:n_clusters]),
    show_legend = TRUE,
    show_annotation_name = FALSE,
    width = unit(3, "mm")
  )

  # Column annotation: zone colors
  col_ha <- HeatmapAnnotation(
    Zone = zone_levels,
    col = list(Zone = zone_colors_used),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    simple_anno_size = unit(3, "mm")
  )

  # Color scale: RdBu diverging
  col_fun <- colorRamp2(c(-2, -1, 0, 1, 2),
                         c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"))

  ht <- Heatmap(
    heatmap_z,
    name = "Z-score",
    col = col_fun,
    cluster_rows = FALSE,  # already ordered by cluster
    cluster_columns = FALSE,
    show_row_names = nrow(heatmap_z) <= 60,
    row_names_gp = gpar(fontsize = 5, fontfamily = "Helvetica"),
    column_names_gp = gpar(fontsize = 7, fontfamily = "Helvetica"),
    column_names_rot = 45,
    left_annotation = row_ha,
    top_annotation = col_ha,
    column_title = "Spatial gradient: Core -> Nest",
    column_title_gp = gpar(fontsize = 9, fontface = "bold", fontfamily = "Helvetica"),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 7, fontfamily = "Helvetica"),
      labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica"),
      legend_height = unit(2, "cm")
    ),
    border = TRUE,
    use_raster = FALSE
  )

  pdf(file.path(OUTPUT_DIR, "Panel_1C_GradientHeatmap.pdf"),
      width = 4, height = max(4, min(12, nrow(heatmap_z) * 0.08)))
  draw(ht, padding = unit(c(2, 2, 2, 10), "mm"))
  dev.off()
  cat("    Saved: Panel_1C_GradientHeatmap.pdf\n")
}

# ---- Panel 1D: Topic score line plots ----
cat("  Panel 1D: Topic score line plots...\n")

# Compute topic gradients
topic_gradient <- data.frame(
  topic = topic_cols,
  rho = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(topic_cols)) {
  tc <- topic_cols[i]
  vals <- topic_by_zone[tc, ]
  if (sd(vals, na.rm = TRUE) > 0) {
    ct <- cor.test(zone_indices, vals, method = "spearman", exact = FALSE)
    topic_gradient$rho[i] <- ct$estimate
  } else {
    topic_gradient$rho[i] <- 0
  }
}

# Select topics with |rho| > 0.5, plus any key biological topics
selected_topics <- topic_gradient %>%
  filter(abs(rho) > 0.5) %>%
  pull(topic)

# Ensure we have at least some topics
if (length(selected_topics) < 3) {
  # Add top 3 by absolute rho
  extra <- topic_gradient %>%
    arrange(desc(abs(rho))) %>%
    head(3) %>%
    pull(topic)
  selected_topics <- unique(c(selected_topics, extra))
}

cat(sprintf("    Selected topics: %d\n", length(selected_topics)))

# Compute per-cell topic values with zone + SEM
topic_plot_data <- meta_cancer %>%
  select(Zone, all_of(selected_topics)) %>%
  pivot_longer(cols = all_of(selected_topics), names_to = "Topic", values_to = "Score") %>%
  group_by(Zone, Topic) %>%
  summarise(
    mean_score = mean(Score, na.rm = TRUE),
    se_score = sd(Score, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(Zone_Index = as.integer(Zone))

# Add rho annotation
topic_plot_data <- topic_plot_data %>%
  left_join(topic_gradient %>% rename(Topic = topic), by = "Topic") %>%
  mutate(Topic_Label = paste0(Topic, " (rho=", sprintf("%.2f", rho), ")"))

p1d <- ggplot(topic_plot_data, aes(x = Zone_Index, y = mean_score, color = Topic_Label)) +
  geom_ribbon(aes(ymin = mean_score - se_score, ymax = mean_score + se_score,
                  fill = Topic_Label), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.5) +
  scale_x_continuous(breaks = 1:length(zone_levels), labels = zone_levels) +
  scale_color_manual(values = rep(
    c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3"),
    length.out = length(selected_topics))) +
  scale_fill_manual(values = rep(
    c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3"),
    length.out = length(selected_topics))) +
  theme_nature() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 5)
  ) +
  guides(color = guide_legend(ncol = 2), fill = guide_legend(ncol = 2)) +
  labs(
    title = "Topic score gradients across spatial zones",
    x = NULL,
    y = "Mean topic score"
  )

ggsave(file.path(OUTPUT_DIR, "Panel_1D_TopicGradients.pdf"),
       p1d, width = 5.5, height = 4.5)
cat("    Saved: Panel_1D_TopicGradients.pdf\n")

# ---- Panel 1E: Pathway enrichment dot plot ----
cat("  Panel 1E: Pathway enrichment dot plot...\n")

# Combine top pathways from both directions
n_show <- 10
core_plot <- if (nrow(core_enrichment) > 0) {
  head(core_enrichment, n_show) %>%
    mutate(direction = "Core-enriched")
} else {
  data.frame()
}

nest_plot <- if (nrow(nest_enrichment) > 0) {
  head(nest_enrichment, n_show) %>%
    mutate(direction = "Nest-enriched")
} else {
  data.frame()
}

enrichment_plot_data <- bind_rows(core_plot, nest_plot)

if (nrow(enrichment_plot_data) > 0) {
  enrichment_plot_data <- enrichment_plot_data %>%
    mutate(
      pathway_short = gsub("HALLMARK_", "", pathway),
      pathway_short = gsub("_", " ", pathway_short),
      pathway_short = str_to_title(pathway_short),
      neg_log10p = -log10(pvalue)
    )

  p1e <- ggplot(enrichment_plot_data,
                aes(x = direction, y = reorder(pathway_short, neg_log10p))) +
    geom_point(aes(size = overlap, color = neg_log10p)) +
    scale_color_gradient(low = "grey70", high = "#D73027", name = expression(-log[10]*p)) +
    scale_size_continuous(range = c(2, 7), name = "Overlap") +
    theme_nature() +
    theme(
      legend.position = "right",
      axis.text.y = element_text(size = 6)
    ) +
    labs(
      title = "Pathway enrichment: Core vs Nest genes",
      subtitle = paste0("Background: ", length(background_genes), " CosMx panel genes"),
      x = NULL,
      y = NULL
    )

  ggsave(file.path(OUTPUT_DIR, "Panel_1E_PathwayEnrichment.pdf"),
         p1e, width = 5.5, height = 5)
  cat("    Saved: Panel_1E_PathwayEnrichment.pdf\n")
}

# ---- Panel 1F: QC box plot (nCount_RNA per zone) ----
cat("  Panel 1F: QC box plot...\n")

plot_df_1f <- meta_cancer %>%
  select(Zone, nCount_RNA) %>%
  filter(!is.na(Zone))

p1f <- ggplot(plot_df_1f, aes(x = Zone, y = nCount_RNA, fill = Zone)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, linewidth = 0.3) +
  scale_fill_manual(values = zone_colors_used) +
  theme_nature() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1, size = 6)
  ) +
  labs(
    title = "Transcript depth by spatial zone",
    x = NULL,
    y = "nCount_RNA"
  ) +
  ylim(0, NA)

ggsave(file.path(OUTPUT_DIR, "Panel_1F_QC_Depth.pdf"),
       p1f, width = 3.5, height = 3)
cat("    Saved: Panel_1F_QC_Depth.pdf\n")

# ---- Combined multi-panel layout ----
cat("  Combining panels...\n")

# Row 1: 1A + 1B
# Row 2: 1C + 1D
# Row 3: 1E + 1F

# For the combined figure, we use ggplot-based panels (1A, 1B, 1D, 1E, 1F)
# Panel 1C (ComplexHeatmap) is saved separately as it uses grid graphics

# Build combined layout of ggplot panels
combined_top <- plot_grid(p1a, NULL, ncol = 2, rel_widths = c(1, 1.5),
                          labels = c("a", "b"), label_size = 10)
# Note: panel 1B is patchwork, so we'll combine it differently

combined_bottom <- plot_grid(
  p1f + theme(plot.margin = margin(5, 10, 5, 5)),
  NULL,
  ncol = 2, rel_widths = c(0.4, 0.6),
  labels = c("f", ""), label_size = 10
)

# Save what we can combine (1A + 1D + 1E + 1F)
p_combined_gg <- plot_grid(
  p1a, p1d,
  if (exists("p1e")) p1e else NULL, p1f,
  ncol = 2,
  labels = c("a", "d", "e", "f"),
  label_size = 10,
  rel_widths = c(1, 1),
  rel_heights = c(1, 0.8)
)

ggsave(file.path(OUTPUT_DIR, "Combined_Panels_ggplot.pdf"),
       p_combined_gg, width = 10, height = 8)
cat("    Saved: Combined_Panels_ggplot.pdf\n")

# =============================================================================
# VERIFICATION
# =============================================================================
cat("\n--- VERIFICATION ---\n")

# 1. Check zone cell counts
cat("  1. Zone cell counts:\n")
print(table(meta_cancer$Zone))
if (any(table(meta_cancer$Zone) == 0)) {
  cat("  FAIL: Empty zones detected!\n")
} else {
  cat("  PASS: All zones have cells\n")
}

# 2. QC depth bias
cat(sprintf("  2. Depth ratio: %.2f ", depth_ratio))
if (depth_ratio <= 2) cat("PASS\n") else cat("WARNING\n")

# 3. Check expected markers
expected_core <- c("ESR1", "GATA3")
expected_nest <- c()  # will depend on data
for (g in expected_core) {
  if (g %in% gradient_results$gene) {
    g_rho <- gradient_results$rho[gradient_results$gene == g]
    cat(sprintf("  3. %s rho=%.3f (expected core-enriched) ", g, g_rho))
    if (g_rho < 0) cat("PASS\n") else cat("CHECK\n")
  }
}

# 4. Background size
cat(sprintf("  4. Pathway background size: %d ", length(background_genes)))
if (length(background_genes) < 2000) cat("PASS (CosMx panel)\n") else cat("WARNING (too large?)\n")

# 5. Output files
cat("  5. Output files:\n")
output_files <- list.files(OUTPUT_DIR, recursive = FALSE)
for (f in output_files) {
  fsize <- file.size(file.path(OUTPUT_DIR, f))
  cat(sprintf("    %s (%.1f KB)\n", f, fsize / 1024))
}

cat("\n=", rep("=", 69), "\n", sep = "")
cat("DONE - All outputs in:", OUTPUT_DIR, "\n")
cat("=", rep("=", 69), "\n", sep = "")
