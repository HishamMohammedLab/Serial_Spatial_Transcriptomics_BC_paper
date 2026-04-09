# Figure 2 Patient A: Spatial Maps with POLYGONS - REVISED GROUPINGS
# Purpose: Visualize cell polygons with 4 functional niches

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(data.table)

# Output directory
output_dir <- "Figure2_PatientA"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# POLYGON FILE PATHS
# ============================================================================

polygon_dir <- "supplementary_input_data/Polygons/RNA"

polygon_files <- list(
  PatA_Bx2 = file.path(polygon_dir, "R1134_265303-polygons.csv"),
  PatA_Bx3 = file.path(polygon_dir, "R1124_321955-polygons.csv"),
  PatA_Bx4 = file.path(polygon_dir, "R1134_322078-polygons.csv")
)

# Verify files exist
for(name in names(polygon_files)) {
  if(!file.exists(polygon_files[[name]])) {
    stop("File not found: ", polygon_files[[name]])
  }
  message("Found: ", name, " -> ", basename(polygon_files[[name]]))
}

# ============================================================================
# DEFINE FUNCTIONAL NICHE GROUPINGS (REVISED - 4 GROUPS)
# ============================================================================

# Merged malignant groups (Treatment Responsive was artifact)

domain_groupings <- data.frame(
  spatial_domain = as.character(0:9),
  group_name = c(
    "Malignant",                  # 0 - was Clone3+ Cancer
    "Malignant",                  # 1 - was Treatment Responsive (merged)
    "Diploid Epithelial",         # 2
    "Malignant",                  # 3 - was Treatment Responsive (merged)
    "Tumor-Stroma Interface",     # 4
    "Malignant",                  # 5 - was Clone3+ Cancer
    "Malignant",                  # 6 - was Clone3+ Cancer
    "Fibroblast-Dominant Niche",  # 7 - was Liver Stroma
    "Fibroblast-Dominant Niche",  # 8 - was Liver Stroma
    "Diploid Epithelial"          # 9
  ),
  stringsAsFactors = FALSE
)

# Classy color palette (4 colors)
group_colors <- c(
  "Malignant" = "#C44E52",                 # Muted brick red
  "Diploid Epithelial" = "#DD8452",        # Warm terracotta
  "Tumor-Stroma Interface" = "#8172B3",    # Dusty purple
  "Fibroblast-Dominant Niche" = "#55A868"  # Forest green
)

# Define order for legend
group_order <- c("Malignant", "Diploid Epithelial",
                 "Tumor-Stroma Interface", "Fibroblast-Dominant Niche")

clone_colors <- c(
  "PatA_Clone1" = "#FF7F0E",
  "PatA_Clone2" = "#2CA02C",
  "PatA_Clone3" = "#D62728",
  "PatA_Diploid" = "#7F7F7F"
)

stagate_colors <- c(
  "0" = "#E41A1C", "1" = "#377EB8", "2" = "#4DAF4A", "3" = "#984EA3",
  "4" = "#FF7F00", "5" = "#FFFF33", "6" = "#A65628", "7" = "#F781BF",
  "8" = "#999999", "9" = "#66C2A5"
)

# ============================================================================
# LOAD SEURAT METADATA
# ============================================================================

message("\nLoading Seurat object for metadata...")
obj <- readRDS("data/CosMx_SMMART_345k_clean.rds")

# Extract Patient A metadata
patA <- subset(obj, Patient == "Patient_A")
meta <- patA@meta.data
meta$cell_barcode <- rownames(meta)

# Parse cell ID to match polygon format (fov_cellID)
meta <- meta %>%
  mutate(
    fov = as.integer(gsub("c_[0-9]+_([0-9]+)_.*", "\\1", cell_barcode)),
    cellID = as.integer(gsub("c_[0-9]+_[0-9]+_([0-9]+)", "\\1", cell_barcode)),
    spatial_domain = as.character(spatial_domain)
  ) %>%
  left_join(domain_groupings, by = "spatial_domain")

# Set factor levels for consistent ordering
meta$group_name <- factor(meta$group_name, levels = group_order)

message("Metadata prepared: ", nrow(meta), " cells")
message("Timepoints: ", paste(unique(meta$Timepoint), collapse = ", "))

# Print new group distribution
message("\n=== REVISED NICHE DISTRIBUTION ===")
print(table(meta$group_name, meta$Timepoint, useNA = "ifany"))

# ============================================================================
# FUNCTION TO LOAD AND PROCESS POLYGONS
# ============================================================================

load_polygons <- function(polygon_file, metadata, timepoint) {

  message("\nLoading polygons for ", timepoint, "...")

  polys <- fread(polygon_file)
  message("  Raw polygon vertices: ", nrow(polys))

  meta_tp <- metadata %>% filter(Timepoint == timepoint)
  message("  Cells in metadata: ", nrow(meta_tp))

  polys$match_key <- paste(polys$fov, polys$cellID, sep = "_")
  meta_tp$match_key <- paste(meta_tp$fov, meta_tp$cellID, sep = "_")

  message("  Matching cells: ", length(unique(polys$match_key[polys$match_key %in% meta_tp$match_key])))

  polys_annotated <- polys %>%
    left_join(
      meta_tp %>% select(match_key, spatial_domain, group_name, PatA_Clone, predicted.id),
      by = "match_key"
    )

  polys_annotated <- polys_annotated %>%
    group_by(fov, cellID) %>%
    mutate(vertex_order = row_number()) %>%
    ungroup()

  # Ensure factor levels are preserved
  polys_annotated$group_name <- factor(polys_annotated$group_name, levels = group_order)

  message("  Annotated vertices: ", sum(!is.na(polys_annotated$group_name)))

  return(polys_annotated)
}

# ============================================================================
# LOAD ALL POLYGON DATA
# ============================================================================

polys_bx2 <- load_polygons(polygon_files$PatA_Bx2, meta, "Bx2")
polys_bx3 <- load_polygons(polygon_files$PatA_Bx3, meta, "Bx3")
polys_bx4 <- load_polygons(polygon_files$PatA_Bx4, meta, "Bx4")

# ============================================================================
# CLEAN POLYGON PLOTTING FUNCTION (NO LABELS)
# ============================================================================

plot_polygons_clean <- function(poly_data, color_var, colors,
                                fill_alpha = 0.8, line_size = 0.02,
                                show_legend = FALSE) {

  plot_data <- poly_data %>% filter(!is.na(.data[[color_var]]))

  if(nrow(plot_data) == 0) {
    return(NULL)
  }

  message("Plotting ", length(unique(paste(plot_data$fov, plot_data$cellID))), " cells")

  p <- ggplot(plot_data, aes(x = x_global_px, y = y_global_px,
                              group = interaction(fov, cellID),
                              fill = .data[[color_var]])) +
    geom_polygon(color = NA, linewidth = line_size, alpha = fill_alpha) +
    scale_fill_manual(values = colors, na.value = "grey80", drop = FALSE) +
    coord_fixed() +
    theme_void() +
    theme(
      legend.position = if(show_legend) "right" else "none",
      plot.margin = margin(5, 5, 5, 5)
    )

  return(p)
}

# ============================================================================
# GENERATE CLEAN INDIVIDUAL PLOTS (NO LABELS)
# ============================================================================

message("\n=== GENERATING CLEAN INDIVIDUAL PLOTS ===\n")

# Bx2 (Liver - different met site)
p_bx2 <- plot_polygons_clean(polys_bx2, "group_name", group_colors)
if(!is.null(p_bx2)) {
  ggsave(file.path(output_dir, "PatientA_Bx2_SpatialNiches_Clean_v2.pdf"),
         p_bx2, width = 8, height = 8)
  message("Saved: PatientA_Bx2_SpatialNiches_Clean_v2.pdf")
}

# Bx3 (Bone Marrow)
p_bx3 <- plot_polygons_clean(polys_bx3, "group_name", group_colors)
if(!is.null(p_bx3)) {
  ggsave(file.path(output_dir, "PatientA_Bx3_SpatialNiches_Clean_v2.pdf"),
         p_bx3, width = 8, height = 8)
  message("Saved: PatientA_Bx3_SpatialNiches_Clean_v2.pdf")
}

# Bx4 (Liver - same site as Bx2)
p_bx4 <- plot_polygons_clean(polys_bx4, "group_name", group_colors)
if(!is.null(p_bx4)) {
  ggsave(file.path(output_dir, "PatientA_Bx4_SpatialNiches_Clean_v2.pdf"),
         p_bx4, width = 8, height = 8)
  message("Saved: PatientA_Bx4_SpatialNiches_Clean_v2.pdf")
}

# ============================================================================
# CREATE SEPARATE LEGEND
# ============================================================================

message("\nCreating legend...")

# Create a dummy plot just for the legend
legend_data <- data.frame(
  x = 1:4,
  y = 1:4,
  group = factor(group_order, levels = group_order)
)

p_legend <- ggplot(legend_data, aes(x = x, y = y, fill = group)) +
  geom_point(shape = 22, size = 8, stroke = 0.5) +
  scale_fill_manual(values = group_colors, name = "Spatial Niche") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1.2, "cm")
  ) +
  guides(fill = guide_legend(override.aes = list(size = 6)))

# Extract just the legend
library(cowplot)
legend_only <- get_legend(p_legend)

# Save legend as separate file
pdf(file.path(output_dir, "PatientA_SpatialNiches_Legend_v2.pdf"), width = 4, height = 4)
grid::grid.draw(legend_only)
dev.off()
message("Saved: PatientA_SpatialNiches_Legend_v2.pdf")

# ============================================================================
# GENERATE COMBINED FIGURES (WITH LABELS)
# ============================================================================

message("\n=== GENERATING COMBINED FIGURES ===\n")

# Function for plots with legend
plot_polygons_labeled <- function(poly_data, color_var, colors, title,
                                   fill_alpha = 0.8, line_size = 0.02) {

  plot_data <- poly_data %>% filter(!is.na(.data[[color_var]]))

  if(nrow(plot_data) == 0) return(NULL)

  p <- ggplot(plot_data, aes(x = x_global_px, y = y_global_px,
                              group = interaction(fov, cellID),
                              fill = .data[[color_var]])) +
    geom_polygon(color = NA, linewidth = line_size, alpha = fill_alpha) +
    scale_fill_manual(values = colors, na.value = "grey80",
                      name = "Spatial Niche", drop = FALSE) +
    coord_fixed() +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    ) +
    labs(title = title)

  return(p)
}

# All timepoints comparison
p_groups_bx2 <- plot_polygons_labeled(polys_bx2, "group_name", group_colors, "Bx2 (Liver)")
p_groups_bx3 <- plot_polygons_labeled(polys_bx3, "group_name", group_colors, "Bx3 (Bone Marrow)")
p_groups_bx4 <- plot_polygons_labeled(polys_bx4, "group_name", group_colors, "Bx4 (Liver)")

if(!is.null(p_groups_bx2) && !is.null(p_groups_bx3) && !is.null(p_groups_bx4)) {
  fig_all <- (p_groups_bx2 | p_groups_bx3 | p_groups_bx4) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Patient A: Spatial Niches Across Treatment",
      subtitle = "4 functional niches based on cell composition and molecular signatures"
    )
  ggsave(file.path(output_dir, "PatientA_SpatialNiches_AllTimepoints_v2.pdf"),
         fig_all, width = 22, height = 8)
  message("Saved: PatientA_SpatialNiches_AllTimepoints_v2.pdf")
}

# Liver comparison (Bx2 vs Bx4 - same tissue type)
if(!is.null(p_groups_bx2) && !is.null(p_groups_bx4)) {
  fig_liver <- (p_groups_bx2 | p_groups_bx4) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Patient A: Liver Metastases Comparison (Bx2 vs Bx4)"
    )
  ggsave(file.path(output_dir, "PatientA_SpatialNiches_Bx2vsBx4_Liver_v2.pdf"),
         fig_liver, width = 16, height = 8)
  message("Saved: PatientA_SpatialNiches_Bx2vsBx4_Liver_v2.pdf")
}

# ============================================================================
# SUMMARY STATISTICS FOR NEW GROUPINGS
# ============================================================================

message("\n=== NICHE SUMMARY STATISTICS ===\n")

# Cell counts by niche and timepoint
niche_summary <- meta %>%
  filter(!is.na(group_name)) %>%
  group_by(Timepoint, group_name) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  pivot_wider(names_from = Timepoint, values_from = n_cells, values_fill = 0)

message("Cell counts by niche:")
print(niche_summary)

# Calculate percentages
niche_pct <- meta %>%
  filter(!is.na(group_name)) %>%
  group_by(Timepoint, group_name) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Timepoint) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  select(-n) %>%
  pivot_wider(names_from = Timepoint, values_from = pct, values_fill = 0)

message("\nPercentages by niche:")
print(niche_pct)

# Clone distribution within Malignant niche
clone_in_malignant <- meta %>%
  filter(group_name == "Malignant", !is.na(PatA_Clone)) %>%
  group_by(Timepoint, PatA_Clone) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Timepoint) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup()

message("\nClone distribution within Malignant niche:")
print(clone_in_malignant %>% filter(PatA_Clone == "PatA_Clone3"))

# Save summaries
write.csv(niche_summary,
          file.path(output_dir, "PatientA_SpatialNiche_Summary_v2.csv"),
          row.names = FALSE)

write.csv(niche_pct,
          file.path(output_dir, "PatientA_SpatialNiche_Percentages_v2.csv"),
          row.names = FALSE)

message("\n=== SCRIPT COMPLETE ===")
message("Output directory: ", output_dir)
message("\nFiles generated:")
message("  - PatientA_Bx2_SpatialNiches_Clean_v2.pdf")
message("  - PatientA_Bx3_SpatialNiches_Clean_v2.pdf")
message("  - PatientA_Bx4_SpatialNiches_Clean_v2.pdf")
message("  - PatientA_SpatialNiches_Legend_v2.pdf")
message("  - PatientA_SpatialNiches_AllTimepoints_v2.pdf")
message("  - PatientA_SpatialNiches_Bx2vsBx4_Liver_v2.pdf")
