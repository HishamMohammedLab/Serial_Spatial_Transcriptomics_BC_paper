# ==============================================================================
# SPATIAL DOMAIN VISUALIZATION: BX1 vs BX2
# Two versions: Script-style colors and custom zone colors
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

# Load SEPARATE polygon files
poly_bx1 <- read.csv("supplementary_input_data/Polygons/RNA/R1134_272830-polygons.csv")
poly_bx2 <- read.csv("supplementary_input_data/Polygons/RNA/R1134_303148-polygons.csv")

# Output directory
out_dir <- "PatientC_Analysis/CAF_Subtypes/Spatial_Viz/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# PREP FUNCTION
# ==============================================================================

prep_biopsy <- function(obj, pat_id, timepoint, poly_df) {
  
  sub_obj <- subset(obj, Patient == pat_id & Timepoint == timepoint)
  meta <- sub_obj@meta.data %>% mutate(cell_id = rownames(.))
  
  # Parse IDs
  parts <- strsplit(meta$cell_id, "_")
  meta$seurat_fov    <- sapply(parts, function(x) as.integer(x[3]))
  meta$seurat_cellID <- sapply(parts, function(x) as.integer(x[4]))
  meta$match_key     <- paste0(meta$seurat_fov, "_", meta$seurat_cellID)
  
  # Create zone classification (two groups)
  meta$Zone <- case_when(
    meta$spatial_domain %in% c("2", "3") ~ "Cancer Distal",
    TRUE ~ "Cancer Adjacent"
  )
  
  # Merge with polygon centroids
  poly_df$match_key <- paste0(poly_df$fov, "_", poly_df$cellID)
  
  centroids <- poly_df %>%
    group_by(match_key) %>%
    summarise(x = mean(x_global_px, na.rm = TRUE), 
              y = mean(y_global_px, na.rm = TRUE),
              .groups = "drop")
  
  spatial_data <- meta %>%
    inner_join(centroids, by = "match_key")
  
  return(spatial_data)
}

df_bx1 <- prep_biopsy(sobj, "Patient_C", "Bx1", poly_bx1)
df_bx2 <- prep_biopsy(sobj, "Patient_C", "Bx2", poly_bx2)

message("Bx1 cells with coords: ", nrow(df_bx1))
message("Bx2 cells with coords: ", nrow(df_bx2))

rm(sobj)
gc()

# ==============================================================================
# VERSION 1: SCRIPT-STYLE COLORS (Viridis/Magma aesthetic)
# ==============================================================================

message("\n=== VERSION 1: Script-style colors ===")

# Using a warm palette similar to the purity maps
plot_zones_v1 <- function(df, title, show_legend = FALSE) {
  
  # Color scheme matching script style
  zone_colors_v1 <- c("Cancer Adjacent" = "#E64B35", "Cancer Distal" = "#3C5488")
  
  ggplot(df, aes(x = x, y = y, color = Zone)) +
    geom_point(size = 0.4, alpha = 0.7) +
    scale_color_manual(values = zone_colors_v1, name = "Spatial Zone") +
    coord_fixed() +
    labs(title = title) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      legend.position = ifelse(show_legend, "right", "none")
    )
}

p_v1_bx1 <- plot_zones_v1(df_bx1, "Bx1: Mixed Zones")
p_v1_bx2 <- plot_zones_v1(df_bx2, "Bx2: Distinct Zones", show_legend = TRUE)

p_version1 <- p_v1_bx1 | p_v1_bx2
p_version1 <- p_version1 + plot_annotation(
  title = "Spatial Zones: Cancer Adjacent vs Cancer Distal",
  subtitle = "Adjacent (red) = D6 + Interface | Distal (blue) = D2/3",
  theme = theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11, color = "grey30")
  )
)

pdf(paste0(out_dir, "Fig_Zones_ScriptStyle.pdf"), width = 14, height = 6)
print(p_version1)
dev.off()
message("Saved: Fig_Zones_ScriptStyle.pdf")

# ==============================================================================
# VERSION 2: CUSTOM ZONE COLORS (Brown/Teal from earlier analysis)
# ==============================================================================

message("\n=== VERSION 2: Custom zone colors ===")

plot_zones_v2 <- function(df, title, show_legend = FALSE) {
  
  # Your established color scheme
 zone_colors_v2 <- c("Cancer Adjacent" = "#9E503C", "Cancer Distal" = "#2D5A5A")
  
  ggplot(df, aes(x = x, y = y, color = Zone)) +
    geom_point(size = 0.4, alpha = 0.7) +
    scale_color_manual(values = zone_colors_v2, name = "Spatial Zone") +
    coord_fixed() +
    labs(title = title) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      legend.position = ifelse(show_legend, "right", "none")
    )
}

p_v2_bx1 <- plot_zones_v2(df_bx1, "Bx1: Mixed Zones")
p_v2_bx2 <- plot_zones_v2(df_bx2, "Bx2: Distinct Zones", show_legend = TRUE)

p_version2 <- p_v2_bx1 | p_v2_bx2
p_version2 <- p_version2 + plot_annotation(
  title = "Spatial Zones: Cancer Adjacent vs Cancer Distal",
  subtitle = "Adjacent (brown) = D6 + Interface | Distal (teal) = D2/3",
  theme = theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11, color = "grey30")
  )
)

pdf(paste0(out_dir, "Fig_Zones_CustomColors.pdf"), width = 14, height = 6)
print(p_version2)
dev.off()
message("Saved: Fig_Zones_CustomColors.pdf")

# ==============================================================================
# VERSION 3: WITH DENSITY CONTOURS (More polished)
# ==============================================================================

message("\n=== VERSION 3: With density contours ===")

plot_zones_contours <- function(df, title, show_legend = FALSE) {
  
  zone_colors <- c("Cancer Adjacent" = "#9E503C", "Cancer Distal" = "#2D5A5A")
  
  ggplot() +
    # Background points
    geom_point(data = df, aes(x = x, y = y, color = Zone), 
               size = 0.3, alpha = 0.5) +
    # Density contours for each zone
    stat_density_2d(data = df %>% filter(Zone == "Cancer Adjacent"),
                    aes(x = x, y = y), color = "#6D3A2B", linewidth = 0.6) +
    stat_density_2d(data = df %>% filter(Zone == "Cancer Distal"),
                    aes(x = x, y = y), color = "#1A3D3D", linewidth = 0.6) +
    scale_color_manual(values = zone_colors, name = "Spatial Zone") +
    coord_fixed() +
    labs(title = title) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      legend.position = ifelse(show_legend, "right", "none")
    )
}

p_v3_bx1 <- plot_zones_contours(df_bx1, "Bx1: Overlapping Zones")
p_v3_bx2 <- plot_zones_contours(df_bx2, "Bx2: Separated Zones", show_legend = TRUE)

p_version3 <- p_v3_bx1 | p_v3_bx2
p_version3 <- p_version3 + plot_annotation(
  title = "Spatial Zone Territories with Density Contours",
  subtitle = "Contour lines show zone density boundaries",
  theme = theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11, color = "grey30")
  )
)

pdf(paste0(out_dir, "Fig_Zones_WithContours.pdf"), width = 14, height = 6)
print(p_version3)
dev.off()
message("Saved: Fig_Zones_WithContours.pdf")

# ==============================================================================
# VERSION 4: MINIMAL PUBLICATION STYLE (No legend, clean)
# ==============================================================================

message("\n=== VERSION 4: Minimal publication style ===")

plot_zones_minimal <- function(df, title) {
  
  zone_colors <- c("Cancer Adjacent" = "#9E503C", "Cancer Distal" = "#2D5A5A")
  
  ggplot(df, aes(x = x, y = y, color = Zone)) +
    geom_point(size = 0.3, alpha = 0.8) +
    scale_color_manual(values = zone_colors) +
    coord_fixed() +
    labs(title = title) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      legend.position = "none",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

p_v4_bx1 <- plot_zones_minimal(df_bx1, "Bx1")
p_v4_bx2 <- plot_zones_minimal(df_bx2, "Bx2")

p_version4 <- p_v4_bx1 | p_v4_bx2

pdf(paste0(out_dir, "Fig_Zones_Minimal.pdf"), width = 12, height = 5)
print(p_version4)
dev.off()
message("Saved: Fig_Zones_Minimal.pdf")

# ==============================================================================
# QUANTIFICATION: ZONE COMPOSITION SHIFT
# ==============================================================================

message("\n=== Zone Composition Summary ===")

zone_summary <- bind_rows(
  df_bx1 %>% mutate(Sample = "Bx1"),
  df_bx2 %>% mutate(Sample = "Bx2")
) %>%
  group_by(Sample, Zone) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(pct = n_cells / sum(n_cells) * 100)

message("\nZone composition:")
print(zone_summary)

# Bar plot of composition
p_comp <- ggplot(zone_summary, aes(x = Sample, y = pct, fill = Zone)) +
  geom_col(color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Cancer Adjacent" = "#9E503C", "Cancer Distal" = "#2D5A5A")) +
  labs(
    title = "Zone Composition Shift",
    y = "% of Cells",
    x = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "right"
  )

pdf(paste0(out_dir, "Fig_ZoneComposition_Bar.pdf"), width = 6, height = 5)
print(p_comp)
dev.off()
message("Saved: Fig_ZoneComposition_Bar.pdf")

# ==============================================================================
# COMBINED FIGURE
# ==============================================================================

message("\n=== Creating combined figure ===")

p_combined <- (p_version2) / (p_version3) / (p_comp + plot_spacer() + plot_layout(widths = c(1, 2)))

pdf(paste0(out_dir, "Fig_Zones_Combined.pdf"), width = 14, height = 16)
print(p_combined)
dev.off()
message("Saved: Fig_Zones_Combined.pdf")

# ==============================================================================
# SUMMARY
# ==============================================================================

message("\n")
message(strrep("=", 70))
message("OUTPUTS SAVED TO: ", out_dir)
message(strrep("=", 70))
message("\nFigures:")
message("  - Fig_Zones_ScriptStyle.pdf (Red/Blue - matching your script)")
message("  - Fig_Zones_CustomColors.pdf (Brown/Teal - your CAF colors)")
message("  - Fig_Zones_WithContours.pdf (With density contours)")
message("  - Fig_Zones_Minimal.pdf (Clean publication style)")
message("  - Fig_ZoneComposition_Bar.pdf (Quantification)")
message("  - Fig_Zones_Combined.pdf (All panels)")

message("\nZone summary:")
print(zone_summary)
