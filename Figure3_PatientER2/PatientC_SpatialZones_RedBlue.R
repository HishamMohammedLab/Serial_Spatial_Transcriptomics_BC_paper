# ==============================================================================
# SPATIAL ZONE - RED/BLUE SHADE VARIATIONS
# ==============================================================================

library(Seurat)
library(tidyverse)
library(patchwork)

set.seed(42)

out_dir <- "PatientC_SpatialZones/"

# ==============================================================================
# LOAD DATA
# ==============================================================================

message("Loading data...")

sobj <- readRDS("data/CosMx_SMMART_345k_clean.rds")
poly_bx1 <- read.csv("supplementary_input_data/Polygons/RNA/R1134_272830-polygons.csv")
poly_bx2 <- read.csv("supplementary_input_data/Polygons/RNA/R1134_303148-polygons.csv")

prep_biopsy <- function(obj, pat_id, timepoint, poly_df) {
  sub_obj <- subset(obj, Patient == pat_id & Timepoint == timepoint)
  meta <- sub_obj@meta.data %>% mutate(cell_id = rownames(.))
  parts <- strsplit(meta$cell_id, "_")
  meta$seurat_fov    <- sapply(parts, function(x) as.integer(x[3]))
  meta$seurat_cellID <- sapply(parts, function(x) as.integer(x[4]))
  meta$match_key     <- paste0(meta$seurat_fov, "_", meta$seurat_cellID)
  meta$Zone <- case_when(
    meta$spatial_domain %in% c("2", "3") ~ "Cancer Distal",
    TRUE ~ "Cancer Adjacent"
  )
  poly_df$match_key <- paste0(poly_df$fov, "_", poly_df$cellID)
  centroids <- poly_df %>%
    group_by(match_key) %>%
    summarise(x = mean(x_global_px, na.rm = TRUE),
              y = mean(y_global_px, na.rm = TRUE),
              .groups = "drop")
  meta %>% inner_join(centroids, by = "match_key")
}

df_bx1 <- prep_biopsy(sobj, "Patient_C", "Bx1", poly_bx1)
df_bx2 <- prep_biopsy(sobj, "Patient_C", "Bx2", poly_bx2)

rm(sobj); gc()

# ==============================================================================
# RED/BLUE PALETTE VARIATIONS
# ==============================================================================

palettes <- list(
  # Original
  p1 = list(Adjacent = "#E64B35", Distal = "#3C5488", name = "RB1_Original"),

  # Darker/more saturated
  p2 = list(Adjacent = "#C41E3A", Distal = "#1B3A6D", name = "RB2_Dark_Cardinal_Navy"),

  # Lighter/softer
  p3 = list(Adjacent = "#E57373", Distal = "#64B5F6", name = "RB3_Light_Coral_SkyBlue"),

  # Muted/desaturated
  p4 = list(Adjacent = "#CD5C5C", Distal = "#4A6FA5", name = "RB4_Muted_IndianRed_SlateBlue"),

  # Vibrant/pure
  p5 = list(Adjacent = "#FF0000", Distal = "#0000FF", name = "RB5_Pure_Red_Blue")
)

# ==============================================================================
# PLOT FUNCTION
# ==============================================================================

plot_minimal <- function(df, title, colors) {
  ggplot(df, aes(x = x, y = y, color = Zone)) +
    geom_point(size = 0.3, alpha = 0.8) +
    scale_color_manual(values = c("Cancer Adjacent" = colors$Adjacent,
                                   "Cancer Distal" = colors$Distal)) +
    coord_fixed() +
    labs(title = title) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      legend.position = "none",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

# ==============================================================================
# GENERATE INDIVIDUAL PDFs
# ==============================================================================

message("\n=== Generating Red/Blue variations ===\n")

for (pal in palettes) {
  p_bx1 <- plot_minimal(df_bx1, "Bx1", pal)
  p_bx2 <- plot_minimal(df_bx2, "Bx2", pal)

  p_combined <- (p_bx1 | p_bx2) +
    plot_annotation(
      title = pal$name,
      subtitle = paste0("Adjacent: ", pal$Adjacent, " | Distal: ", pal$Distal),
      theme = theme(
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 10, color = "grey40")
      )
    )

  filename <- paste0(out_dir, "ColorTest_", pal$name, ".pdf")
  pdf(filename, width = 12, height = 5)
  print(p_combined)
  dev.off()
  message("Saved: ", pal$name, ".pdf")
}

# ==============================================================================
# COMPARISON GRID - ALL 5 ON ONE PAGE
# ==============================================================================

message("\n=== Creating comparison grid ===\n")

plot_small <- function(df, pal) {
  ggplot(df, aes(x = x, y = y, color = Zone)) +
    geom_point(size = 0.15, alpha = 0.7) +
    scale_color_manual(values = c("Cancer Adjacent" = pal$Adjacent,
                                   "Cancer Distal" = pal$Distal)) +
    coord_fixed() +
    labs(title = pal$name,
         subtitle = paste0(pal$Adjacent, " / ", pal$Distal)) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 9, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 7, hjust = 0.5, color = "grey50"),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

# Create grid with Bx2 (shows better zone separation)
plots_list <- lapply(palettes, function(pal) plot_small(df_bx2, pal))

grid_plot <- wrap_plots(plots_list, ncol = 5) +
  plot_annotation(
    title = "Red/Blue Palette Comparison (Bx2)",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

pdf(paste0(out_dir, "ColorTest_RedBlue_Grid.pdf"), width = 20, height = 5)
print(grid_plot)
dev.off()
message("Saved: ColorTest_RedBlue_Grid.pdf")

# ==============================================================================
# SUMMARY
# ==============================================================================

message("\n")
message(strrep("=", 70))
message("RED/BLUE VARIATIONS SAVED TO: ", out_dir)
message(strrep("=", 70))
message("\nPalettes:")
for (pal in palettes) {
  message("  ", pal$name, ": Adjacent=", pal$Adjacent, " Distal=", pal$Distal)
}
message("\nGrid comparison: ColorTest_RedBlue_Grid.pdf")

