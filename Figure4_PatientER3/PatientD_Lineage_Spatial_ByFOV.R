# PatientD_Lineage_Spatial_ByFOV.R
# Spatial polygon plots colored by lineage for Patient D Bx3 — one per FOV

library(Seurat)
library(tidyverse)
library(sf)

# ============================================================
# FILE PATHS
# ============================================================
SEURAT_PATH <- "data/CosMx_SMMART_345k_clean.rds"
POLYGON_DIR <- "supplementary_input_data/Polygons/RNA/"
POLYGON_FILE <- "R1124_322091-polygons.csv"
OUTPUT_DIR <- "PatientD_Analysis/LR_Spatial/PatD_Bx3_Lineage_ByFOV/"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# COLOR SCHEME
# ============================================================
LINEAGE_COLORS <- c(
  "Cancer"                  = "#82B0AF",
  "Epithelial"              = "#E69138",
  "T Lymphocyte"            = "#E9C44D",
  "B Lymphocyte"            = "#D95F5F",
  "Myeloid"                 = "#B7B0A8",
  "Plasma"                  = "#8B426F",
  "Fibroblast"              = "#A47C66",
  "Endothelial - Vascular"  = "#5D85AD",
  "Endothelial - Lymphatic" = "#F198A4",
  "Adipocyte"               = "#5D8E4D",
  "Pericyte"                = "#A5799D"
)

# ============================================================
# LOAD SEURAT OBJECT
# ============================================================
cat("Loading Seurat object...\n")
so <- readRDS(SEURAT_PATH)

meta <- so@meta.data
meta$cell_id <- rownames(meta)

d_bx3 <- meta %>% filter(Patient == "Patient_D", Timepoint == "Bx3")
cat(sprintf("Patient D Bx3: %d cells\n", nrow(d_bx3)))

# Parse cell IDs
parts <- strsplit(d_bx3$cell_id, "_")
d_bx3$seurat_fov <- sapply(parts, function(x) as.integer(x[3]))
d_bx3$seurat_cellID <- sapply(parts, function(x) as.integer(x[4]))
d_bx3$match_key <- paste0(d_bx3$seurat_fov, "_", d_bx3$seurat_cellID)

# ============================================================
# LOAD AND CREATE POLYGONS
# ============================================================
cat("Loading polygon data...\n")
poly_df <- read.csv(file.path(POLYGON_DIR, POLYGON_FILE))

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

  st_sf(fov = cell_polys$fov, cellID = cell_polys$cellID,
        geometry = st_sfc(cell_polys$geometry))
}

sf_polys <- create_sf_global(poly_df)
sf_polys$match_key <- paste0(sf_polys$fov, "_", sf_polys$cellID)
cat(sprintf("Created %d polygons\n", nrow(sf_polys)))

# Merge
sf_merged <- sf_polys %>%
  left_join(d_bx3 %>% select(match_key, Lineage),
            by = "match_key") %>%
  filter(!is.na(Lineage))
cat(sprintf("Matched %d cells\n", nrow(sf_merged)))

# ============================================================
# PLOT EACH FOV
# ============================================================
fovs <- sort(unique(sf_merged$fov))
cat(sprintf("\n%d FOVs to plot: %s\n", length(fovs), paste(fovs, collapse = ", ")))

for (f in fovs) {
  cat(sprintf("\n--- FOV %d ---\n", f))

  sf_fov <- sf_merged %>% filter(fov == f)
  cat(sprintf("  Cells: %d\n", nrow(sf_fov)))

  p <- ggplot() +
    geom_sf(data = sf_fov, aes(fill = Lineage),
            color = "black", linewidth = 0.01) +
    scale_fill_manual(values = LINEAGE_COLORS, name = "Lineage") +
    theme_void() +
    theme(
      plot.title = element_text(family = "Helvetica", size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(family = "Helvetica", size = 10, hjust = 0.5, color = "grey40"),
      legend.text = element_text(family = "Helvetica", size = 8),
      legend.title = element_text(family = "Helvetica", size = 9, face = "bold"),
      legend.position = "right"
    ) +
    labs(
      title = sprintf("Patient D Bx3 | FOV %d", f),
      subtitle = "Cell Lineage"
    )

  output_path <- file.path(OUTPUT_DIR, sprintf("FOV%02d_Lineage.pdf", f))
  ggsave(output_path, p, width = 10, height = 8)
  cat(sprintf("  Saved: %s\n", output_path))
}

cat("\nAll FOVs done!\n")
