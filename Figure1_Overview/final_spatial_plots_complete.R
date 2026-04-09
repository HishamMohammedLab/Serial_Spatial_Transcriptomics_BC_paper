# ==============================================================================
# SMMART: FINAL SPATIAL PLOTS - ALL PATIENTS
# Separate files per biopsy + Zoom plots per FOV
# ==============================================================================

library(Seurat)
library(tidyverse)
library(sf)
library(patchwork)

set.seed(42)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Base output directory
BASE_OUTPUT_DIR <- "Topic_UMAP_Publication/LR_Analysis/Spatial_Plots/"
dir.create(BASE_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Polygon file directory
POLYGON_DIR <- "supplementary_input_data/Polygons/RNA/"

# Polygon file mapping for ALL patients
POLYGON_MAPPING <- list(
  # Patient A
  "Patient_A_Bx2" = "R1134_265303-polygons.csv",
  "Patient_A_Bx3" = "R1124_321955-polygons.csv",
  "Patient_A_Bx4" = "R1134_322078-polygons.csv",
  # Patient B
  "Patient_B_Bx1" = "R1124_265607-polygons.csv",
  "Patient_B_Bx4" = "R1124_322118-polygons.csv",
  # Patient C
  "Patient_C_Bx1" = "R1134_272830-polygons.csv",
  "Patient_C_Bx2" = "R1134_303148-polygons.csv",
  # Patient D
  "Patient_D_Bx1" = "R1134_272840-polygons.csv",
  "Patient_D_Bx2" = "R1134_321920-polygons.csv",
  "Patient_D_Bx3" = "R1124_322091-polygons.csv"
)

# Seurat object path
SEURAT_PATH <- "data/CosMx_SMMART_345k_clean.rds"

# ==============================================================================
# COLOR PALETTES
# ==============================================================================

# Lineage colors (Retro Muted)
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
  "Pericyte"                = "#A5799D",
  "Adipocyte"               = "#5D8E4D"
)

# STAGATE Domain colors
DOMAIN_COLORS <- c(
  "0" = "#D51F26",
  "1" = "#A6CEE3",
  "2" = "#FB9A99",
  "3" = "#D95F02",
  "4" = "#FDBF6F",
  "5" = "#7570B3",
  "6" = "#CAB2D6",
  "7" = "#B2DF8A",
  "8" = "#FFFF99",
  "9" = "#B15928"
)

# ==============================================================================
# COMMON THEME
# ==============================================================================

common_theme <- theme_void(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# ==============================================================================
# LOAD SEURAT DATA
# ==============================================================================

cat("========== LOADING DATA ==========\n\n")

cat("Loading Seurat object...\n")
obj_all <- readRDS(SEURAT_PATH)
obj_all <- UpdateSeuratObject(obj_all)

cat("Total cells:", ncol(obj_all), "\n")

# Extract metadata
meta <- obj_all@meta.data
meta$cell_id <- rownames(meta)

# Parse cell IDs
parts <- strsplit(meta$cell_id, "_")
meta$seurat_sample <- sapply(parts, function(x) as.integer(x[2]))
meta$seurat_fov <- sapply(parts, function(x) as.integer(x[3]))
meta$seurat_cellID <- sapply(parts, function(x) as.integer(x[4]))

# Format columns
meta$Lineage <- as.character(meta$Lineage)
meta$Domain <- as.character(meta$spatial_domain)

cat("\nData loaded successfully.\n")

# ==============================================================================
# FUNCTION: Create SF polygons (GLOBAL coordinates for whole-slide)
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
# FUNCTION: Create SF polygons (LOCAL coordinates for zoom)
# ==============================================================================

create_sf_local <- function(poly_df, fov_num) {
  
  fov_data <- poly_df %>% filter(fov == fov_num)
  if (nrow(fov_data) == 0) return(NULL)
  
  cell_polys <- fov_data %>%
    group_by(cellID) %>%
    filter(n() >= 3) %>%
    summarise(
      geometry = list(tryCatch({
        coords <- cbind(x_local_px, y_local_px)
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
  
  st_sf(cellID = cell_polys$cellID, geometry = st_sfc(cell_polys$geometry))
}

# ==============================================================================
# FUNCTION: Create Lineage Plot (whole-slide)
# ==============================================================================

create_lineage_plot_global <- function(sf_polys, biopsy_meta, title_text) {
  
  sf_polys$match_key <- paste0(sf_polys$fov, "_", sf_polys$cellID)
  meta_lookup <- biopsy_meta %>%
    mutate(match_key = paste0(seurat_fov, "_", seurat_cellID)) %>%
    select(match_key, Lineage)
  
  sf_matched <- sf_polys %>% 
    left_join(meta_lookup, by = "match_key")
  
  sf_with_data <- sf_matched %>% filter(!is.na(Lineage))
  sf_no_data <- sf_matched %>% filter(is.na(Lineage))
  
  p <- ggplot() +
    geom_sf(data = sf_no_data, fill = "grey85", color = "grey30", linewidth = 0.01) +
    geom_sf(data = sf_with_data, aes(fill = Lineage), color = "grey30", linewidth = 0.01) +
    scale_fill_manual(values = LINEAGE_COLORS, name = "Lineage", drop = FALSE) +
    labs(title = title_text) +
    common_theme +
    guides(fill = guide_legend(override.aes = list(size = 3)))
  
  return(p)
}

# ==============================================================================
# FUNCTION: Create Domain Plot (whole-slide)
# ==============================================================================

create_domain_plot_global <- function(sf_polys, biopsy_meta, title_text) {
  
  sf_polys$match_key <- paste0(sf_polys$fov, "_", sf_polys$cellID)
  meta_lookup <- biopsy_meta %>%
    mutate(match_key = paste0(seurat_fov, "_", seurat_cellID)) %>%
    select(match_key, Domain)
  
  sf_matched <- sf_polys %>% 
    left_join(meta_lookup, by = "match_key")
  
  sf_with_data <- sf_matched %>% filter(!is.na(Domain))
  sf_no_data <- sf_matched %>% filter(is.na(Domain))
  
  p <- ggplot() +
    geom_sf(data = sf_no_data, fill = "grey85", color = "grey30", linewidth = 0.01) +
    geom_sf(data = sf_with_data, aes(fill = Domain), color = "grey30", linewidth = 0.01) +
    scale_fill_manual(values = DOMAIN_COLORS, name = "Domain", drop = FALSE) +
    labs(title = title_text) +
    common_theme +
    guides(fill = guide_legend(override.aes = list(size = 3)))
  
  return(p)
}

# ==============================================================================
# FUNCTION: Create Lineage Plot (zoom - local coordinates)
# ==============================================================================

create_lineage_plot_zoom <- function(sf_polys, fov_meta, title_text) {
  
  sf_polys$match_key <- as.character(sf_polys$cellID)
  meta_lookup <- fov_meta %>%
    mutate(match_key = as.character(seurat_cellID)) %>%
    select(match_key, Lineage)
  
  sf_matched <- sf_polys %>% 
    left_join(meta_lookup, by = "match_key")
  
  sf_with_data <- sf_matched %>% filter(!is.na(Lineage))
  sf_no_data <- sf_matched %>% filter(is.na(Lineage))
  
  p <- ggplot() +
    geom_sf(data = sf_no_data, fill = "grey85", color = "grey30", linewidth = 0.05) +
    geom_sf(data = sf_with_data, aes(fill = Lineage), color = "grey30", linewidth = 0.05) +
    scale_fill_manual(values = LINEAGE_COLORS, name = "Lineage", drop = FALSE) +
    labs(title = title_text) +
    common_theme +
    guides(fill = guide_legend(override.aes = list(size = 3)))
  
  return(p)
}

# ==============================================================================
# FUNCTION: Create Domain Plot (zoom - local coordinates)
# ==============================================================================

create_domain_plot_zoom <- function(sf_polys, fov_meta, title_text) {
  
  sf_polys$match_key <- as.character(sf_polys$cellID)
  meta_lookup <- fov_meta %>%
    mutate(match_key = as.character(seurat_cellID)) %>%
    select(match_key, Domain)
  
  sf_matched <- sf_polys %>% 
    left_join(meta_lookup, by = "match_key")
  
  sf_with_data <- sf_matched %>% filter(!is.na(Domain))
  sf_no_data <- sf_matched %>% filter(is.na(Domain))
  
  p <- ggplot() +
    geom_sf(data = sf_no_data, fill = "grey85", color = "grey30", linewidth = 0.05) +
    geom_sf(data = sf_with_data, aes(fill = Domain), color = "grey30", linewidth = 0.05) +
    scale_fill_manual(values = DOMAIN_COLORS, name = "Domain", drop = FALSE) +
    labs(title = title_text) +
    common_theme +
    guides(fill = guide_legend(override.aes = list(size = 3)))
  
  return(p)
}

# ==============================================================================
# PROCESS EACH PATIENT
# ==============================================================================

cat("\n========== PROCESSING PATIENTS ==========\n\n")

patients <- c("Patient_A", "Patient_B", "Patient_C", "Patient_D")

for (patient in patients) {
  
  cat("\n######################################################\n")
  cat("Processing", patient, "\n")
  cat("######################################################\n")
  
  # Create patient subfolder
  patient_dir <- file.path(BASE_OUTPUT_DIR, patient)
  dir.create(patient_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create zoom subfolder
  zoom_dir <- file.path(patient_dir, "Zoom_FOVs")
  dir.create(zoom_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Get biopsies for this patient
  patient_biopsies <- names(POLYGON_MAPPING)[grepl(patient, names(POLYGON_MAPPING))]
  patient_biopsies <- sort(patient_biopsies)
  
  if (length(patient_biopsies) == 0) {
    cat("  No biopsies found for", patient, "\n")
    next
  }
  
  cat("  Biopsies:", paste(patient_biopsies, collapse = ", "), "\n")
  cat("  Output folder:", patient_dir, "\n")
  
  # Filter metadata to this patient
  patient_meta <- meta %>% filter(Patient == patient)
  cat("  Total cells:", nrow(patient_meta), "\n")
  
  # Process each biopsy
  for (biopsy_key in patient_biopsies) {
    
    # Parse biopsy info
    biopsy_parts <- strsplit(biopsy_key, "_")[[1]]
    timepoint <- biopsy_parts[3]  # e.g., "Bx2"
    
    cat("\n  ================================================\n")
    cat("  Processing", biopsy_key, "\n")
    cat("  ================================================\n")
    
    # Load polygon file
    poly_path <- file.path(POLYGON_DIR, POLYGON_MAPPING[[biopsy_key]])
    
    if (!file.exists(poly_path)) {
      cat("    WARNING: Polygon file not found:", poly_path, "\n")
      next
    }
    
    poly_df <- read.csv(poly_path)
    cat("    Polygon file loaded:", nrow(poly_df), "vertices\n")
    
    # Get metadata for this biopsy
    biopsy_meta <- patient_meta %>% filter(Timepoint == timepoint)
    cat("    Cells in metadata:", nrow(biopsy_meta), "\n")
    
    # Get available FOVs
    fovs_available <- sort(unique(poly_df$fov))
    cat("    FOVs available:", length(fovs_available), "\n")
    
    # =========================================================================
    # WHOLE-SLIDE PLOTS
    # =========================================================================
    
    cat("\n    --- Creating whole-slide plots ---\n")
    
    sf_whole <- create_sf_global(poly_df)
    
    if (!is.null(sf_whole) && nrow(sf_whole) > 0) {
      
      # Lineage whole-slide
      p_lineage <- create_lineage_plot_global(
        sf_whole, 
        biopsy_meta, 
        title_text = paste(patient, timepoint, "- Lineage")
      )
      
      ggsave(
        file.path(patient_dir, paste0(timepoint, "_Lineage.pdf")),
        p_lineage,
        width = 10, height = 8, bg = "white"
      )
      cat("    ✓", paste0(timepoint, "_Lineage.pdf"), "\n")
      
      # Domain whole-slide
      p_domain <- create_domain_plot_global(
        sf_whole, 
        biopsy_meta, 
        title_text = paste(patient, timepoint, "- STAGATE Domain")
      )
      
      ggsave(
        file.path(patient_dir, paste0(timepoint, "_Domain.pdf")),
        p_domain,
        width = 10, height = 8, bg = "white"
      )
      cat("    ✓", paste0(timepoint, "_Domain.pdf"), "\n")
      
    } else {
      cat("    WARNING: Could not create whole-slide polygons\n")
    }
    
    # =========================================================================
    # ZOOM PLOTS (per FOV)
    # =========================================================================
    
    cat("\n    --- Creating zoom plots for", length(fovs_available), "FOVs ---\n")
    
    for (fov_num in fovs_available) {
      
      # Create local polygons for this FOV
      sf_zoom <- create_sf_local(poly_df, fov_num)
      
      if (is.null(sf_zoom) || nrow(sf_zoom) == 0) {
        cat("      FOV", fov_num, ": No polygons (skipped)\n")
        next
      }
      
      # Get FOV-specific metadata
      fov_meta <- biopsy_meta %>% filter(seurat_fov == fov_num)
      
      # Lineage zoom
      p_lineage_zoom <- create_lineage_plot_zoom(
        sf_zoom, 
        fov_meta, 
        title_text = paste(patient, timepoint, "FOV", fov_num, "- Lineage")
      )
      
      ggsave(
        file.path(zoom_dir, paste0(timepoint, "_FOV", fov_num, "_Lineage.pdf")),
        p_lineage_zoom,
        width = 8, height = 8, bg = "white"
      )
      
      # Domain zoom
      p_domain_zoom <- create_domain_plot_zoom(
        sf_zoom, 
        fov_meta, 
        title_text = paste(patient, timepoint, "FOV", fov_num, "- Domain")
      )
      
      ggsave(
        file.path(zoom_dir, paste0(timepoint, "_FOV", fov_num, "_Domain.pdf")),
        p_domain_zoom,
        width = 8, height = 8, bg = "white"
      )
    }
    
    cat("    ✓ Zoom plots saved to:", zoom_dir, "\n")
    
    # Clean up
    rm(poly_df, sf_whole)
    gc()
  }
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n\n========== COMPLETE ==========\n")
cat("\nOutput structure:\n")
cat(BASE_OUTPUT_DIR, "\n")

for (patient in patients) {
  patient_biopsies <- names(POLYGON_MAPPING)[grepl(patient, names(POLYGON_MAPPING))]
  patient_biopsies <- sort(patient_biopsies)
  
  cat("  └──", patient, "/\n")
  
  for (bx in patient_biopsies) {
    timepoint <- strsplit(bx, "_")[[1]][3]
    cat("      ├──", paste0(timepoint, "_Lineage.pdf"), "\n")
    cat("      ├──", paste0(timepoint, "_Domain.pdf"), "\n")
  }
  
  cat("      └── Zoom_FOVs/\n")
  cat("          ├── Bx#_FOV#_Lineage.pdf\n")
  cat("          └── Bx#_FOV#_Domain.pdf\n")
}

cat("\nColor schemes:\n")
cat("  Lineage: Retro Muted palette\n")
cat("  Domain: Striking (0,3,5) vs Muted (1,7) vs Neutral (2,4,6,8,9)\n")
