# STAGATE Domain Lineage Composition - Stacked Bar Plots
# One figure per patient showing lineage proportions by domain for each Bx

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(grid)

# =============================================================================
# OUTPUT DIRECTORY
# =============================================================================

output_dir <- "STAGATE_Domain_Lineage"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# =============================================================================
# DEFINE LINEAGE COLORS (Retro Muted Palette)
# =============================================================================

lineage_colors <- c(
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

# =============================================================================
# LOAD DATA
# =============================================================================

message("Loading Seurat object...")
obj <- readRDS("data/CosMx_SMMART_345k_clean.rds")

meta <- obj@meta.data
meta$cell_id <- rownames(meta)

# Check available columns
message("\nChecking metadata columns...")
message("Patients: ", paste(unique(meta$Patient), collapse = ", "))

# Check for spatial_domain or STAGATE_Domain column
domain_col <- if ("spatial_domain" %in% names(meta)) "spatial_domain" else "STAGATE_Domain"
message("Using domain column: ", domain_col)

# Check for Timepoint column
timepoint_col <- if ("Timepoint" %in% names(meta)) "Timepoint" else "Biopsy"
message("Using timepoint column: ", timepoint_col)

# Check for Lineage column
lineage_col <- if ("Lineage" %in% names(meta)) "Lineage" else "predicted.id"
message("Using lineage column: ", lineage_col)

# Standardize column names
meta$Domain <- as.character(meta[[domain_col]])
meta$Bx <- meta[[timepoint_col]]
meta$Lineage <- meta[[lineage_col]]

# Check lineage categories
message("\nLineage categories found:")
print(table(meta$Lineage))

# =============================================================================
# STANDARDIZE LINEAGE NAMES (if needed)
# =============================================================================

# Create mapping for common variations
lineage_mapping <- c(
  "cancer" = "Cancer",
  "Cancer Epithelial" = "Cancer",
  "epithelial" = "Epithelial",
  "Normal Epithelial" = "Epithelial",
  "tcell" = "T Lymphocyte",
  "T-cells" = "T Lymphocyte",
  "T Lymphocyte" = "T Lymphocyte",
  "bcell" = "B Lymphocyte",
  "B-cells" = "B Lymphocyte",
  "B Lymphocyte" = "B Lymphocyte",
  "myeloid" = "Myeloid",
  "Myeloid" = "Myeloid",
  "plasma" = "Plasma",
  "Plasma" = "Plasma",
  "Plasmablasts" = "Plasma",
  "fibroblast" = "Fibroblast",
  "Fibroblast" = "Fibroblast",
  "CAFs" = "Fibroblast",
  "endothelial_vascular" = "Endothelial - Vascular",
  "Endothelial - Vascular" = "Endothelial - Vascular",
  "Endothelial" = "Endothelial - Vascular",
  "endothelial_lymphatic" = "Endothelial - Lymphatic",
  "Endothelial - Lymphatic" = "Endothelial - Lymphatic",
  "adipocyte" = "Adipocyte",
  "Adipocyte" = "Adipocyte",
  "pericyte" = "Pericyte",
  "Pericyte" = "Pericyte",
  "PVL" = "Pericyte"
)

# Apply mapping
meta$Lineage_Std <- ifelse(meta$Lineage %in% names(lineage_mapping),
                            lineage_mapping[meta$Lineage],
                            meta$Lineage)

# Check standardized lineages
message("\nStandardized lineage categories:")
print(table(meta$Lineage_Std))

# Filter to only known lineages
known_lineages <- names(lineage_colors)
meta_filtered <- meta %>% filter(Lineage_Std %in% known_lineages)

message("\nCells retained after filtering: ", nrow(meta_filtered), " / ", nrow(meta))

# =============================================================================
# CALCULATE LINEAGE PROPORTIONS BY DOMAIN AND BX
# =============================================================================

calc_proportions <- function(data) {
  data %>%
    filter(!is.na(Domain), !is.na(Bx), !is.na(Lineage_Std)) %>%
    group_by(Patient, Bx, Domain, Lineage_Std) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(Patient, Bx, Domain) %>%
    mutate(
      total = sum(n),
      pct = 100 * n / total
    ) %>%
    ungroup()
}

all_props <- calc_proportions(meta_filtered)

# =============================================================================
# CREATE PLOTTING FUNCTION
# =============================================================================

create_domain_barplot <- function(data, patient_name, bx_name, lineage_cols, show_legend = FALSE) {

  plot_data <- data %>%
    filter(Patient == patient_name, Bx == bx_name)

  if (nrow(plot_data) == 0) return(NULL)

  # Order domains numerically
  plot_data$Domain <- factor(plot_data$Domain, levels = as.character(0:9))

  # Order lineages by overall frequency (most common at bottom)
  lineage_order <- plot_data %>%
    group_by(Lineage_Std) %>%
    summarise(total = sum(n)) %>%
    arrange(desc(total)) %>%
    pull(Lineage_Std)

  plot_data$Lineage_Std <- factor(plot_data$Lineage_Std, levels = rev(lineage_order))

  p <- ggplot(plot_data, aes(x = Domain, y = pct, fill = Lineage_Std)) +
    geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
    scale_fill_manual(values = lineage_cols, name = "Lineage") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
    labs(
      title = bx_name,
      x = "STAGATE Domain",
      y = "Percentage of Cells"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.title = element_text(size = 22),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12)
    )

  if (show_legend) {
    p <- p + theme(legend.position = "right")
  } else {
    p <- p + theme(legend.position = "none")
  }

  return(p)
}

# =============================================================================
# HELPER: Extract legend from a ggplot
# =============================================================================

get_legend <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if (length(leg) > 0) {
    return(tmp$grobs[[leg]])
  }
  return(NULL)
}

# =============================================================================
# GENERATE PLOTS FOR EACH PATIENT
# =============================================================================

patients <- unique(all_props$Patient)

for (pat in patients) {

  message("\n=== Processing ", pat, " ===")

  pat_data <- all_props %>% filter(Patient == pat)
  timepoints <- sort(unique(pat_data$Bx))

  message("Timepoints: ", paste(timepoints, collapse = ", "))

  # Create a plot for each timepoint (without legends)
  plot_list <- list()
  for (tp in timepoints) {
    p <- create_domain_barplot(all_props, pat, tp, lineage_colors, show_legend = FALSE)
    if (!is.null(p)) {
      plot_list[[tp]] <- p
    }
  }

  if (length(plot_list) == 0) {
    message("No plots generated for ", pat)
    next
  }

  # Combine plots
  n_plots <- length(plot_list)

  if (n_plots == 1) {
    plots_combined <- plot_list[[1]]
    fig_width <- 10
    fig_height <- 6
  } else if (n_plots == 2) {
    plots_combined <- (plot_list[[1]] | plot_list[[2]])
    fig_width <- 14
    fig_height <- 6
  } else if (n_plots == 3) {
    plots_combined <- (plot_list[[1]] | plot_list[[2]] | plot_list[[3]])
    fig_width <- 18
    fig_height <- 6
  } else {
    # 4 or more - arrange in 2 rows
    plots_combined <- wrap_plots(plot_list, ncol = 2)
    fig_width <- 14
    fig_height <- 10
  }

  # For Patient A and D: no legend; for others: add shared legend
  if (pat %in% c("Patient_A", "Patient_D")) {
    combined <- plots_combined +
      plot_annotation(
        title = paste0(pat, ": Lineage Composition by STAGATE Domain"),
        theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
      )
  } else {
    # Create one plot with legend to extract it
    legend_plot <- create_domain_barplot(all_props, pat, timepoints[1], lineage_colors, show_legend = TRUE)
    shared_legend <- get_legend(legend_plot)

    # Add shared legend on the right
    combined <- plots_combined + wrap_elements(shared_legend) +
      plot_layout(widths = c(rep(1, n_plots), 0.3)) +
      plot_annotation(
        title = paste0(pat, ": Lineage Composition by STAGATE Domain"),
        theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
      )

    # Adjust width to account for legend
    fig_width <- fig_width + 2
  }

  # Save PDF only
  filename <- gsub("Patient_", "", pat)
  ggsave(file.path(output_dir, paste0(filename, "_STAGATE_Lineage_ByDomain.pdf")),
         combined, width = fig_width, height = fig_height)
  message("Saved: ", filename, "_STAGATE_Lineage_ByDomain.pdf")
}

# =============================================================================
# SAVE DATA TABLE
# =============================================================================

write.csv(all_props, file.path(output_dir, "All_Patients_Lineage_by_Domain_Data.csv"),
          row.names = FALSE)
message("\nSaved: All_Patients_Lineage_by_Domain_Data.csv")

message("\n=== COMPLETE ===")
message("Output directory: ", output_dir)
