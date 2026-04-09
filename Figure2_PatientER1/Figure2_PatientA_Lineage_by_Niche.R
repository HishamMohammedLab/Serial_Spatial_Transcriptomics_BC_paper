# Figure 2: Lineage Composition by Spatial Niche
# Bar plots showing percent lineages in each niche for Bx2 and Bx4 (both Liver)

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Output directory
output_dir <- "Figure2_PatientA"

# ============================================================================
# DEFINE NICHE GROUPINGS (4 GROUPS)
# ============================================================================

domain_groupings <- data.frame(
  spatial_domain = as.character(0:9),
  group_name = c(
    "Malignant",                  # 0
    "Malignant",                  # 1
    "Diploid Epithelial",         # 2
    "Malignant",                  # 3
    "Tumor-Stroma Interface",     # 4
    "Malignant",                  # 5
    "Malignant",                  # 6
    "Fibroblast-Dominant Niche",  # 7
    "Fibroblast-Dominant Niche",  # 8
    "Diploid Epithelial"          # 9
  ),
  stringsAsFactors = FALSE
)

niche_order <- c("Malignant", "Diploid Epithelial",
                 "Tumor-Stroma Interface", "Fibroblast-Dominant Niche")

niche_colors <- c(
  "Malignant" = "#C44E52",
  "Diploid Epithelial" = "#DD8452",
  "Tumor-Stroma Interface" = "#8172B3",
  "Fibroblast-Dominant Niche" = "#55A868"
)

# ============================================================================
# LOAD DATA
# ============================================================================

message("Loading Seurat object...")
obj <- readRDS("data/CosMx_SMMART_345k_clean.rds")

# Extract Patient A metadata
patA <- subset(obj, Patient == "Patient_A")
meta <- patA@meta.data
meta$cell_id <- rownames(meta)
meta$spatial_domain <- as.character(meta$spatial_domain)

# Add niche groupings
meta <- meta %>%
  left_join(domain_groupings, by = "spatial_domain")

meta$group_name <- factor(meta$group_name, levels = niche_order)

# Check what lineage column exists
message("\nAvailable lineage columns:")
lineage_cols <- grep("lineage|predicted|cell.?type", names(meta), ignore.case = TRUE, value = TRUE)
print(lineage_cols)

# Use predicted.id as lineage
message("\nLineage categories (predicted.id):")
print(table(meta$predicted.id))

# ============================================================================
# CALCULATE LINEAGE PERCENTAGES BY NICHE
# ============================================================================

calc_lineage_pct <- function(data, timepoint_name) {
  data %>%
    filter(Timepoint == timepoint_name, !is.na(group_name), !is.na(predicted.id)) %>%
    group_by(group_name, predicted.id) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(group_name) %>%
    mutate(
      total = sum(n),
      pct = 100 * n / total
    ) %>%
    ungroup() %>%
    mutate(Timepoint = timepoint_name)
}

lineage_bx2 <- calc_lineage_pct(meta, "Bx2")
lineage_bx4 <- calc_lineage_pct(meta, "Bx4")

message("\n=== Bx2 Lineage by Niche ===")
print(lineage_bx2 %>% arrange(group_name, desc(pct)))

message("\n=== Bx4 Lineage by Niche ===")
print(lineage_bx4 %>% arrange(group_name, desc(pct)))

# ============================================================================
# DEFINE LINEAGE COLORS
# ============================================================================

# Get all unique lineages
all_lineages <- unique(c(lineage_bx2$predicted.id, lineage_bx4$predicted.id))
message("\nAll lineages: ", paste(all_lineages, collapse = ", "))

# Define a color palette for lineages
lineage_colors <- c(
  "Cancer Epithelial" = "#E41A1C",
  "Normal Epithelial" = "#FF7F00",
  "CAFs" = "#4DAF4A",
  "Endothelial" = "#377EB8",
  "PVL" = "#984EA3",
  "Myeloid" = "#FFFF33",
  "T-cells" = "#A65628",
  "B-cells" = "#F781BF",
  "Plasmablasts" = "#999999"
)

# Add any missing lineages with default colors
missing_lineages <- setdiff(all_lineages, names(lineage_colors))
if(length(missing_lineages) > 0) {
  extra_colors <- scales::hue_pal()(length(missing_lineages))
  names(extra_colors) <- missing_lineages
  lineage_colors <- c(lineage_colors, extra_colors)
}

# ============================================================================
# CREATE BAR PLOTS
# ============================================================================

create_lineage_barplot <- function(data, title) {

  # Order lineages by overall frequency
  lineage_order <- data %>%
    group_by(predicted.id) %>%
    summarise(total = sum(n)) %>%
    arrange(desc(total)) %>%
    pull(predicted.id)

  data$predicted.id <- factor(data$predicted.id, levels = rev(lineage_order))

  ggplot(data, aes(x = group_name, y = pct, fill = predicted.id)) +
    geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
    scale_fill_manual(values = lineage_colors, name = "Cell Lineage") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
    labs(
      title = title,
      x = "Spatial Niche",
      y = "Percentage of Cells"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9)
    )
}

# Create individual plots
p_bx2 <- create_lineage_barplot(lineage_bx2, "Bx2 (Liver)")
p_bx4 <- create_lineage_barplot(lineage_bx4, "Bx4 (Liver)")

# Save individual plots
ggsave(file.path(output_dir, "PatientA_Bx2_Lineage_by_Niche.pdf"),
       p_bx2, width = 8, height = 6)
message("Saved: PatientA_Bx2_Lineage_by_Niche.pdf")

ggsave(file.path(output_dir, "PatientA_Bx4_Lineage_by_Niche.pdf"),
       p_bx4, width = 8, height = 6)
message("Saved: PatientA_Bx4_Lineage_by_Niche.pdf")

# Combined figure
combined <- (p_bx2 | p_bx4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Patient A: Cell Lineage Composition by Spatial Niche (Liver Biopsies)",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

ggsave(file.path(output_dir, "PatientA_Lineage_by_Niche_Bx2_Bx4.pdf"),
       combined, width = 14, height = 6)
message("Saved: PatientA_Lineage_by_Niche_Bx2_Bx4.pdf")

# ============================================================================
# SAVE DATA TABLE
# ============================================================================

lineage_all <- bind_rows(lineage_bx2, lineage_bx4)
write.csv(lineage_all,
          file.path(output_dir, "PatientA_Lineage_by_Niche_Data.csv"),
          row.names = FALSE)
message("Saved: PatientA_Lineage_by_Niche_Data.csv")

message("\n=== COMPLETE ===")
