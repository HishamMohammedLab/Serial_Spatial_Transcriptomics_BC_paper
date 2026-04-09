# Figure 2 Patient A: Spatial Maps with Functional Domain Groupings
# Created: January 20, 2026
# Purpose: Visualize STAGATE domains colored by new functional groupings

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(viridis)

# Output directory
output_dir <- "Figure2_PatientA"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# DEFINE FUNCTIONAL DOMAIN GROUPINGS
# ============================================================================

domain_groupings <- data.frame(
  spatial_domain = as.character(0:9),
  group_name = c(
    "Clone3+ Cancer",        # 0: Was 66% Clone3 at Bx2, treatment response
    "Treatment Responsive",  # 1: Immune_Rich, dramatically shrank
    "Diploid Epithelial",    # 2: 75% Diploid, may be normal tissue
    "Treatment Responsive",  # 3: Essentially eliminated (43 cells at Bx4)
    "Tumor-Stroma Interface", # 4: Highest immune checkpoint
    "Clone3+ Cancer",        # 5: Highest ECM, IFN, L-R activity
    "Clone3+ Cancer",        # 6: Myeloid-dominated, high Clone3
    "Liver Stroma",          # 7: High stromal, FOV 27-30
    "Liver Stroma",          # 8: High stromal, FOV 27-30
    "Diploid Epithelial"     # 9: 69% Diploid
  ),
  stringsAsFactors = FALSE
)

# Define colors
group_colors <- c(
  "Clone3+ Cancer" = "#D62728",           # Red - true malignancy
  "Diploid Epithelial" = "#FF9896",       # Light red/pink
  "Tumor-Stroma Interface" = "#9467BD",   # Purple
  "Treatment Responsive" = "#1F77B4",     # Blue
  "Liver Stroma" = "#2CA02C"              # Green
)

clone_colors <- c(
  "PatA_Clone1" = "#FF7F0E",   # Orange - 8q gain
  "PatA_Clone2" = "#2CA02C",   # Green - 1q gain
  "PatA_Clone3" = "#D62728",   # Red - 11q13 amp (TRUE CANCER)
  "PatA_Diploid" = "#7F7F7F"   # Grey - normal
)

stagate_colors <- c(
  "0" = "#E41A1C", "1" = "#377EB8", "2" = "#4DAF4A", "3" = "#984EA3",
  "4" = "#FF7F00", "5" = "#FFFF33", "6" = "#A65628", "7" = "#F781BF",
  "8" = "#999999", "9" = "#66C2A5"
)

# ============================================================================
# LOAD DATA
# ============================================================================

message("Loading Seurat object...")
obj <- readRDS("data/CosMx_SMMART_345k_clean.rds")
message("Total cells: ", ncol(obj))

# Subset to Patient A
patA <- subset(obj, Patient == "Patient_A")
message("Patient A cells: ", ncol(patA))

# Extract metadata
meta <- patA@meta.data
meta$cell_id <- rownames(meta)
meta$spatial_domain <- as.character(meta$spatial_domain)

# Map to functional groups
meta <- meta %>%
  left_join(domain_groupings, by = "spatial_domain")

message("\nGroup distribution:")
print(table(meta$group_name, meta$Timepoint, useNA = "ifany"))

# ============================================================================
# CREATE SPATIAL PLOTS
# ============================================================================

# Function for spatial plot
plot_spatial <- function(data, color_var, colors, title, point_size = 0.1) {
  ggplot(data, aes(x = x_global_px, y = y_global_px, color = .data[[color_var]])) +
    geom_point(size = point_size, alpha = 0.6) +
    scale_color_manual(values = colors, na.value = "grey80") +
    coord_fixed() +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    ) +
    labs(title = title, color = gsub("_", " ", color_var))
}

# ============================================================================
# GENERATE PLOTS BY TIMEPOINT
# ============================================================================

message("\nGenerating spatial plots...")

# Split by timepoint
bx2_data <- meta %>% filter(Timepoint == "Bx2")
bx3_data <- meta %>% filter(Timepoint == "Bx3")
bx4_data <- meta %>% filter(Timepoint == "Bx4")

# --- FUNCTIONAL GROUPS ---

p_groups_bx2 <- plot_spatial(bx2_data, "group_name", group_colors,
                              "Bx2 - Functional Groups\n(Different Met Site)")
p_groups_bx3 <- plot_spatial(bx3_data, "group_name", group_colors,
                              "Bx3 - Functional Groups")
p_groups_bx4 <- plot_spatial(bx4_data, "group_name", group_colors,
                              "Bx4 - Functional Groups")

# --- STAGATE ORIGINAL ---

p_stagate_bx2 <- plot_spatial(bx2_data, "spatial_domain", stagate_colors,
                               "Bx2 - STAGATE Domains")
p_stagate_bx3 <- plot_spatial(bx3_data, "spatial_domain", stagate_colors,
                               "Bx3 - STAGATE Domains")
p_stagate_bx4 <- plot_spatial(bx4_data, "spatial_domain", stagate_colors,
                               "Bx4 - STAGATE Domains")

# --- CLONE DISTRIBUTION ---

p_clone_bx2 <- plot_spatial(bx2_data, "PatA_Clone", clone_colors,
                             "Bx2 - Clone (11q13)")
p_clone_bx3 <- plot_spatial(bx3_data, "PatA_Clone", clone_colors,
                             "Bx3 - Clone (11q13)")
p_clone_bx4 <- plot_spatial(bx4_data, "PatA_Clone", clone_colors,
                             "Bx4 - Clone (11q13)")

# ============================================================================
# COMBINED FIGURES
# ============================================================================

# 1. STAGATE vs Functional Groups - Bx4 (side by side comparison)
comparison_bx4 <- (p_stagate_bx4 | p_groups_bx4) +
  plot_annotation(
    title = "Patient A Bx4: STAGATE Domains vs Functional Groupings",
    subtitle = "Left: Original 10 STAGATE domains | Right: 5 Functional groups based on Clone3, L-R, signatures"
  ) &
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(output_dir, "PatientA_Bx4_STAGATE_vs_FunctionalGroups.pdf"),
       comparison_bx4, width = 16, height = 8)
message("Saved: PatientA_Bx4_STAGATE_vs_FunctionalGroups.pdf")

# 2. All timepoints - Functional Groups
all_groups <- (p_groups_bx2 | p_groups_bx3 | p_groups_bx4) +
  plot_annotation(
    title = "Patient A: Functional Domain Groups Across Treatment",
    subtitle = "Bx2 (different met site) | Bx3 | Bx4 (same site, post-treatment)"
  )

ggsave(file.path(output_dir, "PatientA_FunctionalGroups_AllTimepoints.pdf"),
       all_groups, width = 20, height = 7)
message("Saved: PatientA_FunctionalGroups_AllTimepoints.pdf")

# 3. Same site comparison (Bx3 vs Bx4)
same_site <- (p_groups_bx3 | p_groups_bx4) +
  plot_annotation(
    title = "Patient A: Same Metastatic Site Comparison (Bx3 → Bx4)",
    subtitle = "Functional domain groups showing treatment response"
  )

ggsave(file.path(output_dir, "PatientA_FunctionalGroups_Bx3vsBx4.pdf"),
       same_site, width = 14, height = 7)
message("Saved: PatientA_FunctionalGroups_Bx3vsBx4.pdf")

# 4. Clone distribution vs Functional Groups (Bx4)
clone_vs_groups <- (p_clone_bx4 | p_groups_bx4) +
  plot_annotation(
    title = "Patient A Bx4: Clone Distribution vs Functional Groups",
    subtitle = "Clone3 (red) = 11q13 amplification marker | Mapped to Clone3+ Cancer domains"
  )

ggsave(file.path(output_dir, "PatientA_Bx4_Clone_vs_FunctionalGroups.pdf"),
       clone_vs_groups, width = 16, height = 8)
message("Saved: PatientA_Bx4_Clone_vs_FunctionalGroups.pdf")

# 5. Three-row comparison: STAGATE | Groups | Clone (Bx4)
triple_comparison <- (p_stagate_bx4 / p_groups_bx4 / p_clone_bx4) +
  plot_annotation(
    title = "Patient A Bx4: Three Views of Spatial Organization",
    subtitle = "Top: STAGATE domains | Middle: Functional groups | Bottom: Clone (11q13 amp)"
  )

ggsave(file.path(output_dir, "PatientA_Bx4_TripleComparison.pdf"),
       triple_comparison, width = 10, height = 18)
message("Saved: PatientA_Bx4_TripleComparison.pdf")

# 6. All timepoints - STAGATE original
all_stagate <- (p_stagate_bx2 | p_stagate_bx3 | p_stagate_bx4) +
  plot_annotation(title = "Patient A: Original STAGATE Domains")

ggsave(file.path(output_dir, "PatientA_STAGATE_AllTimepoints.pdf"),
       all_stagate, width = 20, height = 7)
message("Saved: PatientA_STAGATE_AllTimepoints.pdf")

# 7. All timepoints - Clone distribution
all_clone <- (p_clone_bx2 | p_clone_bx3 | p_clone_bx4) +
  plot_annotation(
    title = "Patient A: Clone Distribution Across Treatment",
    subtitle = "Clone3 (red) = TRUE cancer with 11q13 amplification"
  )

ggsave(file.path(output_dir, "PatientA_Clone_AllTimepoints.pdf"),
       all_clone, width = 20, height = 7)
message("Saved: PatientA_Clone_AllTimepoints.pdf")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

# Group summary by timepoint
group_summary <- meta %>%
  filter(!is.na(group_name)) %>%
  group_by(Timepoint, group_name) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  pivot_wider(names_from = Timepoint, values_from = n_cells, values_fill = 0)

message("\n=== DOMAIN GROUP CELL COUNTS ===")
print(group_summary)

write.csv(group_summary,
          file.path(output_dir, "PatientA_FunctionalGroup_Summary.csv"),
          row.names = FALSE)

# Clone by group
clone_by_group <- meta %>%
  filter(!is.na(group_name), !is.na(PatA_Clone)) %>%
  group_by(group_name, PatA_Clone) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group_name) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup()

message("\n=== CLONE DISTRIBUTION BY GROUP ===")
print(clone_by_group %>% filter(PatA_Clone == "PatA_Clone3") %>%
        select(group_name, n, pct) %>% arrange(desc(pct)))

write.csv(clone_by_group,
          file.path(output_dir, "PatientA_Clone_by_FunctionalGroup.csv"),
          row.names = FALSE)

message("\n=== SCRIPT COMPLETE ===")
message("Output directory: ", output_dir)
