# ==============================================================================
# SMMART: Topic-Based UMAP - Patient & Biopsy Colors
# Each patient has a base color, biopsies are shades of that color
# ==============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(uwot)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

TOPICS_OBJ_PATH <- "data/CosMx_SMMART_345k_clean.rds"  

OUTPUT_DIR <- "Topic_UMAP_Publication/LR_Analysis/Epithelial_Analysis/"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Topic columns
TOPIC_COLS <- paste0("Topic_", 1:30)

# UMAP parameters
UMAP_NEIGHBORS <- 30
UMAP_MIN_DIST <- 0.3
UMAP_SEED <- 42

# ==============================================================================
# COLOR PALETTE - Patient/Biopsy Shades
# ==============================================================================

# Patient A = Red shades
# Patient B = Green shades
# Patient C = Blue shades
# Patient D = Yellow/Gold shades

patient_biopsy_colors <- c(
  # Patient A - Red shades (Bx1, Bx2, Bx3, Bx4)
  "Patient_A_Bx1" = "#FF6B6B",  # Light red
  "Patient_A_Bx2" = "#E63946",  # Medium red
  "Patient_A_Bx3" = "#C1121F",  # Dark red
  "Patient_A_Bx4" = "#780000",  # Darkest red
  
  # Patient B - Green shades
  "Patient_B_Bx1" = "#95D5B2",  # Light green
  "Patient_B_Bx2" = "#52B788",  # Medium green
  "Patient_B_Bx3" = "#2D6A4F",  # Dark green
  "Patient_B_Bx4" = "#1B4332",  # Darkest green
  
  # Patient C - Blue shades
  "Patient_C_Bx1" = "#90E0EF",  # Light blue
  "Patient_C_Bx2" = "#00B4D8",  # Medium blue
  "Patient_C_Bx3" = "#0077B6",  # Dark blue
  "Patient_C_Bx4" = "#03045E",  # Darkest blue
  
  # Patient D - Yellow/Gold shades
  "Patient_D_Bx1" = "#FFE566",  # Light yellow
  "Patient_D_Bx2" = "#FFD60A",  # Medium yellow
  "Patient_D_Bx3" = "#E6AC00",  # Gold
  "Patient_D_Bx4" = "#B38600"   # Dark gold
)

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
# LOAD DATA
# ==============================================================================

message("=== Loading topic object ===")
obj <- readRDS(TOPICS_OBJ_PATH)
obj <- UpdateSeuratObject(obj)
message(paste("Topic object cells:", ncol(obj)))

# ==============================================================================
# CREATE PATIENT_BIOPSY COLUMN
# ==============================================================================

message("\n=== Creating Patient_Biopsy column ===")

# Check existing columns
message("Patient values:")
print(table(obj@meta.data$Patient, useNA = "ifany"))

message("\nTimepoint values:")
print(table(obj@meta.data$Timepoint, useNA = "ifany"))

# Create combined Patient_Biopsy column
obj@meta.data$Patient_Biopsy <- paste0(obj@meta.data$Patient, "_", obj@meta.data$Timepoint)

message("\nPatient_Biopsy distribution:")
print(table(obj@meta.data$Patient_Biopsy, useNA = "ifany"))

# Set factor levels in order (Patient A first, then B, C, D; Bx1-4 within each)
all_combinations <- c(
  "Patient_A_Bx1", "Patient_A_Bx2", "Patient_A_Bx3", "Patient_A_Bx4",
  "Patient_B_Bx1", "Patient_B_Bx2", "Patient_B_Bx3", "Patient_B_Bx4",
  "Patient_C_Bx1", "Patient_C_Bx2", "Patient_C_Bx3", "Patient_C_Bx4",
  "Patient_D_Bx1", "Patient_D_Bx2", "Patient_D_Bx3", "Patient_D_Bx4"
)

# Only keep levels that exist in the data
existing_levels <- all_combinations[all_combinations %in% unique(obj@meta.data$Patient_Biopsy)]
obj@meta.data$Patient_Biopsy <- factor(obj@meta.data$Patient_Biopsy, levels = existing_levels)

message("\nFinal Patient_Biopsy levels:")
print(levels(obj@meta.data$Patient_Biopsy))

# ==============================================================================
# PREPARE TOPIC MATRIX & RUN UMAP
# ==============================================================================

message("\n=== Preparing topic matrix ===")

existing_topics <- TOPIC_COLS[TOPIC_COLS %in% colnames(obj@meta.data)]
topic_data <- obj@meta.data[, existing_topics]

valid_cells <- rowSums(is.na(topic_data)) == 0
nonzero_cells <- rowSums(topic_data, na.rm = TRUE) != 0
cells_with_topics <- valid_cells & nonzero_cells
cells_use <- rownames(obj@meta.data)[cells_with_topics]

message(paste("Cells with valid topic scores:", length(cells_use)))

topic_matrix <- as.matrix(topic_data[cells_use, ])

message("\n=== Running UMAP ===")
set.seed(UMAP_SEED)
umap_result <- uwot::umap(
  topic_matrix,
  n_neighbors = UMAP_NEIGHBORS,
  min_dist = UMAP_MIN_DIST,
  metric = "cosine",
  n_components = 2,
  verbose = TRUE
)
message("UMAP complete.")

# ==============================================================================
# CREATE PLOTTING DATA FRAME
# ==============================================================================

obj_topics <- subset(obj, cells = cells_use)

umap_df <- data.frame(
  UMAP_1 = umap_result[, 1],
  UMAP_2 = umap_result[, 2],
  Patient = obj_topics@meta.data$Patient,
  Timepoint = obj_topics@meta.data$Timepoint,
  Patient_Biopsy = obj_topics@meta.data$Patient_Biopsy,
  Domain = obj_topics@meta.data$spatial_domain
)

message("\n=== Patient_Biopsy in plotting dataframe ===")
print(table(umap_df$Patient_Biopsy, useNA = "ifany"))

# ==============================================================================
# PLOT: UMAP BY PATIENT & BIOPSY
# ==============================================================================

message("\n=== Generating Patient/Biopsy UMAP ===")

# Filter colors to only those present in data
colors_to_use <- patient_biopsy_colors[names(patient_biopsy_colors) %in% levels(umap_df$Patient_Biopsy)]

message("Colors being used:")
print(colors_to_use)

p_patient_biopsy <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = Patient_Biopsy)) +
  geom_point(size = 0.3, alpha = 0.6) +
  scale_color_manual(values = colors_to_use, na.value = "grey80") +
  labs(
    x = "Topic UMAP 1",
    y = "Topic UMAP 2",
    title = paste0("Epithelial Cells by Patient & Biopsy (n = ", nrow(umap_df), ")"),
    color = "Patient / Biopsy"
  ) +
  coord_equal() +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggsave(
  file.path(OUTPUT_DIR, "TopicUMAP_Patient_Biopsy.pdf"),
  p_patient_biopsy,
  width = 11, height = 8
)
message(paste("Saved:", file.path(OUTPUT_DIR, "TopicUMAP_Patient_Biopsy.pdf")))

p_patient_domain <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = Domain)) +
  geom_point(size = 0.3, alpha = 0.6) +
  scale_color_manual(values = DOMAIN_COLORS, na.value = "grey80") +
  labs(
    x = "Topic UMAP 1",
    y = "Topic UMAP 2",
    title = paste0("Epithelial Cells by Spatial Domain (n = ", nrow(umap_df), ")"),
    color = "Domain"
  ) +
  coord_equal() +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggsave(
  file.path(OUTPUT_DIR, "TopicUMAP_Patient_Domain.pdf"),
  p_patient_domain,
  width = 11, height = 8
)
message(paste("Saved:", file.path(OUTPUT_DIR, "TopicUMAP_Patient_Domain.pdf")))

# ==============================================================================
# ADDITIONAL: Split by Patient
# ==============================================================================

message("\n=== Generating per-patient plots ===")

pdf(file.path(OUTPUT_DIR, "TopicUMAP_Patient_Biopsy_Split.pdf"), width = 12, height = 10)

# All patients together but faceted
p_facet <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = Patient_Biopsy)) +
  geom_point(size = 0.3, alpha = 0.6) +
  scale_color_manual(values = colors_to_use, na.value = "grey80") +
  facet_wrap(~Patient, ncol = 2) +
  labs(
    x = "Topic UMAP 1",
    y = "Topic UMAP 2",
    title = "Epithelial Cells by Patient & Biopsy",
    color = "Patient / Biopsy"
  ) +
  coord_equal() +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

print(p_facet)

dev.off()

message(paste("Saved:", file.path(OUTPUT_DIR, "TopicUMAP_Patient_Biopsy_Split.pdf")))

# ==============================================================================
# SUMMARY
# ==============================================================================

message("\n")
message("==============================================================================")
message("OUTPUTS")
message("==============================================================================")
message(paste("1.", file.path(OUTPUT_DIR, "TopicUMAP_Patient_Biopsy.pdf")))
message(paste("2.", file.path(OUTPUT_DIR, "TopicUMAP_Patient_Biopsy_Split.pdf")))

message("\nColor scheme:")
message("  Patient A (Red):    Bx1=#FF6B6B, Bx2=#E63946, Bx3=#C1121F, Bx4=#780000")
message("  Patient B (Green):  Bx1=#95D5B2, Bx2=#52B788, Bx3=#2D6A4F, Bx4=#1B4332")
message("  Patient C (Blue):   Bx1=#90E0EF, Bx2=#00B4D8, Bx3=#0077B6, Bx4=#03045E")
message("  Patient D (Yellow): Bx1=#FFE566, Bx2=#FFD60A, Bx3=#E6AC00, Bx4=#B38600")

message("\n=== COMPLETE ===")
