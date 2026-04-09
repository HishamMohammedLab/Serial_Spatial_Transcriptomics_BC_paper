# ==============================================================================
# TOPICS 7, 13, 32 PROJECTION ON UMAP
# ==============================================================================
# Same UMAP as publication script, colored by each topic
# Grey for low values, red gradient for signal
# ==============================================================================

library(Seurat)
library(ggplot2)
library(uwot)
library(scales)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

OBJ_PATH <- "data/CosMx_SMMART_345k_clean.rds"
OUTPUT_DIR <- "Topic_UMAP_Publication"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

TOPIC_COLS <- paste0("Topic_", 1:30)
TARGET_TOPICS <- c("Topic_27", "Topic_24")

# UMAP parameters 
UMAP_NEIGHBORS <- 30
UMAP_MIN_DIST <- 0.3
UMAP_SEED <- 42

# ==============================================================================
# LOAD DATA
# ==============================================================================

message("Loading Seurat object...")
obj_all <- readRDS(OBJ_PATH)

# Identify cells with valid topic scores (same as publication script)
existing_topics <- TOPIC_COLS[TOPIC_COLS %in% colnames(obj_all@meta.data)]
topic_data <- obj_all@meta.data[, existing_topics]
valid_cells <- rowSums(is.na(topic_data)) == 0
nonzero_cells <- rowSums(topic_data, na.rm = TRUE) > 0
cells_with_topics <- valid_cells & nonzero_cells
cells_use <- rownames(obj_all@meta.data)[cells_with_topics]

message(sprintf("Cells with valid topic scores: %d", length(cells_use)))

# Subset
topic_matrix <- as.matrix(topic_data[cells_use, ])

# ==============================================================================
# RUN UMAP (same as publication script)
# ==============================================================================

message("Running UMAP...")
set.seed(UMAP_SEED)
umap_result <- uwot::umap(
  topic_matrix,
  n_neighbors = UMAP_NEIGHBORS,
  min_dist = UMAP_MIN_DIST,
  metric = "cosine",
  n_components = 2,
  verbose = TRUE
)

message("UMAP complete.\n")

# ==============================================================================
# CREATE PLOTS FOR EACH TOPIC
# ==============================================================================

for (TARGET_TOPIC in TARGET_TOPICS) {
  
  message(sprintf("Processing %s...", TARGET_TOPIC))
  
  # Create data frame
  umap_df <- data.frame(
    UMAP_1 = umap_result[, 1],
    UMAP_2 = umap_result[, 2],
    topic_value = topic_matrix[, TARGET_TOPIC]
  )
  
  # Order by topic value so high values plot on top
  umap_df <- umap_df[order(umap_df$topic_value), ]
  
  # Set up color scale
  min_val <- min(umap_df$topic_value)
  max_val <- max(umap_df$topic_value)
  
  # Set threshold at 25th percentile - below this is "background" grey
  threshold <- quantile(umap_df$topic_value, 0.25)
  
  message(sprintf("  Range: %.4f to %.4f", min_val, max_val))
  message(sprintf("  Threshold (25th percentile): %.4f", threshold))
  
  # Clamp below-threshold values for display
  umap_df$topic_display <- umap_df$topic_value
  umap_df$topic_display[umap_df$topic_display < threshold] <- threshold
  
  # Create plot
  p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = topic_display)) +
    geom_point(size = 0.3, alpha = 0.6) +
    scale_color_gradient(
      low = "grey90",
      high = "#BD0026",
      limits = c(threshold, max_val),
      name = "Weight"
    ) +
    labs(
      x = "Topic UMAP 1",
      y = "Topic UMAP 2",
      title = sprintf("Topic UMAP - %s", TARGET_TOPIC)
    ) +
    coord_equal() +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 9),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
  # Save
  filename <- sprintf("UMAP_%s.pdf", TARGET_TOPIC)
  ggsave(
    file.path(OUTPUT_DIR, filename),
    p,
    width = 8, height = 7
  )
  
  message(sprintf("  Saved: %s\n", filename))
}

message("Done!")
