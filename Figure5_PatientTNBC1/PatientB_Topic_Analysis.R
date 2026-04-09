# ==============================================================================
# Patient B: Topic Analysis
# Key Topics: 24 (Basal Stress), 12 (Mesenchymal/Basal), 28 (Immune Crosstalk),
#             1 (Proliferation), 2 (unknown)
# ==============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

output_dir <- "PatientB_Analysis/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load cancer cells with topics
cat("Loading Seurat object with topics...\n")
obj <- subset(readRDS("data/CosMx_SMMART_345k_clean.rds"), broad_lineage == "Cancer")

cat("Total cancer cells:", ncol(obj), "\n")
cat("Patients:", paste(unique(obj$Patient), collapse = ", "), "\n")
cat("Available metadata:", paste(head(colnames(obj@meta.data), 20), collapse = ", "), "\n")

# ==============================================================================
# Identify Topic Columns
# ==============================================================================

topic_cols <- grep("^Topic_", colnames(obj@meta.data), value = TRUE)
cat("\nTopic columns found:", length(topic_cols), "\n")
cat("Topics:", paste(head(topic_cols, 10), collapse = ", "), "...\n")

# Key topics of interest
key_topics <- c("Topic_24", "Topic_12", "Topic_1", "Topic_28", "Topic_2")
key_topics <- key_topics[key_topics %in% topic_cols]
cat("Key topics available:", paste(key_topics, collapse = ", "), "\n")

# Topic annotations (from user)
topic_annotations <- c(
  "Topic_24" = "Basal Stress Adaptation",
  "Topic_12" = "Mesenchymal/Basal",
  "Topic_28" = "Immune Crosstalk",
  "Topic_1" = "Proliferation",
  "Topic_2" = "Unknown"
)

# ==============================================================================
# Extract Topic Data
# ==============================================================================

cat("\nExtracting topic data...\n")

topic_df <- obj@meta.data %>%
  select(Patient, Timepoint, any_of(topic_cols)) %>%
  mutate(cell_id = rownames(obj@meta.data))

# Cell counts
cat("\nCell counts by Patient-Timepoint:\n")
print(table(topic_df$Patient, topic_df$Timepoint))

# ==============================================================================
# Calculate Mean Topic Scores Per Patient-Timepoint
# ==============================================================================

cat("\nCalculating mean topic scores...\n")

topic_means <- topic_df %>%
  group_by(Patient, Timepoint) %>%
  summarise(
    across(all_of(topic_cols), ~mean(.x, na.rm = TRUE)),
    N_cells = n(),
    .groups = "drop"
  )

# Patient order
patient_order <- c("Patient_B", "Patient_A", "Patient_C", "Patient_D")
topic_means$Patient <- factor(topic_means$Patient, levels = patient_order)

# ==============================================================================
# ANALYSIS 1: Key Topics - Patient B vs ER+ Comparison
# ==============================================================================

cat("\nAnalyzing key topics...\n")

# Get early and late for each patient
topic_early_late <- topic_means %>%
  group_by(Patient) %>%
  summarise(
    across(all_of(key_topics), list(
      Early = ~first(.x),
      Late = ~last(.x),
      FC = ~log2((last(.x) + 0.1) / (first(.x) + 0.1))
    )),
    .groups = "drop"
  )

# Print key topic summary
cat("\n=== KEY TOPICS SUMMARY ===\n")
for(topic in key_topics) {
  annotation <- topic_annotations[topic]
  cat("\n", topic, " (", annotation, "):\n", sep = "")

  for(pat in patient_order) {
    early_col <- paste0(topic, "_Early")
    late_col <- paste0(topic, "_Late")
    fc_col <- paste0(topic, "_FC")

    if(all(c(early_col, late_col, fc_col) %in% colnames(topic_early_late))) {
      pat_data <- topic_early_late %>% filter(Patient == pat)
      if(nrow(pat_data) > 0) {
        cat(sprintf("  %s: Early=%.3f, Late=%.3f, Log2FC=%.2f\n",
                    pat, pat_data[[early_col]], pat_data[[late_col]], pat_data[[fc_col]]))
      }
    }
  }
}

# ==============================================================================
# PLOT 1: Key Topics Barplot - All Patients
# ==============================================================================

cat("\nGenerating key topics barplot...\n")

# Prepare data for key topics
key_topic_long <- topic_means %>%
  select(Patient, Timepoint, all_of(key_topics)) %>%
  pivot_longer(
    cols = all_of(key_topics),
    names_to = "Topic",
    values_to = "Score"
  ) %>%
  mutate(
    Topic_Label = case_when(
      Topic == "Topic_24" ~ "Topic 24\n(Basal Stress)",
      Topic == "Topic_12" ~ "Topic 12\n(Mesench./Basal)",
      Topic == "Topic_28" ~ "Topic 28\n(Immune Crosstalk)",
      Topic == "Topic_1" ~ "Topic 1\n(Proliferation)",
      Topic == "Topic_2" ~ "Topic 2",
      TRUE ~ Topic
    ),
    Patient_Label = factor(Patient,
                           levels = c("Patient_B", "Patient_A", "Patient_C", "Patient_D"),
                           labels = c("B (TNBC)", "A (ER+)", "C (ER+)", "D (ER+)"))
  )

# Get early and late only
key_topic_summary <- key_topic_long %>%
  group_by(Patient, Patient_Label, Topic, Topic_Label) %>%
  summarise(
    Early = first(Score),
    Late = last(Score),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c("Early", "Late"),
    names_to = "Stage",
    values_to = "Score"
  )

# Patient colors
patient_colors <- c(
  "B (TNBC)" = "#377EB8",
  "A (ER+)" = "#E41A1C",
  "C (ER+)" = "#4DAF4A",
  "D (ER+)" = "#984EA3"
)

p1 <- ggplot(key_topic_summary, aes(x = Patient_Label, y = Score, fill = Patient_Label)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_grid(Topic_Label ~ Stage, scales = "free_y") +
  scale_fill_manual(values = patient_colors, guide = "none") +
  labs(
    title = "Key Topics Across Patients: TNBC vs ER+",
    subtitle = "Topics: Basal Stress (24), Mesenchymal (12), Immune (28), Proliferation (1)",
    x = NULL,
    y = "Mean Topic Score"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey40", size = 9),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "PatientB_KeyTopics_AllPatients.pdf"),
       p1, width = 10, height = 12)

# ==============================================================================
# PLOT 2: Topic Change Heatmap - Patient B Focus
# ==============================================================================

cat("Generating topic change heatmap...\n")

# Calculate fold changes for all topics
topic_fc <- topic_means %>%
  group_by(Patient) %>%
  summarise(
    across(all_of(topic_cols), ~log2((last(.x) + 0.1) / (first(.x) + 0.1))),
    .groups = "drop"
  )

# Pivot for heatmap
topic_fc_long <- topic_fc %>%
  pivot_longer(
    cols = all_of(topic_cols),
    names_to = "Topic",
    values_to = "Log2FC"
  ) %>%
  mutate(
    Topic_Num = as.numeric(sub("Topic_", "", Topic)),
    Patient_Label = factor(Patient,
                           levels = c("Patient_B", "Patient_A", "Patient_C", "Patient_D"),
                           labels = c("B (TNBC)", "A (ER+)", "C (ER+)", "D (ER+)"))
  )

# Highlight key topics
topic_fc_long$Is_Key <- topic_fc_long$Topic %in% key_topics

# Order topics by Patient B FC
topic_order <- topic_fc_long %>%
  filter(Patient == "Patient_B") %>%
  arrange(desc(Log2FC)) %>%
  pull(Topic)

topic_fc_long$Topic <- factor(topic_fc_long$Topic, levels = topic_order)

p2 <- ggplot(topic_fc_long, aes(x = Patient_Label, y = Topic, fill = Log2FC)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_point(data = filter(topic_fc_long, Is_Key),
             aes(x = Patient_Label, y = Topic),
             shape = 8, size = 2, color = "black") +
  scale_fill_gradient2(
    low = "#3575B5", mid = "white", high = "#D73027",
    midpoint = 0, name = "Log2FC",
    limits = c(-3, 3), oob = scales::squish
  ) +
  labs(
    title = "Topic Changes Across All Patients",
    subtitle = "Stars (*) mark key topics of interest | Topics ordered by Patient B FC",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "grey40", size = 9),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "PatientB_AllTopics_Heatmap.pdf"),
       p2, width = 7, height = 14)

# ==============================================================================
# ANALYSIS 2: Find Topics That Distinguish Patient B
# ==============================================================================

cat("\nFinding topics that distinguish Patient B...\n")

# Calculate divergence (Patient B - mean of ER+)
topic_divergence <- topic_fc_long %>%
  pivot_wider(names_from = Patient, values_from = Log2FC) %>%
  mutate(
    Mean_ERplus = (Patient_A + Patient_C + Patient_D) / 3,
    PatientB_Divergence = Patient_B - Mean_ERplus
  ) %>%
  arrange(desc(abs(PatientB_Divergence)))

cat("\nTop 10 Topics Most Different in Patient B vs ER+:\n")
print(head(topic_divergence %>% select(Topic, Patient_B, Mean_ERplus, PatientB_Divergence), 10))

cat("\nBottom 10 Topics Most Similar in Patient B vs ER+:\n")
print(tail(topic_divergence %>% select(Topic, Patient_B, Mean_ERplus, PatientB_Divergence), 10))

# ==============================================================================
# PLOT 3: Top Divergent Topics
# ==============================================================================

cat("Generating top divergent topics plot...\n")

# Get top 10 most divergent
top_divergent <- head(topic_divergence, 10)

top_div_long <- top_divergent %>%
  select(Topic, Patient_B, Mean_ERplus) %>%
  pivot_longer(
    cols = c("Patient_B", "Mean_ERplus"),
    names_to = "Group",
    values_to = "Log2FC"
  ) %>%
  mutate(
    Group = factor(Group,
                   levels = c("Patient_B", "Mean_ERplus"),
                   labels = c("TNBC (B)", "ER+ Mean")),
    Topic = factor(Topic, levels = top_divergent$Topic)
  )

p3 <- ggplot(top_div_long, aes(x = Topic, y = Log2FC, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("TNBC (B)" = "#377EB8", "ER+ Mean" = "#E41A1C"),
                    name = "Subtype") +
  labs(
    title = "Top 10 Most Divergent Topics: TNBC vs ER+",
    subtitle = "Topics showing largest difference between Patient B and ER+ patients",
    x = NULL,
    y = "Log2 Fold Change (Late vs Early)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey40", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "PatientB_TopDivergent_Topics.pdf"),
       p3, width = 10, height = 6)

# ==============================================================================
# PLOT 4: Key Topic Trajectories
# ==============================================================================

cat("Generating topic trajectories...\n")

# All timepoints for key topics
trajectory_data <- topic_means %>%
  select(Patient, Timepoint, all_of(key_topics)) %>%
  pivot_longer(
    cols = all_of(key_topics),
    names_to = "Topic",
    values_to = "Score"
  ) %>%
  mutate(
    Topic_Label = case_when(
      Topic == "Topic_24" ~ "Topic 24 (Basal Stress)",
      Topic == "Topic_12" ~ "Topic 12 (Mesenchymal)",
      Topic == "Topic_28" ~ "Topic 28 (Immune)",
      Topic == "Topic_1" ~ "Topic 1 (Proliferation)",
      Topic == "Topic_2" ~ "Topic 2",
      TRUE ~ Topic
    ),
    Patient_Label = factor(Patient,
                           levels = c("Patient_B", "Patient_A", "Patient_C", "Patient_D"),
                           labels = c("B (TNBC)", "A (ER+)", "C (ER+)", "D (ER+)"))
  ) %>%
  group_by(Patient) %>%
  mutate(TimeOrder = row_number()) %>%
  ungroup()

p4 <- ggplot(trajectory_data, aes(x = TimeOrder, y = Score,
                                   color = Patient_Label, group = Patient_Label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~Topic_Label, scales = "free_y", ncol = 3) +
  scale_color_manual(values = patient_colors, name = "Patient") +
  scale_x_continuous(breaks = 1:4, labels = c("1", "2", "3", "4")) +
  labs(
    title = "Key Topic Trajectories Over Time",
    subtitle = "TNBC (blue) vs ER+ patients (red/green/purple)",
    x = "Timepoint",
    y = "Mean Topic Score"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey40", size = 10),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "PatientB_KeyTopic_Trajectories.pdf"),
       p4, width = 12, height = 8)

# ==============================================================================
# ANALYSIS 3: Correlation Between Topics and L-R Patterns
# ==============================================================================

cat("\nAnalyzing topic-LR correlations...\n")

# Save topic divergence data
write.csv(topic_divergence, file.path(output_dir, "PatientB_Topic_Divergence.csv"),
          row.names = FALSE)

# Save all topic means
write.csv(topic_means, file.path(output_dir, "PatientB_Topic_Means_AllPatients.csv"),
          row.names = FALSE)

# ==============================================================================
# Print Summary
# ==============================================================================

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("PATIENT B TOPIC ANALYSIS SUMMARY\n")
cat("=", rep("=", 70), "\n", sep = "")

cat("\nKey Topics of Interest:\n")
for(topic in key_topics) {
  pat_b_fc <- topic_fc_long %>%
    filter(Topic == topic, Patient == "Patient_B") %>%
    pull(Log2FC)
  er_mean <- topic_divergence %>%
    filter(Topic == topic) %>%
    pull(Mean_ERplus)

  cat(sprintf("  %s (%s): B=%.2f, ER+=%.2f, Divergence=%.2f\n",
              topic, topic_annotations[topic], pat_b_fc, er_mean, pat_b_fc - er_mean))
}

cat("\nMost TNBC-Specific Topics (B >> ER+):\n")
top5_up <- head(topic_divergence %>% filter(PatientB_Divergence > 0), 5)
for(i in 1:nrow(top5_up)) {
  cat(sprintf("  %s: Divergence = +%.2f\n",
              top5_up$Topic[i], top5_up$PatientB_Divergence[i]))
}

cat("\nMost ER+-Specific Topics (B << ER+):\n")
top5_down <- head(topic_divergence %>% filter(PatientB_Divergence < 0) %>%
                    arrange(PatientB_Divergence), 5)
for(i in 1:nrow(top5_down)) {
  cat(sprintf("  %s: Divergence = %.2f\n",
              top5_down$Topic[i], top5_down$PatientB_Divergence[i]))
}

cat("\n\nOUTPUT FILES:\n")
cat("  1. PatientB_KeyTopics_AllPatients.pdf - Key topics comparison\n")
cat("  2. PatientB_AllTopics_Heatmap.pdf - All topics heatmap\n")
cat("  3. PatientB_TopDivergent_Topics.pdf - Most divergent topics\n")
cat("  4. PatientB_KeyTopic_Trajectories.pdf - Topic trajectories\n")
cat("  5. PatientB_Topic_Divergence.csv - Divergence data\n")
cat("  6. PatientB_Topic_Means_AllPatients.csv - All topic means\n")
cat("\n")
