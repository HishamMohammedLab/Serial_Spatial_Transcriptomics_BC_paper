# ==============================================================================
# Patient B Context: Luminal vs Basal Identity Comparison Across All Patients
# Shows TNBC vs ER+ divergence in epithelial identity
# ==============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

output_dir <- "/LR_Analysis/PatientB_Analysis/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading Seurat object...\n")
obj <- subset(readRDS("data/CosMx_SMMART_345k_clean.rds"), broad_lineage == "Cancer")

cat("Total cancer cells:", ncol(obj), "\n")
cat("Patients:", paste(unique(obj$Patient), collapse = ", "), "\n")

# ==============================================================================
# Define Key Marker Genes - Luminal vs Basal Identity
# ==============================================================================

luminal_genes <- c("ESR1", "GATA3", "FOXA1", "KRT8", "KRT18", "KRT19",
                   "XBP1", "ERBB3", "S100P", "PGR", "AR", "SPDEF")

basal_genes <- c("KRT5", "KRT14", "KRT15", "KRT16", "KRT17", "KRT23",
                 "VIM", "EGFR", "CDH3", "CAV1", "LAMB3", "LAMC2")

emt_genes <- c("VIM", "CD44", "MMP7", "SERPINA3", "TM4SF1", "SNAI2",
               "ZEB1", "TWIST1", "FN1", "CDH2")

stress_genes <- c("CRYAB", "HSPA1A", "HSP90AA1", "HSPB1", "DNAJB1",
                  "BAG3", "HMGB2", "ATF4", "DDIT3")

# Check availability
all_genes <- unique(c(luminal_genes, basal_genes, emt_genes, stress_genes))
available <- all_genes[all_genes %in% rownames(obj)]
missing <- all_genes[!all_genes %in% rownames(obj)]

cat("\nAvailable genes:", length(available), "/", length(all_genes), "\n")
if(length(missing) > 0) {
  cat("Missing:", paste(missing, collapse = ", "), "\n")
}

# Update gene lists to available only
luminal_genes <- luminal_genes[luminal_genes %in% available]
basal_genes <- basal_genes[basal_genes %in% available]
emt_genes <- emt_genes[emt_genes %in% available]
stress_genes <- stress_genes[stress_genes %in% available]

# ==============================================================================
# Extract Expression Data
# ==============================================================================

cat("\nExtracting expression data...\n")
expr_matrix <- GetAssayData(obj, layer = "data")

# Create data frame
genes_to_use <- unique(c(luminal_genes, basal_genes, emt_genes, stress_genes))
expr_df <- as.data.frame(t(as.matrix(expr_matrix[genes_to_use, , drop = FALSE])))
expr_df$Patient <- obj$Patient
expr_df$Timepoint <- obj$Timepoint
expr_df$cell_id <- rownames(expr_df)

# ==============================================================================
# Calculate Mean Expression Per Patient-Timepoint
# ==============================================================================

cat("Calculating mean expression per patient-timepoint...\n")

mean_expr <- expr_df %>%
  group_by(Patient, Timepoint) %>%
  summarise(
    across(all_of(genes_to_use), ~mean(.x, na.rm = TRUE)),
    N_cells = n(),
    .groups = "drop"
  )

# Create patient order (B first, then A, C, D)
patient_order <- c("Patient_B", "Patient_A", "Patient_C", "Patient_D")
mean_expr$Patient <- factor(mean_expr$Patient, levels = patient_order)

# ==============================================================================
# Calculate Module Scores
# ==============================================================================

# Helper function to calculate module score (mean of z-scored genes)
calc_module_score <- function(data, genes) {
  if(length(genes) == 0) return(rep(NA, nrow(data)))
  gene_cols <- genes[genes %in% colnames(data)]
  if(length(gene_cols) == 0) return(rep(NA, nrow(data)))
  rowMeans(data[, gene_cols, drop = FALSE], na.rm = TRUE)
}

mean_expr$Luminal_Score <- calc_module_score(mean_expr, luminal_genes)
mean_expr$Basal_Score <- calc_module_score(mean_expr, basal_genes)
mean_expr$EMT_Score <- calc_module_score(mean_expr, emt_genes)
mean_expr$Stress_Score <- calc_module_score(mean_expr, stress_genes)

# ==============================================================================
# PLOT 1: Luminal vs Basal Barplot (All Patients)
# ==============================================================================

cat("\nGenerating Luminal vs Basal comparison...\n")

# Prepare data for plotting
module_data <- mean_expr %>%
  select(Patient, Timepoint, N_cells, Luminal_Score, Basal_Score, EMT_Score, Stress_Score) %>%
  pivot_longer(
    cols = ends_with("_Score"),
    names_to = "Module",
    values_to = "Score"
  ) %>%
  mutate(
    Module = gsub("_Score", "", Module),
    Module = factor(Module, levels = c("Luminal", "Basal", "EMT", "Stress")),
    Patient_Label = factor(Patient,
                           levels = c("Patient_B", "Patient_A", "Patient_C", "Patient_D"),
                           labels = c("B (TNBC)", "A (ER+)", "C (ER+)", "D (ER+)"))
  )

# Take early and late timepoints for each patient
module_summary <- module_data %>%
  group_by(Patient, Patient_Label, Module) %>%
  summarise(
    Early = first(Score),
    Late = last(Score),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c("Early", "Late"),
    names_to = "Stage",
    values_to = "Score"
  ) %>%
  mutate(Stage = factor(Stage, levels = c("Early", "Late")))

# Patient colors
patient_colors <- c(
  "B (TNBC)" = "#377EB8",  # Blue for TNBC
  "A (ER+)" = "#E41A1C",   # Red for ER+
  "C (ER+)" = "#4DAF4A",   # Green for ER+
  "D (ER+)" = "#984EA3"    # Purple for ER+
)

p1 <- ggplot(module_summary, aes(x = Patient_Label, y = Score, fill = Patient_Label)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_grid(Module ~ Stage, scales = "free_y") +
  scale_fill_manual(values = patient_colors, guide = "none") +
  labs(
    title = "Epithelial Identity Programs Across Patients",
    subtitle = "TNBC (Patient B) vs ER+ (Patients A, C, D)",
    x = NULL,
    y = "Mean Expression Score"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 10),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    panel.grid = element_blank()
  )

ggsave(file.path(output_dir, "PatientB_LuminalBasal_AllPatients.pdf"),
       p1, width = 10, height = 10)

# ==============================================================================
# PLOT 2: Key Individual Genes Comparison
# ==============================================================================

cat("Generating key gene comparison...\n")

# Select key genes for visualization
key_luminal <- c("ESR1", "GATA3", "KRT8", "KRT19")
key_basal <- c("KRT5", "KRT17", "VIM", "EGFR")

key_genes <- c(key_luminal[key_luminal %in% genes_to_use],
               key_basal[key_basal %in% genes_to_use])

key_gene_data <- mean_expr %>%
  select(Patient, Timepoint, all_of(key_genes)) %>%
  pivot_longer(
    cols = all_of(key_genes),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  mutate(
    Category = case_when(
      Gene %in% key_luminal ~ "Luminal",
      Gene %in% key_basal ~ "Basal",
      TRUE ~ "Other"
    ),
    Patient_Label = factor(Patient,
                           levels = c("Patient_B", "Patient_A", "Patient_C", "Patient_D"),
                           labels = c("B (TNBC)", "A (ER+)", "C (ER+)", "D (ER+)"))
  )

# Get early and late
key_gene_summary <- key_gene_data %>%
  group_by(Patient, Patient_Label, Gene, Category) %>%
  summarise(
    Early = first(Expression),
    Late = last(Expression),
    .groups = "drop"
  ) %>%
  mutate(FC = log2((Late + 0.01) / (Early + 0.01)))

# Order genes
key_gene_summary$Gene <- factor(key_gene_summary$Gene, levels = key_genes)

p2 <- ggplot(key_gene_summary, aes(x = Gene, y = Late, fill = Patient_Label)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~Category, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = patient_colors, name = "Patient") +
  labs(
    title = "Key Luminal and Basal Markers at Late Timepoint",
    subtitle = "Expression comparison: TNBC vs ER+ patients",
    x = NULL,
    y = "Mean Normalized Expression"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 10),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "right",
    panel.grid = element_blank()
  )

ggsave(file.path(output_dir, "PatientB_KeyGenes_LuminalBasal.pdf"),
       p2, width = 11, height = 5)

# ==============================================================================
# PLOT 3: Luminal/Basal Ratio (Key Insight)
# ==============================================================================

cat("Generating Luminal/Basal ratio plot...\n")

ratio_data <- mean_expr %>%
  mutate(
    LB_Ratio = (Luminal_Score + 0.1) / (Basal_Score + 0.1),
    Log2_LB_Ratio = log2(LB_Ratio),
    Patient_Label = factor(Patient,
                           levels = c("Patient_B", "Patient_A", "Patient_C", "Patient_D"),
                           labels = c("B (TNBC)", "A (ER+)", "C (ER+)", "D (ER+)"))
  )

# Create x-axis positions
ratio_data <- ratio_data %>%
  group_by(Patient) %>%
  mutate(TimeOrder = row_number()) %>%
  ungroup()

p3 <- ggplot(ratio_data, aes(x = TimeOrder, y = Log2_LB_Ratio,
                              color = Patient_Label, group = Patient_Label)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = patient_colors, name = "Patient") +
  scale_x_continuous(breaks = 1:4, labels = c("Early", "Mid", "Late", "Final")) +
  labs(
    title = "Luminal/Basal Identity Ratio Over Time",
    subtitle = "Log2(Luminal Score / Basal Score) | Above 0 = Luminal-dominant",
    x = "Timepoint",
    y = "Log2(Luminal/Basal Ratio)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 10),
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  annotate("text", x = 0.7, y = 0.5, label = "Luminal\ndominant",
           hjust = 0, size = 3, color = "grey50") +
  annotate("text", x = 0.7, y = -0.5, label = "Basal\ndominant",
           hjust = 0, size = 3, color = "grey50")

ggsave(file.path(output_dir, "PatientB_LuminalBasal_Ratio.pdf"),
       p3, width = 8, height = 5)

# ==============================================================================
# PLOT 4: Heatmap of All Key Genes
# ==============================================================================

cat("Generating heatmap of key genes...\n")

# Prepare data for heatmap
heatmap_genes <- c(luminal_genes, basal_genes)
heatmap_genes <- heatmap_genes[heatmap_genes %in% genes_to_use]

heatmap_data <- mean_expr %>%
  select(Patient, Timepoint, all_of(heatmap_genes)) %>%
  unite("Sample", Patient, Timepoint, sep = "_") %>%
  column_to_rownames("Sample")

# Z-score normalize
heatmap_scaled <- scale(heatmap_data)

# Convert to long format for ggplot
heatmap_long <- as.data.frame(heatmap_scaled) %>%
  rownames_to_column("Sample") %>%
  pivot_longer(
    cols = -Sample,
    names_to = "Gene",
    values_to = "Z_score"
  ) %>%
  mutate(
    Patient = sub("_Bx.*", "", Sample),
    Timepoint = sub(".*_(Bx[0-9]+)", "\\1", Sample),
    Category = case_when(
      Gene %in% luminal_genes ~ "Luminal",
      Gene %in% basal_genes ~ "Basal",
      TRUE ~ "Other"
    ),
    Patient_Label = factor(Patient,
                           levels = c("Patient_B", "Patient_A", "Patient_C", "Patient_D"),
                           labels = c("B", "A", "C", "D")),
    Gene = factor(Gene, levels = rev(heatmap_genes))
  )

p4 <- ggplot(heatmap_long, aes(x = Sample, y = Gene, fill = Z_score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#3575B5", mid = "white", high = "#D73027",
    midpoint = 0, name = "Z-score",
    limits = c(-2.5, 2.5), oob = scales::squish
  ) +
  facet_grid(Category ~ Patient_Label, scales = "free", space = "free") +
  labs(
    title = "Luminal and Basal Gene Expression Across Patients",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "PatientB_LuminalBasal_Heatmap.pdf"),
       p4, width = 12, height = 10)

# ==============================================================================
# Save Summary Statistics
# ==============================================================================

cat("\nSaving summary statistics...\n")

summary_stats <- mean_expr %>%
  select(Patient, Timepoint, N_cells, Luminal_Score, Basal_Score, EMT_Score, Stress_Score) %>%
  mutate(
    LB_Ratio = (Luminal_Score + 0.1) / (Basal_Score + 0.1),
    Log2_LB_Ratio = log2(LB_Ratio)
  ) %>%
  arrange(Patient, Timepoint)

write.csv(summary_stats, file.path(output_dir, "PatientB_LuminalBasal_Summary.csv"),
          row.names = FALSE)

# ==============================================================================
# Print Summary
# ==============================================================================

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("LUMINAL VS BASAL IDENTITY SUMMARY\n")
cat("=", rep("=", 70), "\n", sep = "")

for(pat in patient_order) {
  pat_data <- summary_stats %>% filter(Patient == pat)
  subtype <- ifelse(pat == "Patient_B", "TNBC", "ER+")

  cat("\n", pat, " (", subtype, "):\n", sep = "")
  cat(rep("-", 50), "\n", sep = "")

  for(i in 1:nrow(pat_data)) {
    cat(sprintf("  %s: Luminal=%.3f, Basal=%.3f, L/B ratio=%.2f\n",
                pat_data$Timepoint[i],
                pat_data$Luminal_Score[i],
                pat_data$Basal_Score[i],
                pat_data$LB_Ratio[i]))
  }
}

cat("\n\n")
cat("OUTPUT FILES:\n")
cat("  1. PatientB_LuminalBasal_AllPatients.pdf - Module score comparison\n")
cat("  2. PatientB_KeyGenes_LuminalBasal.pdf - Key individual genes\n")
cat("  3. PatientB_LuminalBasal_Ratio.pdf - L/B ratio over time\n")
cat("  4. PatientB_LuminalBasal_Heatmap.pdf - Full gene heatmap\n")
cat("  5. PatientB_LuminalBasal_Summary.csv - Summary statistics\n")
cat("\n")
