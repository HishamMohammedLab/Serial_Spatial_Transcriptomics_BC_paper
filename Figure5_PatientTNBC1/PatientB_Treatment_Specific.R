# ==============================================================================
# Patient B: Treatment-Specific Gene Analysis
#
# IMPORTANT CONTEXT:
# - Bx1: Lymph node
# - Bx4: Soft tissue (DIFFERENT SITE - major confound!)
#
# Treatments between biopsies:
# - Olaparib (PARP inhibitor) - long duration
# - Durvalumab (anti-PD-L1) - long duration
# - Sacituzumab govitecan (TROP2-ADC) - short, just before Bx4
#
# ==============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

output_dir <- "LR_Analysis/PatientB_Analysis/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading cancer cells with topics...\n")
obj <- subset(readRDS("data/CosMx_SMMART_345k_clean.rds"), broad_lineage == "Cancer")

# Subset to Patient B
patB <- subset(obj, subset = Patient == "Patient_B")
cat("Patient B cancer cells:", ncol(patB), "\n")
cat("By timepoint:\n")
print(table(patB$Timepoint))

# ==============================================================================
# Treatment-Specific Gene Signatures
# ==============================================================================

# These are cancer cell-intrinsic markers that should be LESS affected by site

treatment_genes <- list(

  # ----- OLAPARIB (PARP inhibitor) Response -----
  # PARP inhibitors cause DNA damage accumulation
  # Resistance involves HR restoration, drug efflux, replication fork protection

  "PARP_Target" = c("PARP1", "PARP2"),

  "DNA_Damage_Response" = c("ATM", "ATR", "CHEK1", "CHEK2", "TP53",
                             "BRCA1", "BRCA2", "RAD51", "PALB2"),

  "Homologous_Recombination" = c("RAD51", "BRCA1", "BRCA2", "PALB2",
                                  "FANCD2", "FANCA"),

  "PARP_Resistance" = c("ABCB1", "ABCG2",  # Drug efflux
                        "ATR", "CHEK1",     # Replication fork protection
                        "TP53BP1"),         # 53BP1 loss = HR restoration

  # ----- DURVALUMAB (anti-PD-L1) Response -----
  # PD-L1 blockade - look for immune evasion adaptations

  "PD_L1_Axis" = c("CD274", "PDCD1LG2"),  # PD-L1, PD-L2

  "Interferon_Response" = c("IFNGR1", "IFNGR2", "STAT1", "IRF1", "IRF9",
                            "ISG15", "MX1", "OAS1"),

  "Antigen_Presentation" = c("B2M", "HLA-A", "HLA-B", "HLA-C",
                              "TAP1", "TAP2", "TAPBP"),

  "Alternative_Checkpoints" = c("LGALS9", "HAVCR2",  # Galectin-9/TIM-3
                                 "CD276",             # B7-H3
                                 "VTCN1",             # B7-H4
                                 "TNFRSF14"),         # HVEM

  # ----- SACITUZUMAB (TROP2-ADC) Response -----
  # Targets TROP2, delivers SN-38 (topoisomerase I inhibitor)

  "TROP2_Target" = c("TACSTD2"),  # TROP2

  "Topoisomerase" = c("TOP1", "TOP2A", "TOP2B"),

  "Drug_Efflux" = c("ABCB1", "ABCG2", "ABCC1", "ABCC2"),

  # ----- General Stress/Resistance -----

  "Stress_Response" = c("HSP90AA1", "HSPA1A", "HSPB1", "CRYAB",
                        "ATF4", "DDIT3", "XBP1"),

  "Apoptosis_Evasion" = c("BCL2", "MCL1", "BIRC5", "XIAP"),

  # ----- Proliferation (treatment-independent control) -----

  "Proliferation" = c("MKI67", "PCNA", "CCND1", "CDK4", "CDK6")
)

# Check availability
all_genes <- unique(unlist(treatment_genes))
available <- all_genes[all_genes %in% rownames(patB)]
missing <- setdiff(all_genes, available)

cat("\nAvailable treatment-related genes:", length(available), "/", length(all_genes), "\n")
if(length(missing) > 0) {
  cat("Missing:", paste(missing, collapse = ", "), "\n")
}

# Update lists
treatment_genes <- lapply(treatment_genes, function(x) x[x %in% available])
treatment_genes <- treatment_genes[sapply(treatment_genes, length) > 0]

cat("\nGene categories available:\n")
for(cat in names(treatment_genes)) {
  cat(sprintf("  %s: %d genes\n", cat, length(treatment_genes[[cat]])))
}

# ==============================================================================
# Extract Expression
# ==============================================================================

cat("\nExtracting expression data...\n")

expr_mat <- GetAssayData(patB, layer = "data")
genes_use <- unique(unlist(treatment_genes))

expr_df <- as.data.frame(t(as.matrix(expr_mat[genes_use, , drop = FALSE])))
expr_df$Timepoint <- patB$Timepoint
expr_df$cell_id <- rownames(expr_df)

# Cell counts
n_bx1 <- sum(expr_df$Timepoint == "Bx1")
n_bx4 <- sum(expr_df$Timepoint == "Bx4")
cat("Cells: Bx1 (Lymph Node) =", n_bx1, ", Bx4 (Soft Tissue) =", n_bx4, "\n")

# ==============================================================================
# Calculate Statistics
# ==============================================================================

cat("\nCalculating statistics...\n")

stats_list <- list()

for(gene in genes_use) {
  bx1_vals <- expr_df[expr_df$Timepoint == "Bx1", gene]
  bx4_vals <- expr_df[expr_df$Timepoint == "Bx4", gene]

  mean_bx1 <- mean(bx1_vals, na.rm = TRUE)
  mean_bx4 <- mean(bx4_vals, na.rm = TRUE)

  pct_bx1 <- mean(bx1_vals > 0, na.rm = TRUE) * 100
  pct_bx4 <- mean(bx4_vals > 0, na.rm = TRUE) * 100

  log2fc <- log2((mean_bx4 + 0.01) / (mean_bx1 + 0.01))

  # Wilcoxon test
  wilcox_p <- tryCatch(
    wilcox.test(bx1_vals, bx4_vals)$p.value,
    error = function(e) NA
  )

  # Find category
  category <- NA
  for(cat_name in names(treatment_genes)) {
    if(gene %in% treatment_genes[[cat_name]]) {
      category <- cat_name
      break
    }
  }

  stats_list[[gene]] <- data.frame(
    Gene = gene,
    Category = category,
    Mean_Bx1 = mean_bx1,
    Mean_Bx4 = mean_bx4,
    Pct_Bx1 = pct_bx1,
    Pct_Bx4 = pct_bx4,
    Log2FC = log2fc,
    P_value = wilcox_p,
    stringsAsFactors = FALSE
  )
}

stats_df <- bind_rows(stats_list)
stats_df$P_adj <- p.adjust(stats_df$P_value, method = "BH")
stats_df$Sig <- ifelse(stats_df$P_adj < 0.05, "*", "")

# ==============================================================================
# PLOT 1: Treatment-Specific Gene Changes
# ==============================================================================

cat("\nGenerating treatment-specific plots...\n")

# Group categories by treatment
olaparib_cats <- c("PARP_Target", "DNA_Damage_Response", "Homologous_Recombination", "PARP_Resistance")
durvalumab_cats <- c("PD_L1_Axis", "Interferon_Response", "Antigen_Presentation", "Alternative_Checkpoints")
sacituzumab_cats <- c("TROP2_Target", "Topoisomerase", "Drug_Efflux")

# Plot function
plot_treatment_genes <- function(data, categories, title, subtitle) {
  plot_data <- data %>%
    filter(Category %in% categories) %>%
    mutate(
      Category = factor(Category, levels = categories),
      Gene = reorder(Gene, Log2FC)
    )

  if(nrow(plot_data) == 0) return(NULL)

  ggplot(plot_data, aes(x = Log2FC, y = Gene, fill = Log2FC > 0)) +
    geom_bar(stat = "identity") +
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +
    geom_text(aes(label = Sig, x = ifelse(Log2FC > 0, Log2FC + 0.1, Log2FC - 0.1)),
              hjust = ifelse(plot_data$Log2FC > 0, 0, 1), size = 4) +
    facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = c("TRUE" = "#D73027", "FALSE" = "#3575B5"), guide = "none") +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Log2 Fold Change (Bx4 vs Bx1)",
      y = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(color = "grey40", size = 9),
      strip.text.y = element_text(angle = 0, face = "bold", hjust = 0, size = 9),
      axis.text.y = element_text(size = 9),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# Olaparib response
p1 <- plot_treatment_genes(stats_df, olaparib_cats,
                            "Olaparib (PARP Inhibitor) Response Genes",
                            "DNA damage response, HR, resistance markers | * = adj.p < 0.05")

if(!is.null(p1)) {
  ggsave(file.path(output_dir, "PatientB_Olaparib_Response.pdf"),
         p1, width = 8, height = 8)
}

# Durvalumab response
p2 <- plot_treatment_genes(stats_df, durvalumab_cats,
                            "Durvalumab (Anti-PD-L1) Response Genes",
                            "PD-L1/immune checkpoint, IFN response, antigen presentation | * = adj.p < 0.05")

if(!is.null(p2)) {
  ggsave(file.path(output_dir, "PatientB_Durvalumab_Response.pdf"),
         p2, width = 8, height = 10)
}

# Sacituzumab response
p3 <- plot_treatment_genes(stats_df, sacituzumab_cats,
                            "Sacituzumab (TROP2-ADC) Response Genes",
                            "TROP2 target, topoisomerase, drug efflux | * = adj.p < 0.05")

if(!is.null(p3)) {
  ggsave(file.path(output_dir, "PatientB_Sacituzumab_Response.pdf"),
         p3, width = 8, height = 5)
}

# ==============================================================================
# PLOT 2: Combined Overview
# ==============================================================================

cat("Generating combined overview...\n")

# Assign treatment category
stats_df <- stats_df %>%
  mutate(
    Treatment = case_when(
      Category %in% olaparib_cats ~ "Olaparib\n(PARP-i)",
      Category %in% durvalumab_cats ~ "Durvalumab\n(anti-PD-L1)",
      Category %in% sacituzumab_cats ~ "Sacituzumab\n(TROP2-ADC)",
      Category %in% c("Stress_Response", "Apoptosis_Evasion") ~ "General\nResistance",
      Category == "Proliferation" ~ "Proliferation\n(Control)",
      TRUE ~ "Other"
    )
  )

# Order by treatment then FC
stats_df <- stats_df %>%
  arrange(Treatment, desc(Log2FC))

p4 <- ggplot(stats_df, aes(x = Log2FC, y = reorder(Gene, Log2FC), fill = Treatment)) +
  geom_bar(stat = "identity") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  facet_grid(Treatment ~ ., scales = "free_y", space = "free_y") +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(
    title = "Patient B: Treatment-Related Gene Changes",
    subtitle = "Bx1 (Lymph Node) → Bx4 (Soft Tissue) | CAUTION: Site difference confounds interpretation",
    x = "Log2 Fold Change",
    y = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "#D73027", size = 9, face = "italic"),
    strip.text.y = element_text(angle = 0, face = "bold", hjust = 0, size = 8),
    axis.text.y = element_text(size = 7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "PatientB_Treatment_Overview.pdf"),
       p4, width = 9, height = 14)

# ==============================================================================
# PLOT 3: Key Findings Summary
# ==============================================================================

cat("Generating key findings summary...\n")

# Select most interesting genes
key_genes <- c(
  # TROP2 target
  "TACSTD2",
  # PD-L1 axis
  "CD274", "PDCD1LG2",
  # Antigen presentation
  "B2M", "HLA-A", "HLA-B",
  # DNA damage
  "ATM", "CHEK1", "RAD51",
  # Stress
  "HSP90AA1", "HSPA1A",
  # Proliferation
  "MKI67", "PCNA"
)

key_genes <- key_genes[key_genes %in% stats_df$Gene]

key_stats <- stats_df %>%
  filter(Gene %in% key_genes) %>%
  mutate(Gene = factor(Gene, levels = key_genes))

p5 <- ggplot(key_stats, aes(x = Gene, y = Log2FC, fill = Log2FC > 0)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_text(aes(label = Sig, y = Log2FC + sign(Log2FC) * 0.1), size = 5) +
  scale_fill_manual(values = c("TRUE" = "#D73027", "FALSE" = "#3575B5"),
                    name = "Direction",
                    labels = c("Decreased", "Increased")) +
  labs(
    title = "Key Treatment-Related Genes: Patient B",
    subtitle = "Bx1 (Lymph Node) → Bx4 (Soft Tissue) | * = adj.p < 0.05",
    x = NULL,
    y = "Log2 Fold Change"
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

ggsave(file.path(output_dir, "PatientB_KeyGenes_Treatment.pdf"),
       p5, width = 10, height = 5)

# ==============================================================================
# Save Data
# ==============================================================================

write.csv(stats_df, file.path(output_dir, "PatientB_Treatment_Gene_Stats.csv"),
          row.names = FALSE)

# ==============================================================================
# Print Summary
# ==============================================================================

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("PATIENT B TREATMENT-SPECIFIC ANALYSIS\n")
cat("=", rep("=", 70), "\n", sep = "")

cat("\n** CRITICAL CAVEAT **\n")
cat("  Bx1 = Lymph Node | Bx4 = Soft Tissue\n")
cat("  Site differences may confound treatment-related changes!\n")

cat("\nTREATMENTS:\n")
cat("  1. Olaparib (PARP-i) - long duration after Bx1\n")
cat("  2. Durvalumab (anti-PD-L1) - long duration after Bx1\n")
cat("  3. Sacituzumab (TROP2-ADC) - short, just before Bx4\n")

cat("\nKEY FINDINGS:\n")

# TROP2
if("TACSTD2" %in% stats_df$Gene) {
  trop2 <- stats_df %>% filter(Gene == "TACSTD2")
  cat(sprintf("\n  TROP2 (Sacituzumab target): Log2FC = %.2f, p = %.3f\n",
              trop2$Log2FC, trop2$P_adj))
  if(trop2$Log2FC < -0.5) {
    cat("    → DECREASED - possible target downregulation/selection\n")
  } else if(trop2$Log2FC > 0.5) {
    cat("    → INCREASED - target still expressed\n")
  } else {
    cat("    → STABLE\n")
  }
}

# PD-L1
if("CD274" %in% stats_df$Gene) {
  pdl1 <- stats_df %>% filter(Gene == "CD274")
  cat(sprintf("\n  PD-L1/CD274 (Durvalumab target): Log2FC = %.2f, p = %.3f\n",
              pdl1$Log2FC, pdl1$P_adj))
}

# Antigen presentation
ap_genes <- c("B2M", "HLA-A", "HLA-B")
ap_genes <- ap_genes[ap_genes %in% stats_df$Gene]
if(length(ap_genes) > 0) {
  ap_data <- stats_df %>% filter(Gene %in% ap_genes)
  cat(sprintf("\n  Antigen Presentation: Mean Log2FC = %.2f\n",
              mean(ap_data$Log2FC)))
  if(mean(ap_data$Log2FC) < -0.3) {
    cat("    → DECREASED - possible immune evasion after Durvalumab\n")
  }
}

# Stress response
stress_genes <- treatment_genes[["Stress_Response"]]
if(length(stress_genes) > 0) {
  stress_data <- stats_df %>% filter(Gene %in% stress_genes)
  cat(sprintf("\n  Stress Response: Mean Log2FC = %.2f\n",
              mean(stress_data$Log2FC, na.rm = TRUE)))
}

cat("\n\nOUTPUT FILES:\n")
cat("  1. PatientB_Olaparib_Response.pdf\n")
cat("  2. PatientB_Durvalumab_Response.pdf\n")
cat("  3. PatientB_Sacituzumab_Response.pdf\n")
cat("  4. PatientB_Treatment_Overview.pdf\n")
cat("  5. PatientB_KeyGenes_Treatment.pdf\n")
cat("  6. PatientB_Treatment_Gene_Stats.csv\n")
cat("\n")
