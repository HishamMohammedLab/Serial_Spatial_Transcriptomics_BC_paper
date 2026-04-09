# ==============================================================================
# Patient D - Pathway Analysis: Cancer Core vs Tumor Nest
# Purpose: Compare pathway activity between domain groups in cancer cells
# ==============================================================================

library(Seurat)
library(tidyverse)
library(patchwork)

# ------------------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------------------

SEURAT_PATH <- "data/CosMx_SMMART_345k_clean.rds"
OUTPUT_DIR <- "PatientD_Analysis/Pathway_Analysis/"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "Supplemental_GeneBarplots"), showWarnings = FALSE)

# Domain mapping
DOMAIN_TO_GROUP <- c(
  "0" = "Cancer_Core",
  "1" = "Diploid_Epithelial",
  "2" = "Cancer_Core",
  "3" = "Stromal_Interface",
  "4" = "Diploid_Epithelial",
  "5" = "Tumor_Nest",
  "6" = "Diploid_Epithelial",
  "7" = "Diploid_Epithelial",
  "8" = "Diploid_Epithelial",
  "9" = "Cancer_Core"
)

# Colors for comparison
DOMAIN_COLORS <- c(
  "Cancer_Core" = "#D51F26",
  "Tumor_Nest"  = "#7570B3"
)

# ------------------------------------------------------------------------------
# DEFINE PATHWAY GENE SETS
# ------------------------------------------------------------------------------

# Comprehensive pathway definitions based on literature + DE genes provided
PATHWAY_GENES <- list(
  
  # Estrogen Signaling
  Estrogen_Response = c("ESR1", "GATA3", "PGR", "TFF1", "TFF3", "XBP1", 
                        "FOXA1", "AGR2", "CA12", "PDZK1", "CCND1", "MYC",
                        "IGFBP4", "IGFBP5", "STC2", "GREB1"),
  
  # Proliferation/Cell Cycle
  Proliferation = c("MKI67", "PCNA", "CCND1", "CCNE1", "CDK4", "CDK2",
                    "STMN1", "TYMS", "TOP2A", "MCM2", "MCM7", "CENPF",
                    "BIRC5", "PLK1", "AURKA", "BUB1"),
  
  # MAPK/ERK Signaling (including negative regulators)
  MAPK_Signaling = c("DUSP1", "DUSP4", "DUSP5", "DUSP6", "SPRY1", "SPRY2", 
                     "SPRY4", "EGR1", "FOS", "FOSB", "JUN", "JUNB",
                     "IER3", "ZFP36", "NR4A1", "ATF3"),
  
  # Hypoxia Response
  Hypoxia = c("VEGFA", "SLC2A1", "LDHA", "PGK1", "ENO1", "ALDOA",
              "BNIP3", "CA9", "HILPDA", "NDRG1", "PDK1", "HK2",
              "GLUT1", "DDIT4", "EGLN3", "ADM"),
  
  # Heat Shock/Stress Response
  Heat_Shock = c("HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B", "HSPB1", 
                 "HSPA8", "HSPH1", "DNAJB1", "BAG3", "HSPA5",
                 "HSP90B1", "TRAP1", "HSPD1", "HSPE1"),
  
  # Interferon Response
  Interferon_Response = c("ISG15", "ISG20", "IFI27", "IFIT1", "IFIT2", "IFIT3",
                          "MX1", "MX2", "OAS1", "OAS2", "STAT1", "IRF1",
                          "IRF7", "IFITM1", "IFITM3", "BST2"),
  
  # Antigen Presentation/MHC
  Antigen_Presentation = c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F",
                           "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8",
                           "PSMB9", "NLRC5", "CD74", "CIITA"),
  
  # Apoptosis
  Apoptosis = c("BCL2", "BCL2L1", "MCL1", "BAX", "BAK1", "BID",
                "BIRC5", "XIAP", "CASP3", "CASP8", "CASP9",
                "PARP1", "TP53", "BBC3", "PMAIP1"),
  
  # PI3K/AKT Signaling
  PI3K_AKT = c("AKT1", "PIK3CA", "PTEN", "MTOR", "RPS6KB1", "EIF4EBP1",
               "FOXO1", "FOXO3", "GSK3B", "TSC1", "TSC2", "RHEB"),
  
  # TGF-beta Signaling
  TGFB_Signaling = c("TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2",
                     "SMAD2", "SMAD3", "SMAD4", "SMAD7", "SERPINE1",
                     "CTGF", "COL1A1", "FN1", "ACTA2"),
  
  # ECM/Adhesion
  ECM_Adhesion = c("FN1", "COL1A1", "COL3A1", "COL4A1", "COL18A1",
                   "ITGB1", "ITGB4", "ITGB6", "ITGA6", "ITGA2",
                   "DDR1", "DDR2", "CD44", "SDC1", "SDC4"),
  
  # Lipid Metabolism
  Lipid_Metabolism = c("FASN", "SREBF1", "SREBF2", "SCD", "ACACA",
                       "HMGCR", "HMGCS1", "LDLR", "FABP5", "FABP4",
                       "ACLY", "DGAT1", "PLIN2"),
  
  # Oxidative Stress
  Oxidative_Stress = c("SOD1", "SOD2", "CAT", "GPX1", "GPX4", "PRDX1",
                       "PRDX2", "TXN", "TXNRD1", "NQO1", "HMOX1",
                       "GCLC", "GCLM", "GSR", "GSTP1"),
  
  # NF-kB Signaling
  NFkB_Signaling = c("NFKB1", "NFKB2", "RELA", "RELB", "NFKBIA",
                     "NFKBIB", "IKBKB", "IKBKG", "TNF", "IL1B",
                     "IL6", "CCL2", "CXCL8", "BIRC3"),
  
  # Wnt Signaling
  Wnt_Signaling = c("WNT2B", "WNT3A", "WNT5A", "CTNNB1", "APC",
                    "AXIN1", "AXIN2", "LRP5", "LRP6", "FZD1",
                    "TCF7", "LEF1", "MYC", "CCND1", "RNF43"),
  
  # Epithelial Identity
  Epithelial_Identity = c("KRT7", "KRT8", "KRT18", "KRT19", "KRT80",
                          "EPCAM", "CDH1", "CLDN3", "CLDN4", "CLDN7",
                          "TJP1", "OCLN", "MUC1", "ELF3", "GRHL2"),
  
  # Luminal Differentiation
  Luminal_Markers = c("ESR1", "GATA3", "FOXA1", "KRT8", "KRT18",
                      "KRT19", "XBP1", "AR", "SPDEF", "ELF5",
                      "CITED1", "PIP", "AGR2", "TFF1")
)

# ------------------------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------------------------

cat("Loading Seurat object...\n")
seu <- readRDS(SEURAT_PATH)
cat("Total cells:", ncol(seu), "\n")

# Extract metadata
meta <- seu@meta.data

# Filter to Patient D cancer cells only
patD_cancer_idx <- which(meta$Patient == "Patient_D" & meta$Lineage == "Cancer")
cat("Patient D cancer cells:", length(patD_cancer_idx), "\n")

# Add domain groups
meta$spatial_domain_chr <- as.character(meta$spatial_domain)
meta$Domain_Group <- DOMAIN_TO_GROUP[meta$spatial_domain_chr]

# Filter to Cancer Core and Tumor Nest only
meta_patD_cancer <- meta[patD_cancer_idx, ]
compare_idx <- which(meta_patD_cancer$Domain_Group %in% c("Cancer_Core", "Tumor_Nest"))
meta_compare <- meta_patD_cancer[compare_idx, ]

cat("Cancer Core cells:", sum(meta_compare$Domain_Group == "Cancer_Core"), "\n")
cat("Tumor Nest cells:", sum(meta_compare$Domain_Group == "Tumor_Nest"), "\n")

# Get cell names for expression extraction
cells_compare <- rownames(meta_compare)

# ------------------------------------------------------------------------------
# CHECK GENE AVAILABILITY
# ------------------------------------------------------------------------------

cat("\nChecking gene availability in dataset...\n")

all_genes <- rownames(seu)

pathway_genes_available <- lapply(PATHWAY_GENES, function(genes) {
  available <- genes[genes %in% all_genes]
  return(available)
})

# Report availability
cat("\nGenes available per pathway:\n")
for (pw in names(pathway_genes_available)) {
  n_total <- length(PATHWAY_GENES[[pw]])
  n_avail <- length(pathway_genes_available[[pw]])
  cat(sprintf("  %s: %d/%d (%.0f%%)\n", pw, n_avail, n_total, 100*n_avail/n_total))
}

# Filter to pathways with at least 3 genes
pathway_genes_available <- pathway_genes_available[sapply(pathway_genes_available, length) >= 3]
cat("\nPathways with >= 3 genes:", length(pathway_genes_available), "\n")

# ------------------------------------------------------------------------------
# EXTRACT EXPRESSION DATA
# ------------------------------------------------------------------------------

cat("\nExtracting expression data...\n")

# Get normalized expression for comparison cells
expr_data <- GetAssayData(seu, layer = "data")[, cells_compare]
cat("Expression matrix:", nrow(expr_data), "genes x", ncol(expr_data), "cells\n")

# ------------------------------------------------------------------------------
# CALCULATE PATHWAY SCORES
# ------------------------------------------------------------------------------

cat("\nCalculating pathway scores...\n")

# Function to calculate pathway score (mean z-scored expression)
calc_pathway_score <- function(expr_mat, genes) {
  genes_present <- genes[genes %in% rownames(expr_mat)]
  if (length(genes_present) < 2) return(rep(NA, ncol(expr_mat)))
  
  # Extract expression for pathway genes
  pathway_expr <- as.matrix(expr_mat[genes_present, , drop = FALSE])
  
  # Z-score each gene across cells
  pathway_z <- t(apply(pathway_expr, 1, function(x) {
    if (sd(x) == 0) return(rep(0, length(x)))
    (x - mean(x)) / sd(x)
  }))
  
  # Mean across genes for each cell
  colMeans(pathway_z, na.rm = TRUE)
}

# Calculate scores for all pathways
pathway_scores <- data.frame(
  cell_id = cells_compare,
  Domain_Group = meta_compare$Domain_Group,
  Timepoint = meta_compare$Timepoint
)

for (pw in names(pathway_genes_available)) {
  pathway_scores[[pw]] <- calc_pathway_score(expr_data, pathway_genes_available[[pw]])
}

# ------------------------------------------------------------------------------
# STATISTICAL COMPARISON
# ------------------------------------------------------------------------------

cat("\nPerforming statistical comparisons...\n")

# Compare Cancer Core vs Tumor Nest for each pathway
comparison_results <- data.frame(
  Pathway = character(),
  Cancer_Core_Mean = numeric(),
  Tumor_Nest_Mean = numeric(),
  Difference = numeric(),
  pvalue = numeric(),
  Higher_In = character(),
  stringsAsFactors = FALSE
)

for (pw in names(pathway_genes_available)) {
  cc_scores <- pathway_scores[[pw]][pathway_scores$Domain_Group == "Cancer_Core"]
  tn_scores <- pathway_scores[[pw]][pathway_scores$Domain_Group == "Tumor_Nest"]
  
  cc_mean <- mean(cc_scores, na.rm = TRUE)
  tn_mean <- mean(tn_scores, na.rm = TRUE)
  
  # Wilcoxon test
  wtest <- wilcox.test(cc_scores, tn_scores)
  
  comparison_results <- rbind(comparison_results, data.frame(
    Pathway = pw,
    Cancer_Core_Mean = cc_mean,
    Tumor_Nest_Mean = tn_mean,
    Difference = tn_mean - cc_mean,
    pvalue = wtest$p.value,
    Higher_In = ifelse(tn_mean > cc_mean, "Tumor_Nest", "Cancer_Core"),
    stringsAsFactors = FALSE
  ))
}

# Adjust p-values
comparison_results$padj <- p.adjust(comparison_results$pvalue, method = "BH")
comparison_results$Significant <- comparison_results$padj < 0.05

# Sort by absolute difference
comparison_results <- comparison_results %>% arrange(desc(abs(Difference)))

cat("\nPathway Comparison Results:\n")
print(comparison_results)

# Save results table
write.csv(comparison_results, 
          file.path(OUTPUT_DIR, "Pathway_Comparison_Results.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# MAIN FIGURE: PATHWAY BAR PLOT
# ------------------------------------------------------------------------------

cat("\nGenerating pathway bar plot...\n")

# Prepare data for plotting
plot_data <- comparison_results %>%
  mutate(
    Pathway_Label = gsub("_", " ", Pathway),
    # Color by direction if significant, grey if not
    Fill_Category = case_when(
      padj >= 0.05 ~ "Not Significant",
      Difference > 0 ~ "Higher in Tumor Nest",
      TRUE ~ "Higher in Cancer Core"
    )
  ) %>%
  arrange(Difference) %>%
  mutate(Pathway_Label = factor(Pathway_Label, levels = Pathway_Label))

# Bar plot showing difference (Tumor Nest - Cancer Core)
p_pathway <- ggplot(plot_data, aes(x = Difference, y = Pathway_Label, fill = Fill_Category)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("Higher in Cancer Core" = "#D51F26", 
                               "Higher in Tumor Nest" = "#7570B3",
                               "Not Significant" = "grey70")) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10)
  ) +
  labs(
    title = "Pathway Activity: Tumor Nest vs Cancer Core",
    subtitle = "Patient D Cancer Cells",
    x = "Mean Score Difference (Tumor Nest - Cancer Core)",
    y = NULL
  )

ggsave(
  filename = file.path(OUTPUT_DIR, "Pathway_Comparison_Barplot.pdf"),
  plot = p_pathway,
  width = 10, height = 8
)
cat("Saved: Pathway_Comparison_Barplot.pdf\n")

# ------------------------------------------------------------------------------
# SUPPLEMENTAL: INDIVIDUAL GENE BAR PLOTS PER PATHWAY
# ------------------------------------------------------------------------------

cat("\nGenerating individual gene bar plots for supplementary...\n")

# Function to create gene bar plot for a pathway
create_gene_barplot <- function(pathway_name, genes, expr_mat, meta_df) {
  
  genes_present <- genes[genes %in% rownames(expr_mat)]
  if (length(genes_present) == 0) return(NULL)
  
  # Calculate mean expression per gene per domain
  gene_means <- data.frame(
    Gene = character(),
    Domain = character(),
    Mean_Expr = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (gene in genes_present) {
    cc_expr <- expr_mat[gene, meta_df$Domain_Group == "Cancer_Core"]
    tn_expr <- expr_mat[gene, meta_df$Domain_Group == "Tumor_Nest"]
    
    gene_means <- rbind(gene_means, data.frame(
      Gene = gene,
      Domain = c("Cancer_Core", "Tumor_Nest"),
      Mean_Expr = c(mean(cc_expr), mean(tn_expr)),
      stringsAsFactors = FALSE
    ))
  }
  
  # Calculate fold change for ordering
  gene_fc <- gene_means %>%
    pivot_wider(names_from = Domain, values_from = Mean_Expr) %>%
    mutate(FC = Tumor_Nest - Cancer_Core) %>%
    arrange(FC)
  
  gene_means$Gene <- factor(gene_means$Gene, levels = gene_fc$Gene)
  gene_means$Domain <- gsub("_", " ", gene_means$Domain)
  
  # Create plot
  p <- ggplot(gene_means, aes(x = Gene, y = Mean_Expr, fill = Domain)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = c("Cancer Core" = "#D51F26", "Tumor Nest" = "#7570B3")) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
    ) +
    labs(
      title = paste0(gsub("_", " ", pathway_name), " - Gene Expression"),
      subtitle = paste0("n=", length(genes_present), " genes"),
      x = NULL,
      y = "Mean Expression (normalized)"
    )
  
  return(p)
}

# Generate and save gene barplots for each pathway
for (pw in names(pathway_genes_available)) {
  cat("  Processing:", pw, "\n")
  
  p_genes <- create_gene_barplot(
    pathway_name = pw,
    genes = pathway_genes_available[[pw]],
    expr_mat = expr_data,
    meta_df = meta_compare
  )
  
  if (!is.null(p_genes)) {
    # Adjust width based on number of genes
    n_genes <- length(pathway_genes_available[[pw]])
    plot_width <- max(6, min(16, 4 + n_genes * 0.5))
    
    ggsave(
      filename = file.path(OUTPUT_DIR, "Supplemental_GeneBarplots", 
                           paste0("GeneBarplot_", pw, ".pdf")),
      plot = p_genes,
      width = plot_width, height = 5
    )
  }
}

cat("Saved: Supplemental gene barplots\n")

# ------------------------------------------------------------------------------
# COMBINED SUPPLEMENTAL FIGURE
# ------------------------------------------------------------------------------

cat("\nCreating combined supplemental figure...\n")

# Select key pathways for combined figure
key_pathways <- c("Estrogen_Response", "Proliferation", "MAPK_Signaling", 
                  "Hypoxia", "Heat_Shock", "Interferon_Response")
key_pathways <- key_pathways[key_pathways %in% names(pathway_genes_available)]

plots_list <- list()
for (pw in key_pathways) {
  plots_list[[pw]] <- create_gene_barplot(
    pathway_name = pw,
    genes = pathway_genes_available[[pw]],
    expr_mat = expr_data,
    meta_df = meta_compare
  )
}

# Combine with patchwork
if (length(plots_list) > 0) {
  p_combined <- wrap_plots(plots_list, ncol = 2) +
    plot_annotation(
      title = "Supplementary: Individual Gene Expression by Pathway",
      subtitle = "Patient D Cancer Cells: Cancer Core vs Tumor Nest",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5)
      )
    )
  
  ggsave(
    filename = file.path(OUTPUT_DIR, "Supplemental_KeyPathways_Combined.pdf"),
    plot = p_combined,
    width = 16, height = 12
  )
  cat("Saved: Supplemental_KeyPathways_Combined.pdf\n")
}

# ------------------------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------------------------

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Cells analyzed:\n")
cat("  Cancer Core:", sum(meta_compare$Domain_Group == "Cancer_Core"), "\n")
cat("  Tumor Nest:", sum(meta_compare$Domain_Group == "Tumor_Nest"), "\n")

cat("\nSignificant pathways (padj < 0.05):\n")
sig_pathways <- comparison_results %>% filter(Significant)
if (nrow(sig_pathways) > 0) {
  for (i in 1:nrow(sig_pathways)) {
    cat(sprintf("  %s: %s (diff=%.3f, padj=%.2e)\n",
                sig_pathways$Pathway[i],
                sig_pathways$Higher_In[i],
                sig_pathways$Difference[i],
                sig_pathways$padj[i]))
  }
} else {
  cat("  None\n")
}

cat("\nOutput files:\n")
cat("  - Pathway_Comparison_Results.csv\n")
cat("  - Pathway_Comparison_Barplot.pdf\n")
cat("  - Supplemental_GeneBarplots/ (individual pathway gene plots)\n")
cat("  - Supplemental_KeyPathways_Combined.pdf\n")

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("DONE! Output in:", OUTPUT_DIR, "\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
