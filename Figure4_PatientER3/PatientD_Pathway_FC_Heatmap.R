# ==============================================================================
# PATIENT D: PATHWAY FOLD CHANGE HEATMAP
# Cancer Core vs Tumor Nest
# ==============================================================================

library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

set.seed(42)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

OUTPUT_DIR <- "PatientD_Analysis/Figure_Panels/"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

SEURAT_PATH <- "data/CosMx_SMMART_345k_clean.rds"

# Domain definitions
CANCER_CORE_DOMAINS <- c("0", "2", "9")
TUMOR_NEST_DOMAINS <- c("5")

# ==============================================================================
# PATHWAY GENE SETS (Selected pathways only)
# ==============================================================================

PATHWAYS <- list(
  "Estrogen_Response" = c("ESR1", "GATA3", "PGR", "XBP1", "CCND1", "AZGP1", "TFF1", "TFF3", "FOXA1", "AGR2"),
  "Proliferation" = c("MKI67", "PCNA", "CCND1", "CCNE1", "CCNB1", "CDK4", "CDK1", "TOP2A", "TYMS", "MCM2"),
  "Heat_Shock" = c("HSP90AA1", "HSP90AB1", "HSP90B1", "HSPA1A", "HSPA1B", "HSPB1", "HSPA8", "DNAJB1", "PTGES3", "STIP1"),
  "MAPK" = c("FOS", "FOSB", "JUN", "JUNB", "DUSP1", "DUSP4", "DUSP5", "EGR1", "ATF3", "ZFP36"),
  "Interferon" = c("STAT1", "IRF1", "IFI27", "IFIT1", "IFIT3", "ISG15", "MX1", "OAS1", "IFITM1", "IFITM3"),
  "Apoptosis" = c("BCL2", "BCL2L1", "MCL1", "BAX", "BAK1", "BIRC5", "XIAP", "CFLAR", "TNFSF10", "FAS"),
  "Hypoxia" = c("HIF1A", "VEGFA", "SLC2A1", "LDHA", "PGK1", "ENO1", "HILPDA", "BNIP3", "CA9", "ADM"),
  "ECM_Adhesion" = c("FN1", "COL1A1", "COL4A1", "ITGB1", "ITGAV", "ITGB6", "CD63", "TIMP1", "MMP2", "MMP9")
)

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading Seurat object...\n")
obj_all <- readRDS(SEURAT_PATH)
obj_all <- UpdateSeuratObject(obj_all)

all_genes <- rownames(obj_all)

# Filter to Patient D cancer cells
obj_all$Domain <- as.character(obj_all$spatial_domain)

cells_use <- colnames(obj_all)[
  obj_all$Patient == "Patient_D" & 
  obj_all$Lineage == "Cancer" &
  (obj_all$Domain %in% c(CANCER_CORE_DOMAINS, TUMOR_NEST_DOMAINS))
]

obj_subset <- subset(obj_all, cells = cells_use)
obj_subset$Domain_Group <- ifelse(
  obj_subset$Domain %in% CANCER_CORE_DOMAINS,
  "Cancer_Core",
  "Tumor_Nest"
)

cat("Cancer Core cells:", sum(obj_subset$Domain_Group == "Cancer_Core"), "\n")
cat("Tumor Nest cells:", sum(obj_subset$Domain_Group == "Tumor_Nest"), "\n")

# ==============================================================================
# CALCULATE MEAN EXPRESSION PER PATHWAY
# ==============================================================================

cat("\nCalculating pathway expression...\n")

pathway_info <- list()

# Get expression matrix
all_pathway_genes <- unique(unlist(lapply(PATHWAYS, function(x) x[x %in% all_genes])))
expr_matrix <- GetAssayData(obj_subset, layer = "data")[all_pathway_genes, , drop = FALSE]
expr_matrix <- as.matrix(expr_matrix)

# Calculate mean expression for each pathway per cell
pathway_scores <- data.frame(
  cell_id = colnames(obj_subset),
  Domain_Group = obj_subset$Domain_Group
)

for (pathway in names(PATHWAYS)) {
  genes <- PATHWAYS[[pathway]]
  available <- genes[genes %in% all_genes]
  
  pathway_info[[pathway]] <- list(
    available = available,
    n_available = length(available),
    n_total = length(genes)
  )
  
  if (length(available) >= 2) {
    pathway_scores[[pathway]] <- colMeans(expr_matrix[available, , drop = FALSE], na.rm = TRUE)
  }
}

# Get pathway columns
pathway_cols <- names(PATHWAYS)[names(PATHWAYS) %in% colnames(pathway_scores)]

# ==============================================================================
# CALCULATE FOLD CHANGE
# ==============================================================================

# Mean per domain
mean_by_domain <- pathway_scores %>%
  group_by(Domain_Group) %>%
  summarise(across(all_of(pathway_cols), mean, na.rm = TRUE), .groups = "drop")

# Extract values
cancer_core_means <- mean_by_domain %>% filter(Domain_Group == "Cancer_Core") %>% select(-Domain_Group) %>% unlist()
tumor_nest_means <- mean_by_domain %>% filter(Domain_Group == "Tumor_Nest") %>% select(-Domain_Group) %>% unlist()

# Calculate fold change (Tumor Nest / Cancer Core)
fc_df <- data.frame(
  Pathway = pathway_cols,
  Cancer_Core = cancer_core_means,
  Tumor_Nest = tumor_nest_means
) %>%
  mutate(
    FC = Tumor_Nest / Cancer_Core,
    Log2FC = log2(FC),
    n_genes = sapply(Pathway, function(p) pathway_info[[p]]$n_available)
  )

# Statistical test
stat_results <- pathway_scores %>%
  pivot_longer(cols = all_of(pathway_cols), names_to = "Pathway", values_to = "Score") %>%
  group_by(Pathway) %>%
  summarise(
    p_value = wilcox.test(Score ~ Domain_Group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

fc_df <- fc_df %>%
  left_join(stat_results, by = "Pathway") %>%
  arrange(desc(Log2FC))

cat("\nFold change summary:\n")
print(fc_df %>% as.data.frame())

# ==============================================================================
# HEATMAP - SINGLE COLUMN LOG2FC (CLEAN VERSION)
# ==============================================================================

cat("\nCreating fold change heatmap...\n")

# Order by Log2FC
pathway_order <- fc_df %>% arrange(desc(Log2FC)) %>% pull(Pathway)

# Create matrix (single column)
fc_mat <- matrix(fc_df$Log2FC[match(pathway_order, fc_df$Pathway)], ncol = 1)
rownames(fc_mat) <- pathway_order
colnames(fc_mat) <- "Log2FC"

# Color scale - diverging
max_fc <- max(abs(fc_mat), na.rm = TRUE) * 1.1
col_fun <- colorRamp2(c(-max_fc, 0, max_fc), c("#7570B3", "white", "#D51F26"))

pdf(file.path(OUTPUT_DIR, "PatientD_Pathway_FC_Heatmap.pdf"), width = 5, height = 5)

ht <- Heatmap(
  fc_mat,
  name = "Log2FC",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_gp = gpar(fontsize = 12, fontface = "bold"),
  column_title = "Tumor Nest vs Cancer Core",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  width = unit(2, "cm"),
  use_raster = FALSE
)

draw(ht)
dev.off()

cat("✓ Saved PatientD_Pathway_FC_Heatmap.pdf\n")

# ==============================================================================
# HEATMAP - TWO COLUMN (BOTH DOMAINS) - CLEAN
# ==============================================================================

# Create matrix with both domains
expr_mat <- as.matrix(fc_df[match(pathway_order, fc_df$Pathway), c("Cancer_Core", "Tumor_Nest")])
rownames(expr_mat) <- pathway_order

# Row annotation with FC only
row_anno2 <- rowAnnotation(
  Log2FC = fc_df$Log2FC[match(pathway_order, fc_df$Pathway)],
  col = list(Log2FC = colorRamp2(c(-0.3, 0, 0.3), c("#7570B3", "white", "#D51F26"))),
  annotation_name_gp = gpar(fontsize = 9)
)

# Column annotation
col_anno <- HeatmapAnnotation(
  Domain = c("Cancer_Core", "Tumor_Nest"),
  col = list(Domain = c("Cancer_Core" = "#7570B3", "Tumor_Nest" = "#D51F26")),
  show_annotation_name = FALSE
)

pdf(file.path(OUTPUT_DIR, "PatientD_Pathway_Expression_Heatmap.pdf"), width = 5, height = 5)

ht2 <- Heatmap(
  expr_mat,
  name = "Mean Expr",
  col = colorRamp2(c(0, max(expr_mat)/2, max(expr_mat)), c("white", "#FED976", "#D51F26")),
  top_annotation = col_anno,
  right_annotation = row_anno2,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_title = "Pathway Expression by Domain",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  width = unit(3, "cm"),
  use_raster = FALSE
)

draw(ht2)
dev.off()

cat("✓ Saved PatientD_Pathway_Expression_Heatmap.pdf\n")

# ==============================================================================
# BAR PLOT VERSION OF FC
# ==============================================================================

fc_plot <- fc_df %>%
  mutate(
    Pathway = factor(Pathway, levels = pathway_order),
    Direction = ifelse(Log2FC > 0, "Tumor_Nest", "Cancer_Core")
  )

p_fc_bar <- ggplot(fc_plot, aes(x = Pathway, y = Log2FC, fill = Direction)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_fill_manual(
    values = c("Tumor_Nest" = "#D51F26", "Cancer_Core" = "#7570B3"),
    name = "Enriched In"
  ) +
  coord_flip() +
  labs(
    title = "Pathway Fold Change: Tumor Nest vs Cancer Core",
    subtitle = "Log2(Tumor Nest / Cancer Core)",
    x = NULL,
    y = "Log2 Fold Change"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.text.y = element_text(size = 10, face = "bold"),
    legend.position = "top",
    panel.grid.major.y = element_blank()
  )

ggsave(
  file.path(OUTPUT_DIR, "PatientD_Pathway_FC_BarPlot.pdf"),
  p_fc_bar,
  width = 8, height = 6
)
cat("✓ Saved PatientD_Pathway_FC_BarPlot.pdf\n")

# ==============================================================================
# SUPPLEMENTAL: GENE-LEVEL FOLD CHANGE HEATMAP
# ==============================================================================

cat("\nCreating gene-level FC heatmap...\n")

# Calculate FC for each gene
gene_fc_list <- list()

for (pathway in names(PATHWAYS)) {
  genes <- PATHWAYS[[pathway]]
  available <- genes[genes %in% all_genes]
  
  for (gene in available) {
    # Mean expression per domain
    expr_cc <- mean(expr_matrix[gene, obj_subset$Domain_Group == "Cancer_Core"], na.rm = TRUE)
    expr_tn <- mean(expr_matrix[gene, obj_subset$Domain_Group == "Tumor_Nest"], na.rm = TRUE)
    
    gene_fc_list[[length(gene_fc_list) + 1]] <- data.frame(
      Pathway = pathway,
      Gene = gene,
      Cancer_Core = expr_cc,
      Tumor_Nest = expr_tn,
      Log2FC = log2((expr_tn + 0.01) / (expr_cc + 0.01))
    )
  }
}

gene_fc_df <- bind_rows(gene_fc_list)

# Order pathways by their overall FC
gene_fc_df$Pathway <- factor(gene_fc_df$Pathway, levels = pathway_order)

# Order genes within pathway by FC
gene_fc_df <- gene_fc_df %>%
  arrange(Pathway, desc(Log2FC))

# Create gene order
gene_order <- gene_fc_df$Gene

# Create matrix
gene_fc_mat <- matrix(gene_fc_df$Log2FC, ncol = 1)
rownames(gene_fc_mat) <- gene_fc_df$Gene
colnames(gene_fc_mat) <- "Log2FC"

# Pathway annotation
pathway_colors <- setNames(
  scales::hue_pal()(length(unique(gene_fc_df$Pathway))),
  levels(gene_fc_df$Pathway)
)

gene_row_anno <- rowAnnotation(
  Pathway = gene_fc_df$Pathway,
  col = list(Pathway = pathway_colors),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 9)
)

# Color scale
max_gene_fc <- max(abs(gene_fc_mat), na.rm = TRUE) * 1.1
gene_col_fun <- colorRamp2(c(-max_gene_fc, 0, max_gene_fc), c("#7570B3", "white", "#D51F26"))

pdf(file.path(OUTPUT_DIR, "PatientD_Pathway_GeneLevel_FC_Heatmap.pdf"), width = 6, height = 12)

ht_gene <- Heatmap(
  gene_fc_mat,
  name = "Log2FC",
  col = gene_col_fun,
  left_annotation = gene_row_anno,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_title = "Gene-Level Fold Change\n(Tumor Nest vs Cancer Core)",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  row_split = gene_fc_df$Pathway,
  row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  row_title_rot = 0,
  width = unit(2, "cm"),
  use_raster = FALSE
)

draw(ht_gene)
dev.off()

cat("✓ Saved PatientD_Pathway_GeneLevel_FC_Heatmap.pdf\n")

# Save gene-level data
write_csv(gene_fc_df, file.path(OUTPUT_DIR, "PatientD_Pathway_GeneLevel_FC.csv"))
cat("✓ Saved PatientD_Pathway_GeneLevel_FC.csv\n")

# ==============================================================================
# SAVE DATA
# ==============================================================================

# Add genes used
fc_df$Genes_Used <- sapply(fc_df$Pathway, function(p) paste(pathway_info[[p]]$available, collapse = ", "))

write_csv(fc_df, file.path(OUTPUT_DIR, "PatientD_Pathway_FC_Summary.csv"))
cat("✓ Saved PatientD_Pathway_FC_Summary.csv\n")

# ==============================================================================
# PRINT SUMMARY
# ==============================================================================

cat("\n========== FOLD CHANGE SUMMARY ==========\n")
cat("Positive = higher in Tumor Nest | Negative = higher in Cancer Core\n\n")

for (i in 1:nrow(fc_df)) {
  row <- fc_df[i, ]
  cat(sprintf("%-20s Log2FC: %+.3f  (%.2f vs %.2f) [%d genes]\n",
              row$Pathway, row$Log2FC, row$Tumor_Nest, row$Cancer_Core, row$n_genes))
}
