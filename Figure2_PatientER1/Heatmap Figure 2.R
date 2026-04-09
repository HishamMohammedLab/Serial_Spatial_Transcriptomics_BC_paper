# ==============================================================================
# FIGURE 2c,e: LOG2FC HEATMAPS — Patient ER1 (Patient A), Cancer Cells Only
# Bx2 vs Bx4 Comparison
# ==============================================================================

library(Seurat)
library(tidyverse)

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================

OBJ_PATH   <- "data/CosMx_SMMART_345k_clean.rds"
OUTPUT_DIR <- "Figure2_PatientER1/"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

TARGET_PATIENT <- "Patient_A"
EARLY_TP <- "Bx2"
LATE_TP  <- "Bx4"
PSEUDOCOUNT <- 0.001

# --- Panel c genes ---
panel_c_genes <- list(
  "Estrogen Signaling" = c("ESR1", "GATA3", "XBP1", "AR"),
  "TGF\u03b2 / EMT"   = c("TGFBR2", "JAG1", "SNAI1", "SNAI2"),
  "Proliferation"      = c("MKI67")
)

# --- Panel e genes ---
panel_e_genes <- list(
  "IFN Signaling"          = c("STAT1", "MX1", "OAS1", "IFITM1", "BST2", "IFI27"),
  "Antigen\nPresentation"  = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2",
                                "HLA-DRA", "HLA-DPA1"),
  "Immune\nCheckpoints"    = c("CD274", "PDCD1LG2", "FASLG", "LAG3", "LGALS9", "TNFSF10"),
  "Chemokines"             = c("CXCL8", "CXCL2", "CXCL5", "CXCL9", "CXCL10", "CXCL1")
)

# ==============================================================================
# 2. LOAD & SUBSET TO CANCER CELLS
# ==============================================================================

message("Loading Seurat object...")
obj <- readRDS(OBJ_PATH)

message("Subsetting to Patient A cancer cells (Bx2, Bx4)...")
obj_sub <- subset(obj, subset = Patient == TARGET_PATIENT &
                    Timepoint %in% c(EARLY_TP, LATE_TP) &
                    broad_lineage == "Cancer")
message(sprintf("  Cancer cells: %d (Bx2=%d, Bx4=%d)",
                ncol(obj_sub),
                sum(obj_sub$Timepoint == EARLY_TP),
                sum(obj_sub$Timepoint == LATE_TP)))
rm(obj); gc()

# ==============================================================================
# 3. COMPUTE LOG2FC
# ==============================================================================

compute_log2fc <- function(gene_lists, seurat_obj, early, late, pseudo) {
  all_genes <- unlist(gene_lists, use.names = FALSE)
  avail <- all_genes[all_genes %in% rownames(seurat_obj)]

  expr <- GetAssayData(seurat_obj, layer = "data")

  early_cells <- which(seurat_obj$Timepoint == early)
  late_cells  <- which(seurat_obj$Timepoint == late)

  results <- data.frame()
  cat_order <- c()
  for (cat_name in names(gene_lists)) {
    genes <- gene_lists[[cat_name]]
    genes <- genes[genes %in% avail]
    cat_order <- c(cat_order, cat_name)

    for (g in genes) {
      mean_early <- mean(expr[g, early_cells])
      mean_late  <- mean(expr[g, late_cells])
      log2fc <- log2((mean_late + pseudo) / (mean_early + pseudo))

      results <- rbind(results, data.frame(
        Gene = g,
        Category = cat_name,
        Mean_Early = mean_early,
        Mean_Late = mean_late,
        Log2FC = log2fc,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Order genes: within each category, keep original order
  results$Category <- factor(results$Category, levels = cat_order)
  gene_order <- results$Gene
  results$Gene <- factor(results$Gene, levels = rev(gene_order))

  return(results)
}

df_c <- compute_log2fc(panel_c_genes, obj_sub, EARLY_TP, LATE_TP, PSEUDOCOUNT)
df_e <- compute_log2fc(panel_e_genes, obj_sub, EARLY_TP, LATE_TP, PSEUDOCOUNT)

message("\nPanel c Log2FC:")
print(df_c %>% select(Gene, Category, Log2FC) %>% mutate(Log2FC = round(Log2FC, 2)))
message("\nPanel e Log2FC:")
print(df_e %>% select(Gene, Category, Log2FC) %>% mutate(Log2FC = round(Log2FC, 2)))

# ==============================================================================
# 4. PLOTTING FUNCTION
# ==============================================================================

plot_log2fc_heatmap <- function(df, fc_limit = 2.5, panel_label = "") {

  # Clamp Log2FC for color scale
  df$Log2FC_clamped <- pmax(pmin(df$Log2FC, fc_limit), -fc_limit)

  p <- ggplot(df, aes(x = 1, y = Gene, fill = Log2FC_clamped)) +
    geom_tile(color = "white", linewidth = 0.4) +

    # Diverging blue-white-red
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0,
      limits = c(-fc_limit, fc_limit),
      name = paste0(EARLY_TP, "\u2192", LATE_TP, "\nLog2FC"),
      guide = guide_colorbar(
        barwidth = 0.6, barheight = 4,
        frame.colour = "black", ticks.colour = "black",
        title.position = "top"
      )
    ) +

    # Facet by category
    facet_grid(Category ~ ., scales = "free_y", space = "free_y", switch = "y") +

    labs(x = NULL, y = NULL) +

    theme_void(base_size = 11) +
    theme(
      text = element_text(family = "sans"),

      # Gene names on right
      axis.text.y = element_text(face = "italic", color = "black", size = 10,
                                 hjust = 0, margin = margin(l = 4)),

      # No x-axis
      axis.text.x = element_blank(),

      # Category strips on left
      strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 1,
                                        size = 10, color = "black",
                                        margin = margin(r = 4)),
      strip.placement = "outside",
      strip.background = element_blank(),

      # Panel spacing
      panel.spacing = unit(0.3, "lines"),

      # Legend
      legend.position = "right",
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7),
      legend.margin = margin(l = 6),

      plot.margin = margin(10, 10, 10, 10)
    )

  # Add panel label
  if (nchar(panel_label) > 0) {
    p <- p + labs(title = panel_label) +
      theme(plot.title = element_text(face = "bold", size = 14, hjust = 0,
                                       margin = margin(b = 4)))
  }

  return(p)
}

# ==============================================================================
# 5. GENERATE AND SAVE
# ==============================================================================

message("\nGenerating plots...")

p_c <- plot_log2fc_heatmap(df_c, fc_limit = 2.5, panel_label = "c")
p_e <- plot_log2fc_heatmap(df_e, fc_limit = 2.5, panel_label = "e")

# Save Panel c
ggsave(file.path(OUTPUT_DIR, "Fig2c_Log2FC_Heatmap.pdf"),
       p_c, width = 3.2, height = 3.0)
ggsave(file.path(OUTPUT_DIR, "Fig2c_Log2FC_Heatmap.png"),
       p_c, width = 3.2, height = 3.0, dpi = 300)

# Save Panel e
ggsave(file.path(OUTPUT_DIR, "Fig2e_Log2FC_Heatmap.pdf"),
       p_e, width = 3.2, height = 5.5)
ggsave(file.path(OUTPUT_DIR, "Fig2e_Log2FC_Heatmap.png"),
       p_e, width = 3.2, height = 5.5, dpi = 300)

message("\nSaved:")
message("  ", file.path(OUTPUT_DIR, "Fig2c_Log2FC_Heatmap.pdf"))
message("  ", file.path(OUTPUT_DIR, "Fig2e_Log2FC_Heatmap.pdf"))

# Also save underlying data
write.csv(rbind(df_c %>% mutate(Panel = "c"), df_e %>% mutate(Panel = "e")),
          file.path(OUTPUT_DIR, "Fig2ce_Log2FC_Data.csv"), row.names = FALSE)
