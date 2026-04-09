#!/usr/bin/env Rscript
# ==============================================================================
# RIPLEY'S L FUNCTION: T CELL SPATIAL CLUSTERING
# Patient ER 1 (Patient A) — Bx2 vs Bx4
# RNA (cell type annotation) and Protein (CosMx CD3) modalities
# ==============================================================================

library(Seurat)
library(tidyverse)
library(spatstat)

OUTPUT_DIR <- "Figure2_PatientER1/"
PX_TO_UM <- 0.18
MAX_R <- 1000  # Maximum distance in um
NSIM <- 99     # Monte Carlo simulations for envelope

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading RNA data...\n")
rna_obj <- readRDS("data/CosMx_SMMART_345k_clean.rds")
rna_A <- subset(rna_obj, Patient == "Patient_A" & Timepoint %in% c("Bx2", "Bx4"))
rm(rna_obj); gc()

cat("Loading protein data...\n")
suppressWarnings({ prot_obj <- readRDS("data/CosMx_Protein_210k.rds") })
prot_meta <- prot_obj@meta.data
prot_expr <- GetAssayData(prot_obj, assay = "RNA", layer = "counts")

# Patient A: slide 1 = Bx2, slide 2 = Bx4
keep <- prot_meta$slide_ID_numeric %in% c(1, 2)
prot_meta <- prot_meta[keep, ]
prot_expr <- prot_expr[, keep]
rm(prot_obj); gc()

prot_meta$Timepoint <- ifelse(prot_meta$slide_ID_numeric == 1, "Bx2", "Bx4")
prot_meta$CD3 <- as.numeric(prot_expr["CD3", ])

# ==============================================================================
# 2. IDENTIFY T CELLS
# ==============================================================================

# RNA: cell type annotation
rna_meta <- rna_A@meta.data
cat("RNA T cells: Bx2 =", sum(rna_meta$broad_lineage == "tcell" & rna_meta$Timepoint == "Bx2"),
    ", Bx4 =", sum(rna_meta$broad_lineage == "tcell" & rna_meta$Timepoint == "Bx4"), "\n")

# Protein: per-slide 95th percentile CD3 threshold
thresholds <- prot_meta %>%
  group_by(Timepoint) %>%
  summarise(thresh = quantile(CD3, 0.95), .groups = "drop")
prot_meta <- prot_meta %>%
  left_join(thresholds, by = "Timepoint") %>%
  mutate(is_tcell = CD3 > thresh)
cat("Protein CD3+ T cells: Bx2 =", sum(prot_meta$is_tcell & prot_meta$Timepoint == "Bx2"),
    ", Bx4 =", sum(prot_meta$is_tcell & prot_meta$Timepoint == "Bx4"), "\n")

# ==============================================================================
# 3. RIPLEY'S L FUNCTION
# ==============================================================================

calculate_ripleys_L <- function(df, label, max_r, nsim) {
  cat("  ", label, "(n =", nrow(df), ")...\n")
  xr <- range(df$x); yr <- range(df$y)
  buf <- 0.01 * max(diff(xr), diff(yr))
  win <- owin(xrange = c(xr[1] - buf, xr[2] + buf),
              yrange = c(yr[1] - buf, yr[2] + buf))
  pp <- ppp(x = df$x, y = df$y, window = win)
  r_seq <- seq(0, max_r, length.out = 100)
  L_env <- envelope(pp, Lest, r = r_seq, correction = "Ripley",
                    nsim = nsim, verbose = FALSE)
  data.frame(r = L_env$r, L_centered = L_env$obs - L_env$r,
             L_lo = L_env$lo - L_env$r, L_hi = L_env$hi - L_env$r,
             Label = label)
}

# RNA
cat("\nRNA Ripley's L:\n")
L_rna_bx2 <- calculate_ripleys_L(
  rna_meta %>% filter(broad_lineage == "tcell", Timepoint == "Bx2") %>%
    mutate(x = x_global_px * PX_TO_UM, y = y_global_px * PX_TO_UM), "Bx2", MAX_R, NSIM)
L_rna_bx4 <- calculate_ripleys_L(
  rna_meta %>% filter(broad_lineage == "tcell", Timepoint == "Bx4") %>%
    mutate(x = x_global_px * PX_TO_UM, y = y_global_px * PX_TO_UM), "Bx4", MAX_R, NSIM)
L_rna <- bind_rows(L_rna_bx2, L_rna_bx4) %>% mutate(Modality = "RNA")

# Protein
cat("\nProtein Ripley's L:\n")
L_prot_bx2 <- calculate_ripleys_L(
  prot_meta %>% filter(is_tcell, Timepoint == "Bx2") %>%
    mutate(x = x_slide_mm * 1000, y = y_slide_mm * 1000), "Bx2", MAX_R, NSIM)
L_prot_bx4 <- calculate_ripleys_L(
  prot_meta %>% filter(is_tcell, Timepoint == "Bx4") %>%
    mutate(x = x_slide_mm * 1000, y = y_slide_mm * 1000), "Bx4", MAX_R, NSIM)
L_prot <- bind_rows(L_prot_bx2, L_prot_bx4) %>% mutate(Modality = "Protein (CD3+)")

# ==============================================================================
# 4. RESULTS
# ==============================================================================

cat("\n=== RESULTS ===\n")
cat("RNA   Bx2: max L(r)-r =", round(max(L_rna_bx2$L_centered), 1), "um\n")
cat("RNA   Bx4: max L(r)-r =", round(max(L_rna_bx4$L_centered), 1), "um\n")
cat("Prot  Bx2: max L(r)-r =", round(max(L_prot_bx2$L_centered), 1), "um\n")
cat("Prot  Bx4: max L(r)-r =", round(max(L_prot_bx4$L_centered), 1), "um\n")

# ==============================================================================
# 5. PLOT
# ==============================================================================

L_all <- bind_rows(L_rna, L_prot)
L_all$Modality <- factor(L_all$Modality, levels = c("RNA", "Protein (CD3+)"))
tp_colors <- c("Bx2" = "#2166AC", "Bx4" = "#B2182B")

p <- ggplot(L_all, aes(x = r, y = L_centered, color = Label)) +
  geom_ribbon(aes(ymin = L_lo, ymax = L_hi, fill = Label), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  facet_wrap(~ Modality) +
  scale_color_manual(values = tp_colors, name = NULL) +
  scale_fill_manual(values = tp_colors, guide = "none") +
  scale_x_continuous(limits = c(0, MAX_R), breaks = seq(0, MAX_R, 250)) +
  labs(x = expression("Distance r (" * mu * "m)"), y = "L(r) - r") +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(family = "sans"),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_blank(),
    legend.position = c(0.92, 0.85),
    legend.background = element_blank(),
    legend.text = element_text(size = 10),
    axis.line = element_line(linewidth = 0.4),
    plot.margin = margin(10, 15, 10, 10)
  )

ggsave(file.path(OUTPUT_DIR, "RipleysL_Tcells_RNA_Protein.pdf"), p, width = 9, height = 4)
ggsave(file.path(OUTPUT_DIR, "RipleysL_Tcells_RNA_Protein.png"), p, width = 9, height = 4, dpi = 300)

# Save data
write.csv(L_all, file.path(OUTPUT_DIR, "RipleysL_Tcells_Data.csv"), row.names = FALSE)
cat("\nDone.\n")
