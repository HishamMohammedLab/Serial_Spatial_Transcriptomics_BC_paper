# ============================================================================
# Patient C: Comprehensive Investigation
# Check for errors and find novel angles
# ============================================================================

library(Seurat)
library(tidyverse)
library(data.table)
library(patchwork)

set.seed(42)

# ============================================================================
# CLINICAL CONTEXT
# ============================================================================
# - Two biopsies (Bx1, Bx2), both ER+ but different ESR1 mutations
# - Germline BRCA2 mutant and MSH6 mutant
# - Bx2 has loss of: CDKN2B, MTAP, BRCA2, CDKN2A
# - Treatment: Tamoxifen + Everolimus + Paclitaxel (pre-Bx1)
#              then Fulvestrant + Abemaciclib (between Bx1 and Bx2)
# ============================================================================

output_dir <- "PatientC_Analysis/Investigation/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading data...\n")
obj <- readRDS("data/CosMx_SMMART_345k_clean.rds")

# Subset to Patient C
patC <- subset(obj, Patient == "Patient_C")
cat("Patient C cells:", ncol(patC), "\n")

# Check lineages
cat("\nLineage distribution:\n")
print(table(patC$Lineage, patC$Timepoint))

# Get expression matrix
expr_mat <- GetAssayData(patC, layer = "data")

# ============================================================================
# ISSUE 1: ECM GENES IN CANCER CELLS?
# ============================================================================

cat("\n========== ISSUE 1: ECM GENES BY CELL TYPE ==========\n\n")

# ECM genes from the analysis
ecm_genes <- c("COL1A1", "COL1A2", "COL3A1", "FN1", "COL6A1", "COL6A2", "THBS1")
ecm_avail <- ecm_genes[ecm_genes %in% rownames(expr_mat)]

# Calculate mean expression by lineage and timepoint
ecm_by_lineage <- data.frame()

for (gene in ecm_avail) {
  for (lin in unique(patC$Lineage)) {
    for (tp in c("Bx1", "Bx2")) {
      cells <- which(patC$Lineage == lin & patC$Timepoint == tp)
      if (length(cells) > 10) {
        mean_expr <- mean(expr_mat[gene, cells], na.rm = TRUE)
        pct_pos <- 100 * mean(expr_mat[gene, cells] > 0, na.rm = TRUE)
        ecm_by_lineage <- rbind(ecm_by_lineage, data.frame(
          Gene = gene, Lineage = lin, Timepoint = tp,
          Mean_Expr = mean_expr, Pct_Positive = pct_pos, n_cells = length(cells)
        ))
      }
    }
  }
}

cat("ECM gene expression by lineage:\n")
ecm_wide <- ecm_by_lineage %>%
  filter(Lineage %in% c("Cancer", "Fibroblast")) %>%
  select(Gene, Lineage, Timepoint, Mean_Expr) %>%
  pivot_wider(names_from = c(Lineage, Timepoint), values_from = Mean_Expr)
print(ecm_wide)

# Check if cancer cells express COL1A1
cancer_col1a1 <- ecm_by_lineage %>% filter(Gene == "COL1A1", Lineage == "Cancer")
fibro_col1a1 <- ecm_by_lineage %>% filter(Gene == "COL1A1", Lineage == "Fibroblast")

cat("\n--- COL1A1 specifically ---\n")
cat("Cancer cells:\n")
print(cancer_col1a1)
cat("\nFibroblasts:\n")
print(fibro_col1a1)

cat("\nVERDICT: ")
if (mean(cancer_col1a1$Mean_Expr) > 0.5 * mean(fibro_col1a1$Mean_Expr)) {
  cat("WARNING - Cancer cells have substantial COL1A1 expression!\n")
} else {
  cat("OK - COL1A1 is mostly in fibroblasts as expected\n")
}

fwrite(ecm_by_lineage, file.path(output_dir, "ECM_by_Lineage.csv"))

# ============================================================================
# ISSUE 2: LUMINAL GENES - CHANGES IN CANCER VS STROMA
# ============================================================================

cat("\n========== ISSUE 2: LUMINAL GENES BY CELL TYPE ==========\n\n")

luminal_genes <- c("ESR1", "GATA3", "KRT8", "KRT18", "KRT19", "CDH1", "EPCAM", "AGR2")
luminal_avail <- luminal_genes[luminal_genes %in% rownames(expr_mat)]

luminal_by_lineage <- data.frame()

for (gene in luminal_avail) {
  for (lin in c("Cancer", "Fibroblast")) {
    for (tp in c("Bx1", "Bx2")) {
      cells <- which(patC$Lineage == lin & patC$Timepoint == tp)
      if (length(cells) > 10) {
        mean_expr <- mean(expr_mat[gene, cells], na.rm = TRUE)
        luminal_by_lineage <- rbind(luminal_by_lineage, data.frame(
          Gene = gene, Lineage = lin, Timepoint = tp, Mean_Expr = mean_expr
        ))
      }
    }
  }
}

# Calculate fold changes
luminal_fc <- luminal_by_lineage %>%
  pivot_wider(names_from = Timepoint, values_from = Mean_Expr) %>%
  mutate(FC = Bx2 / Bx1, Log2FC = log2(Bx2 / Bx1))

cat("Luminal gene fold changes (Bx2/Bx1) by cell type:\n")
print(luminal_fc %>% arrange(Lineage, desc(abs(Log2FC))))

fwrite(luminal_fc, file.path(output_dir, "Luminal_FC_by_Lineage.csv"))

# ============================================================================
# ISSUE 3: MUTATION-RELATED GENES
# ============================================================================

cat("\n========== ISSUE 3: MUTATION-RELATED GENES ==========\n\n")

# Genes related to mutations: BRCA2, CDKN2A, CDKN2B, MTAP, MSH6
# Also DNA damage response genes
mutation_genes <- c("BRCA2", "BRCA1", "CDKN2A", "CDKN2B", "CDKN1A", "CDKN1B",
                    "MTAP", "MSH6", "MSH2", "MLH1", "RAD51", "ATM", "ATR",
                    "CHEK1", "CHEK2", "TP53", "RB1", "E2F1", "CCND1", "CCNE1",
                    "CDK4", "CDK6", "CDK2")
mutation_avail <- mutation_genes[mutation_genes %in% rownames(expr_mat)]

cat("Available mutation-related genes:", paste(mutation_avail, collapse = ", "), "\n\n")

# Focus on cancer cells only
cancer_cells_bx1 <- which(patC$Lineage == "Cancer" & patC$Timepoint == "Bx1")
cancer_cells_bx2 <- which(patC$Lineage == "Cancer" & patC$Timepoint == "Bx2")

mutation_expr <- data.frame()
for (gene in mutation_avail) {
  bx1_expr <- mean(expr_mat[gene, cancer_cells_bx1], na.rm = TRUE)
  bx2_expr <- mean(expr_mat[gene, cancer_cells_bx2], na.rm = TRUE)
  bx1_pct <- 100 * mean(expr_mat[gene, cancer_cells_bx1] > 0, na.rm = TRUE)
  bx2_pct <- 100 * mean(expr_mat[gene, cancer_cells_bx2] > 0, na.rm = TRUE)

  mutation_expr <- rbind(mutation_expr, data.frame(
    Gene = gene, Bx1_Mean = bx1_expr, Bx2_Mean = bx2_expr,
    Bx1_Pct = bx1_pct, Bx2_Pct = bx2_pct,
    FC = bx2_expr / bx1_expr, Log2FC = log2(bx2_expr / bx1_expr)
  ))
}

mutation_expr <- mutation_expr %>% arrange(desc(abs(Log2FC)))
cat("Mutation-related gene changes in CANCER cells:\n")
print(mutation_expr)

fwrite(mutation_expr, file.path(output_dir, "Mutation_Genes_Cancer.csv"))

# Key finding: CDKN2A/B loss should show in expression
cat("\n--- Key: CDKN2A (p16) and CDKN1A (p21) ---\n")
cat("Clinical: Bx2 has loss of CDKN2A, CDKN2B\n")
cat(sprintf("CDKN2A: Bx1=%.2f, Bx2=%.2f, FC=%.2f\n",
            mutation_expr$Bx1_Mean[mutation_expr$Gene == "CDKN2A"],
            mutation_expr$Bx2_Mean[mutation_expr$Gene == "CDKN2A"],
            mutation_expr$FC[mutation_expr$Gene == "CDKN2A"]))
cat(sprintf("CDKN1A: Bx1=%.2f, Bx2=%.2f, FC=%.2f\n",
            mutation_expr$Bx1_Mean[mutation_expr$Gene == "CDKN1A"],
            mutation_expr$Bx2_Mean[mutation_expr$Gene == "CDKN1A"],
            mutation_expr$FC[mutation_expr$Gene == "CDKN1A"]))

# ============================================================================
# ISSUE 4: TREATMENT-RELEVANT GENES
# ============================================================================

cat("\n========== ISSUE 4: TREATMENT-RELEVANT GENES ==========\n\n")

# Treatment: Fulvestrant (ER degrader) + Abemaciclib (CDK4/6 inhibitor)
# Everolimus (mTOR inhibitor) was pre-Bx1

# ER pathway
er_genes <- c("ESR1", "PGR", "GREB1", "TFF1", "MYC", "CCND1")
# CDK4/6 pathway (Abemaciclib targets)
cdk_genes <- c("CDK4", "CDK6", "CCND1", "CCNE1", "RB1", "E2F1", "E2F2", "CDKN2A", "CDKN1A", "CDKN1B")
# mTOR pathway (Everolimus targets)
mtor_genes <- c("MTOR", "RPS6", "EIF4E", "PTEN", "AKT1", "PIK3CA")
# Proliferation
prolif_genes <- c("MKI67", "PCNA", "TOP2A", "MCM2", "MCM6")

all_treatment <- unique(c(er_genes, cdk_genes, mtor_genes, prolif_genes))
treatment_avail <- all_treatment[all_treatment %in% rownames(expr_mat)]

treatment_expr <- data.frame()
for (gene in treatment_avail) {
  bx1_expr <- mean(expr_mat[gene, cancer_cells_bx1], na.rm = TRUE)
  bx2_expr <- mean(expr_mat[gene, cancer_cells_bx2], na.rm = TRUE)

  treatment_expr <- rbind(treatment_expr, data.frame(
    Gene = gene, Bx1_Mean = bx1_expr, Bx2_Mean = bx2_expr,
    FC = bx2_expr / bx1_expr, Log2FC = log2(bx2_expr / bx1_expr),
    Pathway = case_when(
      gene %in% er_genes ~ "ER_pathway",
      gene %in% cdk_genes ~ "CDK4/6_pathway",
      gene %in% mtor_genes ~ "mTOR_pathway",
      gene %in% prolif_genes ~ "Proliferation",
      TRUE ~ "Other"
    )
  ))
}

treatment_expr <- treatment_expr %>% arrange(Pathway, desc(abs(Log2FC)))
cat("Treatment-relevant gene changes in CANCER cells:\n")
print(treatment_expr)

fwrite(treatment_expr, file.path(output_dir, "Treatment_Genes_Cancer.csv"))

# ============================================================================
# ISSUE 5: INTERFERON SIGNATURE DETAILS
# ============================================================================

cat("\n========== ISSUE 5: INTERFERON SIGNATURE ==========\n\n")

ifn_genes <- c("MX1", "MX2", "OAS1", "OAS2", "OAS3", "ISG15", "ISG20", "IFIT1", "IFIT2", "IFIT3",
               "IFITM1", "IFITM3", "IRF1", "IRF7", "STAT1", "STAT2", "DDX58", "IFIH1",
               "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "PSMB8", "PSMB9")
ifn_avail <- ifn_genes[ifn_genes %in% rownames(expr_mat)]

# Check in cancer cells and all cells
ifn_expr <- data.frame()
for (gene in ifn_avail) {
  # Cancer only
  bx1_cancer <- mean(expr_mat[gene, cancer_cells_bx1], na.rm = TRUE)
  bx2_cancer <- mean(expr_mat[gene, cancer_cells_bx2], na.rm = TRUE)
  # All cells
  all_bx1 <- which(patC$Timepoint == "Bx1")
  all_bx2 <- which(patC$Timepoint == "Bx2")
  bx1_all <- mean(expr_mat[gene, all_bx1], na.rm = TRUE)
  bx2_all <- mean(expr_mat[gene, all_bx2], na.rm = TRUE)

  ifn_expr <- rbind(ifn_expr, data.frame(
    Gene = gene,
    Cancer_Bx1 = bx1_cancer, Cancer_Bx2 = bx2_cancer, Cancer_FC = bx2_cancer/bx1_cancer,
    All_Bx1 = bx1_all, All_Bx2 = bx2_all, All_FC = bx2_all/bx1_all
  ))
}

ifn_expr <- ifn_expr %>%
  mutate(Cancer_Log2FC = log2(Cancer_FC), All_Log2FC = log2(All_FC)) %>%
  arrange(Cancer_Log2FC)

cat("Interferon signature changes:\n")
print(ifn_expr %>% select(Gene, Cancer_Bx1, Cancer_Bx2, Cancer_Log2FC))

cat("\nMean IFN log2FC in cancer cells:", mean(ifn_expr$Cancer_Log2FC, na.rm = TRUE), "\n")

fwrite(ifn_expr, file.path(output_dir, "IFN_Signature.csv"))

# ============================================================================
# NOVEL ANGLE 1: BRCA2 LOSS AND DNA DAMAGE RESPONSE
# ============================================================================

cat("\n========== NOVEL ANGLE 1: DNA DAMAGE RESPONSE ==========\n\n")

# Patient is germline BRCA2 mutant, Bx2 has BRCA2 loss (LOH?)
# Check DNA damage markers
ddr_genes <- c("BRCA1", "BRCA2", "RAD51", "ATM", "ATR", "CHEK1", "CHEK2",
               "H2AFX", "TP53BP1", "PARP1", "XRCC1", "FANCA", "FANCD2")
ddr_avail <- ddr_genes[ddr_genes %in% rownames(expr_mat)]

ddr_expr <- data.frame()
for (gene in ddr_avail) {
  bx1_expr <- mean(expr_mat[gene, cancer_cells_bx1], na.rm = TRUE)
  bx2_expr <- mean(expr_mat[gene, cancer_cells_bx2], na.rm = TRUE)

  ddr_expr <- rbind(ddr_expr, data.frame(
    Gene = gene, Bx1 = bx1_expr, Bx2 = bx2_expr,
    FC = bx2_expr/bx1_expr, Log2FC = log2(bx2_expr/bx1_expr)
  ))
}

cat("DNA damage response genes in cancer cells:\n")
print(ddr_expr %>% arrange(desc(abs(Log2FC))))

fwrite(ddr_expr, file.path(output_dir, "DDR_Genes.csv"))

# ============================================================================
# NOVEL ANGLE 2: MSH6 AND MISMATCH REPAIR
# ============================================================================

cat("\n========== NOVEL ANGLE 2: MISMATCH REPAIR ==========\n\n")

mmr_genes <- c("MSH2", "MSH6", "MLH1", "PMS2", "EXO1")
mmr_avail <- mmr_genes[mmr_genes %in% rownames(expr_mat)]

mmr_expr <- data.frame()
for (gene in mmr_avail) {
  bx1_expr <- mean(expr_mat[gene, cancer_cells_bx1], na.rm = TRUE)
  bx2_expr <- mean(expr_mat[gene, cancer_cells_bx2], na.rm = TRUE)

  mmr_expr <- rbind(mmr_expr, data.frame(
    Gene = gene, Bx1 = bx1_expr, Bx2 = bx2_expr,
    FC = bx2_expr/bx1_expr, Log2FC = log2(bx2_expr/bx1_expr)
  ))
}

cat("Mismatch repair genes (MSH6 germline mutant):\n")
print(mmr_expr)

# ============================================================================
# NOVEL ANGLE 3: CDK4/6 INHIBITOR RESPONSE (Abemaciclib)
# ============================================================================

cat("\n========== NOVEL ANGLE 3: CDK4/6i RESPONSE ==========\n\n")

# Abemaciclib between Bx1 and Bx2
# With CDKN2A loss, CDK4/6 should be more active (less inhibited)
# But Abemaciclib should suppress this

cdk_detailed <- data.frame()
cdk_check <- c("CDK4", "CDK6", "CCND1", "CCNE1", "RB1", "E2F1", "CDKN2A", "CDKN1A", "CDKN1B", "MKI67", "TOP2A")
cdk_check <- cdk_check[cdk_check %in% rownames(expr_mat)]

for (gene in cdk_check) {
  bx1_expr <- mean(expr_mat[gene, cancer_cells_bx1], na.rm = TRUE)
  bx2_expr <- mean(expr_mat[gene, cancer_cells_bx2], na.rm = TRUE)
  bx1_pct <- 100 * mean(expr_mat[gene, cancer_cells_bx1] > 0, na.rm = TRUE)
  bx2_pct <- 100 * mean(expr_mat[gene, cancer_cells_bx2] > 0, na.rm = TRUE)

  cdk_detailed <- rbind(cdk_detailed, data.frame(
    Gene = gene, Bx1_Mean = bx1_expr, Bx2_Mean = bx2_expr,
    Bx1_Pct = bx1_pct, Bx2_Pct = bx2_pct,
    FC = bx2_expr/bx1_expr, Log2FC = log2(bx2_expr/bx1_expr)
  ))
}

cat("CDK4/6 pathway detailed (Abemaciclib response):\n")
print(cdk_detailed %>% arrange(desc(Log2FC)))

fwrite(cdk_detailed, file.path(output_dir, "CDK_Pathway_Detailed.csv"))

# ============================================================================
# NOVEL ANGLE 4: HEAT SHOCK / STRESS RESPONSE
# ============================================================================

cat("\n========== NOVEL ANGLE 4: HEAT SHOCK PROTEINS ==========\n\n")

# Topic_16 shows HSP genes going UP - interesting with treatment
hsp_genes <- c("HSPA1A", "HSPA1B", "HSP90AA1", "HSP90AB1", "HSP90B1",
               "HSPB1", "HSPA8", "DNAJB1", "DNAJC3", "HSPH1")
hsp_avail <- hsp_genes[hsp_genes %in% rownames(expr_mat)]

hsp_expr <- data.frame()
for (gene in hsp_avail) {
  bx1_expr <- mean(expr_mat[gene, cancer_cells_bx1], na.rm = TRUE)
  bx2_expr <- mean(expr_mat[gene, cancer_cells_bx2], na.rm = TRUE)

  hsp_expr <- rbind(hsp_expr, data.frame(
    Gene = gene, Bx1 = bx1_expr, Bx2 = bx2_expr,
    FC = bx2_expr/bx1_expr, Log2FC = log2(bx2_expr/bx1_expr)
  ))
}

cat("Heat shock protein changes:\n")
print(hsp_expr %>% arrange(desc(Log2FC)))

fwrite(hsp_expr, file.path(output_dir, "HSP_Genes.csv"))

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("PATIENT C INVESTIGATION SUMMARY\n")
cat(strrep("=", 70), "\n")

cat("\n1. ECM GENES BY CELL TYPE:\n")
cat("   Check if cancer cells genuinely express ECM or if it's signal bleeding\n")

cat("\n2. MUTATION-RELATED:\n")
cat("   - CDKN2A loss: Check if expression is down\n")
cat("   - BRCA2 germline + LOH: Check DDR pathway\n")
cat("   - MSH6: Check MMR genes\n")

cat("\n3. TREATMENT-RELATED:\n")
cat("   - Fulvestrant: ER pathway changes\n")
cat("   - Abemaciclib: CDK4/6/proliferation\n")
cat("   - Post-Everolimus: mTOR pathway\n")

cat("\n4. KEY FINDINGS FROM TOPICS:\n")
cat("   - Topic_19/16 UP: Luminal + HSP genes\n")
cat("   - Topic_25 DOWN: IFN response\n")
cat("   - Topic_15 DOWN: Stress/VEGFA\n")

cat("\nOutput files saved to:", output_dir, "\n")
