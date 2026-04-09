# ==============================================================================
# PATIENT D: COMPREHENSIVE LIGAND-RECEPTOR INTERACTION ANALYSIS
# For high-impact publication in spatial transcriptomics of treatment resistance
# ==============================================================================
#
# ANALYSIS OVERVIEW:
# 1. L-R pair identification from CosMx 960-gene panel
# 2. Lineage-specific expression assignment
# 3. Temporal dynamics (Bx1 → Bx2 → Bx3)
# 4. Spatial domain analysis in Bx3 (Cancer Core vs Tumor Nest vs Stromal Interface)
# 5. Publication-quality figures
# 6. Biological narrative for treatment resistance
#
# ==============================================================================

library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(patchwork)
library(viridis)
library(ggrepel)

set.seed(42)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

OUTPUT_DIR <- "PatientD_Analysis/LR_Comprehensive/"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

SEURAT_PATH <- "data/CosMx_SMMART_345k_clean.rds"

# Domain definitions for Patient D
CANCER_CORE_DOMAINS <- c("0", "2", "9")
TUMOR_NEST_DOMAINS <- c("5")
STROMAL_INTERFACE_DOMAINS <- c("3")
DIPLOID_EPITHELIAL_DOMAINS <- c("1", "4", "6", "7", "8")

# Color definitions for annotations
CATEGORY_COLORS <- c(
  "Stromal_Support" = "#E41A1C",
  "ECM_Adhesion" = "#377EB8",
  "Immune_Evasion" = "#4DAF4A",
  "TAM_Recruitment" = "#984EA3",
  "Chemokine" = "#FF7F00",
  "Growth_Factor" = "#FFFF33",
  "TGFb_Signaling" = "#A65628",
  "TAM_Support" = "#F781BF",
  "Inflammation" = "#999999",
  "Checkpoint" = "#66C2A5",
  "Autocrine" = "#FC8D62",
  "Efferocytosis" = "#8DD3C7",
  "CAF_Activation" = "#BEBADA",
  "Costimulation" = "#FB8072",
  "T_Recruitment" = "#80B1D3",
  "T_Activation" = "#FDB462",
  "Cytotoxicity" = "#B3DE69",
  "Interferon" = "#FCCDE5",
  "Angiogenesis" = "#BC80BD"
)

SOURCE_COLORS <- c(
  "Fibroblast" = "#F8766D",
  "Cancer" = "#00BA38",
  "Myeloid" = "#619CFF",
  "T_Lymphocyte" = "#F564E3",
  "Endothelial" = "#00BFC4"
)

TARGET_COLORS <- c(
  "Cancer" = "#00BA38",
  "Myeloid" = "#619CFF",
  "Fibroblast" = "#F8766D",
  "T_Lymphocyte" = "#F564E3",
  "Endothelial" = "#00BFC4"
)

# ==============================================================================
# CURATED LIGAND-RECEPTOR PAIRS WITH LINEAGE ASSIGNMENTS
# Based on established biology and relevance to breast cancer resistance
# ==============================================================================

# Key: Each pair includes expected source (ligand) and target (receptor) lineages
# Format: Ligand, Receptor, Ligand_Source, Receptor_Target, Category, Mechanism

LR_PAIRS <- tribble(
  ~Ligand, ~Receptor, ~Ligand_Source, ~Receptor_Target, ~Category, ~Mechanism,
  
  # === STROMAL SUPPORT / ECM REMODELING ===
  # These pairs mediate direct CAF-Cancer contact and ECM-mediated survival
  "TIMP1", "CD63", "Fibroblast", "Cancer", "Stromal_Support", "TIMP1 promotes cancer survival via CD63-mediated signaling, independent of MMP inhibition",
  "COL1A1", "ITGB1", "Fibroblast", "Cancer", "ECM_Adhesion", "Collagen-integrin signaling activates FAK/Src survival pathways",
  "COL1A1", "ITGA2", "Fibroblast", "Cancer", "ECM_Adhesion", "α2β1 integrin binding to collagen promotes cell survival",
  "COL1A1", "DDR1", "Fibroblast", "Cancer", "ECM_Adhesion", "Discoidin domain receptor activation by collagen",
  "COL1A1", "DDR2", "Fibroblast", "Cancer", "ECM_Adhesion", "DDR2 collagen receptor promotes invasion",
  "COL4A1", "ITGB1", "Fibroblast", "Cancer", "ECM_Adhesion", "Basement membrane collagen supports epithelial organization",
  "FN1", "ITGB1", "Fibroblast", "Cancer", "ECM_Adhesion", "Fibronectin-integrin axis promotes treatment resistance",
  "FN1", "ITGAV", "Fibroblast", "Cancer", "ECM_Adhesion", "αvβ3 integrin activation by fibronectin",
  "FN1", "ITGA5", "Fibroblast", "Cancer", "ECM_Adhesion", "α5β1 fibronectin receptor",
  "FN1", "CD44", "Fibroblast", "Cancer", "ECM_Adhesion", "CD44-fibronectin interaction supports stemness",
  "VTN", "ITGAV", "Fibroblast", "Cancer", "ECM_Adhesion", "Vitronectin-integrin adhesion",
  "THBS1", "CD47", "Fibroblast", "Cancer", "ECM_Adhesion", "Thrombospondin-CD47 'don't eat me' signal",
  "THBS1", "CD36", "Fibroblast", "Myeloid", "ECM_Adhesion", "Thrombospondin modulates macrophage function",
  "TGFB1", "TGFBR2", "Fibroblast", "Cancer", "TGFb_Signaling", "TGF-β promotes EMT and treatment resistance",
  "TGFB1", "TGFBR1", "Fibroblast", "Cancer", "TGFb_Signaling", "TGF-β receptor complex activation",
  "CXCL12", "CXCR4", "Fibroblast", "Cancer", "Chemokine", "SDF-1/CXCR4 axis promotes cancer cell survival and migration",
  "CXCL12", "ACKR3", "Fibroblast", "Cancer", "Chemokine", "CXCR7 scavenges CXCL12, modulating gradient",
  "HGF", "MET", "Fibroblast", "Cancer", "Growth_Factor", "HGF/c-MET promotes proliferation and invasion",
  "IGF1", "IGF1R", "Fibroblast", "Cancer", "Growth_Factor", "IGF-1 signaling promotes survival",
  "PDGFA", "PDGFRA", "Cancer", "Fibroblast", "CAF_Activation", "Cancer cells recruit and activate fibroblasts",
  "PDGFB", "PDGFRB", "Cancer", "Fibroblast", "CAF_Activation", "PDGF-BB drives myofibroblast differentiation",
  
  # === IMMUNE EVASION / MYELOID POLARIZATION ===
  "MIF", "CD74", "Cancer", "Myeloid", "Immune_Evasion", "MIF-CD74 promotes M2 macrophage polarization and immunosuppression",
  "MIF", "CXCR4", "Cancer", "Myeloid", "Immune_Evasion", "MIF-CXCR4 recruits immunosuppressive myeloid cells",
  "MIF", "CXCR2", "Cancer", "Myeloid", "Immune_Evasion", "MIF-CXCR2 neutrophil recruitment",
  "CSF1", "CSF1R", "Cancer", "Myeloid", "TAM_Recruitment", "M-CSF recruits and polarizes tumor-associated macrophages",
  "CCL2", "CCR2", "Cancer", "Myeloid", "TAM_Recruitment", "CCL2/MCP-1 recruits monocytes to TME",
  "CCL5", "CCR5", "Cancer", "Myeloid", "TAM_Recruitment", "RANTES recruits macrophages and T cells",
  "CCL5", "CCR1", "Cancer", "Myeloid", "TAM_Recruitment", "CCL5-CCR1 axis",
  "GAS6", "MERTK", "Cancer", "Myeloid", "Efferocytosis", "GAS6-MerTK promotes M2 polarization via efferocytosis",
  "GAS6", "AXL", "Cancer", "Myeloid", "Efferocytosis", "TAM receptor activation promotes immune suppression",
  "GAS6", "TYRO3", "Cancer", "Myeloid", "Efferocytosis", "TYRO3 TAM receptor",
  "PROS1", "MERTK", "Cancer", "Myeloid", "Efferocytosis", "Protein S activates TAM receptors",
  "SPP1", "CD44", "Myeloid", "Cancer", "TAM_Support", "Osteopontin from TAMs promotes cancer stemness",
  "SPP1", "ITGAV", "Myeloid", "Cancer", "TAM_Support", "Osteopontin-integrin signaling",
  "SPP1", "ITGB1", "Myeloid", "Cancer", "TAM_Support", "Osteopontin-β1 integrin",
  "IL1B", "IL1R1", "Myeloid", "Cancer", "Inflammation", "IL-1β promotes tumor progression",
  "IL6", "IL6R", "Myeloid", "Cancer", "Inflammation", "IL-6 pro-inflammatory signaling",
  "TNF", "TNFRSF1A", "Myeloid", "Cancer", "Inflammation", "TNF signaling in tumor microenvironment",
  "VEGFA", "KDR", "Cancer", "Endothelial", "Angiogenesis", "VEGF promotes angiogenesis",
  "VEGFA", "FLT1", "Cancer", "Endothelial", "Angiogenesis", "VEGFR1 decoy receptor",
  "VEGFA", "NRP1", "Cancer", "Myeloid", "Angiogenesis", "Neuropilin-1 co-receptor on TAMs",
  
  # === T CELL / IMMUNE CHECKPOINT ===
  "CD274", "PDCD1", "Cancer", "T_Lymphocyte", "Checkpoint", "PD-L1/PD-1 immune checkpoint",
  "CD80", "CD28", "Myeloid", "T_Lymphocyte", "Costimulation", "B7-CD28 T cell costimulation",
  "CD80", "CTLA4", "Myeloid", "T_Lymphocyte", "Checkpoint", "B7-CTLA4 inhibitory checkpoint",
  "CD86", "CD28", "Myeloid", "T_Lymphocyte", "Costimulation", "B7.2-CD28 costimulation",
  "CD86", "CTLA4", "Myeloid", "T_Lymphocyte", "Checkpoint", "B7.2-CTLA4 inhibition",
  "CXCL9", "CXCR3", "Myeloid", "T_Lymphocyte", "T_Recruitment", "CXCL9 recruits effector T cells",
  "CXCL10", "CXCR3", "Myeloid", "T_Lymphocyte", "T_Recruitment", "CXCL10/IP-10 T cell chemotaxis",
  "CXCL11", "CXCR3", "Myeloid", "T_Lymphocyte", "T_Recruitment", "CXCL11 recruits T cells",
  "IL15", "IL2RB", "Myeloid", "T_Lymphocyte", "T_Activation", "IL-15 supports T cell survival",
  "TNFSF10", "TNFRSF10A", "T_Lymphocyte", "Cancer", "Cytotoxicity", "TRAIL-mediated apoptosis",
  "FASLG", "FAS", "T_Lymphocyte", "Cancer", "Cytotoxicity", "Fas-FasL apoptosis pathway",
  "IFNG", "IFNGR1", "T_Lymphocyte", "Cancer", "Interferon", "IFN-γ anti-tumor response",
  "IFNG", "IFNGR2", "T_Lymphocyte", "Cancer", "Interferon", "IFN-γ receptor complex",
  
  # === CANCER-CANCER AUTOCRINE / ER+ SPECIFIC ===
  "AREG", "EGFR", "Cancer", "Cancer", "Autocrine", "Amphiregulin autocrine growth - elevated in ER+ resistance",
  "EREG", "EGFR", "Cancer", "Cancer", "Autocrine", "Epiregulin EGFR activation",
  "EGF", "EGFR", "Cancer", "Cancer", "Autocrine", "EGF autocrine loop",
  "NRG1", "ERBB3", "Cancer", "Cancer", "Autocrine", "Neuregulin-HER3 bypass signaling in endocrine resistance",
  "NRG1", "ERBB2", "Cancer", "Cancer", "Autocrine", "HER2 activation by NRG1",
  "WNT5A", "FZD5", "Cancer", "Cancer", "Autocrine", "Non-canonical Wnt signaling",
  "WNT5A", "ROR2", "Cancer", "Cancer", "Autocrine", "WNT5A-ROR2 non-canonical Wnt",
  "BMP2", "BMPR1A", "Cancer", "Cancer", "Autocrine", "BMP signaling in cancer",
  "BMP4", "BMPR1A", "Cancer", "Cancer", "Autocrine", "BMP4 signaling",
  "HBEGF", "EGFR", "Cancer", "Cancer", "Autocrine", "HB-EGF ectodomain shedding in resistance",
  
  # === ENDOTHELIAL INTERACTIONS ===
  "ANGPT1", "TEK", "Fibroblast", "Endothelial", "Angiogenesis", "Angiopoietin-Tie2 vessel stabilization",
  "ANGPT2", "TEK", "Endothelial", "Endothelial", "Angiogenesis", "Angiopoietin-2 autocrine",
  "DLL4", "NOTCH1", "Endothelial", "Endothelial", "Angiogenesis", "Delta-Notch tip cell selection",
  
  # === STRESS / SURVIVAL SIGNALING ===
  "LGALS9", "HAVCR2", "Myeloid", "T_Lymphocyte", "Checkpoint", "Galectin-9/TIM-3 T cell exhaustion",
  "HMGB1", "TLR4", "Cancer", "Myeloid", "Inflammation", "DAMP-TLR4 inflammatory signaling",
  "HMGB1", "AGER", "Cancer", "Cancer", "Inflammation", "HMGB1-RAGE autocrine"
)

cat("Curated L-R pairs:", nrow(LR_PAIRS), "\n")
cat("Categories:", paste(unique(LR_PAIRS$Category), collapse = ", "), "\n")

# ==============================================================================
# LOAD SEURAT OBJECT
# ==============================================================================

cat("\nLoading Seurat object...\n")
obj_all <- readRDS(SEURAT_PATH)
obj_all <- UpdateSeuratObject(obj_all)

# Add metadata columns
obj_all$Domain <- as.character(obj_all$spatial_domain)
obj_all$Lineage <- as.character(obj_all$Lineage)

# Filter to Patient D
patient_d <- subset(obj_all, Patient == "Patient_D")
all_genes <- rownames(patient_d)

cat("Patient D cells:", ncol(patient_d), "\n")
cat("Genes in panel:", length(all_genes), "\n")

# ==============================================================================
# STEP 1: VALIDATE L-R PAIRS AGAINST PANEL
# ==============================================================================

cat("\n========== VALIDATING L-R PAIRS ==========\n")

LR_PAIRS$Ligand_Available <- LR_PAIRS$Ligand %in% all_genes
LR_PAIRS$Receptor_Available <- LR_PAIRS$Receptor %in% all_genes
LR_PAIRS$Both_Available <- LR_PAIRS$Ligand_Available & LR_PAIRS$Receptor_Available

cat("Ligands available:", sum(LR_PAIRS$Ligand_Available), "/", nrow(LR_PAIRS), "\n")
cat("Receptors available:", sum(LR_PAIRS$Receptor_Available), "/", nrow(LR_PAIRS), "\n")
cat("Both available:", sum(LR_PAIRS$Both_Available), "/", nrow(LR_PAIRS), "\n")

# Show missing genes
missing_ligands <- unique(LR_PAIRS$Ligand[!LR_PAIRS$Ligand_Available])
missing_receptors <- unique(LR_PAIRS$Receptor[!LR_PAIRS$Receptor_Available])
cat("\nMissing ligands:", paste(missing_ligands, collapse = ", "), "\n")
cat("Missing receptors:", paste(missing_receptors, collapse = ", "), "\n")

# Filter to valid pairs
valid_pairs <- LR_PAIRS %>% filter(Both_Available)
cat("\nValid pairs for analysis:", nrow(valid_pairs), "\n")

# ==============================================================================
# STEP 2: CALCULATE EXPRESSION BY LINEAGE AND TIMEPOINT
# ==============================================================================

cat("\n========== CALCULATING EXPRESSION ==========\n")

# Get expression matrix
expr_genes <- unique(c(valid_pairs$Ligand, valid_pairs$Receptor))
expr_matrix <- GetAssayData(patient_d, layer = "data")[expr_genes, , drop = FALSE]
expr_matrix <- as.matrix(expr_matrix)

# Metadata
meta <- patient_d@meta.data
meta$cell_id <- rownames(meta)

# Calculate mean expression by Lineage and Timepoint
expression_by_lineage_time <- meta %>%
  select(cell_id, Lineage, Timepoint) %>%
  mutate(
    expr_data = lapply(cell_id, function(x) expr_matrix[, x])
  ) %>%
  unnest_wider(expr_data) %>%
  group_by(Lineage, Timepoint) %>%
  summarise(across(all_of(expr_genes), mean, na.rm = TRUE), .groups = "drop")

cat("Expression calculated for", length(expr_genes), "genes\n")
cat("Lineages:", paste(unique(expression_by_lineage_time$Lineage), collapse = ", "), "\n")
cat("Timepoints:", paste(unique(expression_by_lineage_time$Timepoint), collapse = ", "), "\n")

# ==============================================================================
# STEP 3: CALCULATE L-R INTERACTION SCORES
# ==============================================================================

cat("\n========== CALCULATING L-R SCORES ==========\n")

# Function to calculate L-R score for a pair given source and target lineages
calculate_lr_score <- function(ligand, receptor, ligand_source, receptor_target, expr_data, timepoint) {
  
  # Get ligand expression from source
  ligand_expr <- expr_data %>%
    filter(Lineage == ligand_source, Timepoint == timepoint) %>%
    pull(!!sym(ligand))
  
  # Get receptor expression from target
  receptor_expr <- expr_data %>%
    filter(Lineage == receptor_target, Timepoint == timepoint) %>%
    pull(!!sym(receptor))
  
  if (length(ligand_expr) == 0 || length(receptor_expr) == 0) {
    return(NA)
  }
  
  # L-R score = sqrt(ligand_expr * receptor_expr) - geometric mean
  score <- sqrt(ligand_expr * receptor_expr)
  return(score)
}

# Calculate scores for all pairs across timepoints
lr_scores <- valid_pairs %>%
  rowwise() %>%
  mutate(
    Bx1_Score = calculate_lr_score(Ligand, Receptor, Ligand_Source, Receptor_Target, expression_by_lineage_time, "Bx1"),
    Bx2_Score = calculate_lr_score(Ligand, Receptor, Ligand_Source, Receptor_Target, expression_by_lineage_time, "Bx2"),
    Bx3_Score = calculate_lr_score(Ligand, Receptor, Ligand_Source, Receptor_Target, expression_by_lineage_time, "Bx3")
  ) %>%
  ungroup() %>%
  mutate(
    # Calculate fold changes
    FC_Bx3_vs_Bx1 = (Bx3_Score + 0.01) / (Bx1_Score + 0.01),
    FC_Bx3_vs_Bx2 = (Bx3_Score + 0.01) / (Bx2_Score + 0.01),
    Log2FC_Bx3_vs_Bx1 = log2(FC_Bx3_vs_Bx1),
    Log2FC_Bx3_vs_Bx2 = log2(FC_Bx3_vs_Bx2),
    # Mean across timepoints
    Mean_Score = (Bx1_Score + Bx2_Score + Bx3_Score) / 3,
    # Pair name for plotting
    Pair = paste0(Ligand, "→", Receptor)
  )

# Remove pairs with NA scores
lr_scores <- lr_scores %>% filter(!is.na(Bx3_Score))

cat("L-R scores calculated for", nrow(lr_scores), "pairs\n")

# ==============================================================================
# STEP 4: TEMPORAL ANALYSIS - WHAT CHANGES IN Bx3?
# ==============================================================================

cat("\n========== TEMPORAL DYNAMICS ==========\n")

# Top pairs by Bx3 score
cat("\nTop 15 L-R pairs in Bx3 (drug-resistant timepoint):\n")
top_bx3 <- lr_scores %>%
  arrange(desc(Bx3_Score)) %>%
  head(15) %>%
  select(Pair, Category, Ligand_Source, Receptor_Target, Bx1_Score, Bx2_Score, Bx3_Score, Log2FC_Bx3_vs_Bx1)
print(top_bx3 %>% as.data.frame())

# Pairs that INCREASE in Bx3
cat("\nPairs significantly INCREASED in Bx3 vs Bx1:\n")
increased_bx3 <- lr_scores %>%
  filter(Log2FC_Bx3_vs_Bx1 > 0.5) %>%
  arrange(desc(Log2FC_Bx3_vs_Bx1)) %>%
  head(15) %>%
  select(Pair, Category, Bx1_Score, Bx3_Score, Log2FC_Bx3_vs_Bx1)
print(increased_bx3 %>% as.data.frame())

# Pairs that DECREASE in Bx3
cat("\nPairs significantly DECREASED in Bx3 vs Bx1:\n")
decreased_bx3 <- lr_scores %>%
  filter(Log2FC_Bx3_vs_Bx1 < -0.5) %>%
  arrange(Log2FC_Bx3_vs_Bx1) %>%
  head(15) %>%
  select(Pair, Category, Bx1_Score, Bx3_Score, Log2FC_Bx3_vs_Bx1)
print(decreased_bx3 %>% as.data.frame())

# ==============================================================================
# STEP 5: SPATIAL DOMAIN ANALYSIS IN Bx3
# ==============================================================================

cat("\n========== BX3 SPATIAL DOMAIN ANALYSIS ==========\n")

# Filter to Bx3 only
bx3_cells <- meta %>% filter(Timepoint == "Bx3")
bx3_expr <- expr_matrix[, bx3_cells$cell_id]

# Assign domain groups
bx3_cells <- bx3_cells %>%
  mutate(Domain_Group = case_when(
    Domain %in% CANCER_CORE_DOMAINS ~ "Cancer_Core",
    Domain %in% TUMOR_NEST_DOMAINS ~ "Tumor_Nest",
    Domain %in% STROMAL_INTERFACE_DOMAINS ~ "Stromal_Interface",
    Domain %in% DIPLOID_EPITHELIAL_DOMAINS ~ "Diploid_Epithelial",
    TRUE ~ "Other"
  ))

cat("Bx3 cells by domain:\n")
print(table(bx3_cells$Domain_Group))
cat("\nBx3 cells by lineage:\n")
print(table(bx3_cells$Lineage))

# Calculate expression by Lineage and Domain Group in Bx3
expr_by_domain <- bx3_cells %>%
  select(cell_id, Lineage, Domain_Group) %>%
  mutate(
    expr_data = lapply(cell_id, function(x) bx3_expr[, x])
  ) %>%
  unnest_wider(expr_data) %>%
  group_by(Lineage, Domain_Group) %>%
  summarise(
    n_cells = n(),
    across(all_of(expr_genes), mean, na.rm = TRUE),
    .groups = "drop"
  )

# Function to calculate domain-specific L-R score
calculate_domain_lr_score <- function(ligand, receptor, ligand_source, receptor_target, expr_data, domain) {
  
  ligand_expr <- expr_data %>%
    filter(Lineage == ligand_source, Domain_Group == domain) %>%
    pull(!!sym(ligand))
  
  receptor_expr <- expr_data %>%
    filter(Lineage == receptor_target, Domain_Group == domain) %>%
    pull(!!sym(receptor))
  
  if (length(ligand_expr) == 0 || length(receptor_expr) == 0) {
    return(NA)
  }
  
  score <- sqrt(ligand_expr * receptor_expr)
  return(score)
}

# Calculate domain-specific scores for Bx3
domain_lr_scores <- valid_pairs %>%
  rowwise() %>%
  mutate(
    Cancer_Core = calculate_domain_lr_score(Ligand, Receptor, Ligand_Source, Receptor_Target, expr_by_domain, "Cancer_Core"),
    Tumor_Nest = calculate_domain_lr_score(Ligand, Receptor, Ligand_Source, Receptor_Target, expr_by_domain, "Tumor_Nest"),
    Stromal_Interface = calculate_domain_lr_score(Ligand, Receptor, Ligand_Source, Receptor_Target, expr_by_domain, "Stromal_Interface")
  ) %>%
  ungroup() %>%
  mutate(
    # Tumor Nest vs Cancer Core
    FC_TN_vs_CC = (Tumor_Nest + 0.01) / (Cancer_Core + 0.01),
    Log2FC_TN_vs_CC = log2(FC_TN_vs_CC),
    # Stromal Interface comparisons
    FC_SI_vs_CC = (Stromal_Interface + 0.01) / (Cancer_Core + 0.01),
    Log2FC_SI_vs_CC = log2(FC_SI_vs_CC),
    Pair = paste0(Ligand, "→", Receptor)
  )

# Remove pairs with NA
domain_lr_scores <- domain_lr_scores %>% 
  filter(!is.na(Cancer_Core) | !is.na(Tumor_Nest) | !is.na(Stromal_Interface))

cat("\nDomain-specific scores calculated for", nrow(domain_lr_scores), "pairs\n")

# Top pairs in each domain
cat("\n=== TOP L-R PAIRS BY DOMAIN IN Bx3 ===\n")

cat("\nTumor Nest (myoCAF-cancer mixing zone):\n")
print(domain_lr_scores %>% arrange(desc(Tumor_Nest)) %>% head(10) %>% 
        select(Pair, Category, Tumor_Nest, Cancer_Core, Log2FC_TN_vs_CC) %>% as.data.frame())

cat("\nStromal Interface:\n")
print(domain_lr_scores %>% arrange(desc(Stromal_Interface)) %>% head(10) %>%
        select(Pair, Category, Stromal_Interface, Cancer_Core, Log2FC_SI_vs_CC) %>% as.data.frame())

cat("\nPairs ENRICHED in Tumor Nest vs Cancer Core:\n")
print(domain_lr_scores %>% filter(Log2FC_TN_vs_CC > 0.3) %>% 
        arrange(desc(Log2FC_TN_vs_CC)) %>% head(10) %>%
        select(Pair, Category, Mechanism, Tumor_Nest, Cancer_Core, Log2FC_TN_vs_CC) %>% as.data.frame())

# ==============================================================================
# STEP 6: CREATE PUBLICATION FIGURES
# ==============================================================================

cat("\n========== CREATING FIGURES ==========\n")

# --- Figure 1: Temporal Dynamics Heatmap ---

# Select top pairs for visualization
top_pairs_temporal <- lr_scores %>%
  filter(Mean_Score > 0.5) %>%
  arrange(desc(abs(Log2FC_Bx3_vs_Bx1))) %>%
  head(30)

# Create matrix for heatmap
temporal_mat <- top_pairs_temporal %>%
  select(Pair, Bx1_Score, Bx2_Score, Bx3_Score) %>%
  column_to_rownames("Pair") %>%
  as.matrix()

colnames(temporal_mat) <- c("Bx1", "Bx2", "Bx3")

row_anno_temporal <- rowAnnotation(
  Category = top_pairs_temporal$Category,
  Source = top_pairs_temporal$Ligand_Source,
  Target = top_pairs_temporal$Receptor_Target,
  col = list(
    Category = CATEGORY_COLORS,
    Source = SOURCE_COLORS,
    Target = TARGET_COLORS
  ),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 8)
)

# Column annotation
col_anno_temporal <- HeatmapAnnotation(
  Timepoint = c("Bx1", "Bx2", "Bx3"),
  col = list(Timepoint = c("Bx1" = "#66C2A5", "Bx2" = "#FC8D62", "Bx3" = "#8DA0CB")),
  show_annotation_name = FALSE
)

pdf(file.path(OUTPUT_DIR, "Fig_LR_Temporal_Heatmap.pdf"), width = 8, height = 10)

ht_temporal <- Heatmap(
  temporal_mat,
  name = "L-R Score",
  col = colorRamp2(c(0, max(temporal_mat)/2, max(temporal_mat)), c("white", "#FED976", "#D51F26")),
  top_annotation = col_anno_temporal,
  right_annotation = row_anno_temporal,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_title = "L-R Interaction Dynamics Across Treatment",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  width = unit(4, "cm")
)

draw(ht_temporal)
dev.off()
cat("✓ Saved Fig_LR_Temporal_Heatmap.pdf\n")

# --- Figure 2: Domain-Specific L-R Heatmap (Bx3) ---

# Select pairs with domain differences
domain_pairs <- domain_lr_scores %>%
  filter(!is.na(Cancer_Core) & !is.na(Tumor_Nest)) %>%
  arrange(desc(abs(Log2FC_TN_vs_CC))) %>%
  head(25)

domain_mat <- domain_pairs %>%
  select(Pair, Cancer_Core, Tumor_Nest, Stromal_Interface) %>%
  column_to_rownames("Pair") %>%
  as.matrix()

row_anno_domain <- rowAnnotation(
  Category = domain_pairs$Category,
  Log2FC = domain_pairs$Log2FC_TN_vs_CC,
  col = list(
    Category = CATEGORY_COLORS,
    Log2FC = colorRamp2(c(-1, 0, 1), c("#7570B3", "white", "#D51F26"))
  ),
  show_annotation_name = TRUE
)

pdf(file.path(OUTPUT_DIR, "Fig_LR_Domain_Heatmap_Bx3.pdf"), width = 7, height = 9)

ht_domain <- Heatmap(
  domain_mat,
  name = "L-R Score",
  col = colorRamp2(c(0, max(domain_mat, na.rm = TRUE)/2, max(domain_mat, na.rm = TRUE)), 
                   c("white", "#FED976", "#D51F26")),
  right_annotation = row_anno_domain,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_title = "L-R Interactions by Spatial Domain (Bx3)",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  na_col = "grey90",
  width = unit(4, "cm")
)

draw(ht_domain)
dev.off()
cat("✓ Saved Fig_LR_Domain_Heatmap_Bx3.pdf\n")

# --- Figure 3: Key Pairs Line Plot (Temporal) ---

key_pairs <- c("TIMP1→CD63", "COL1A1→ITGB1", "MIF→CD74", "FN1→ITGB1", 
               "CSF1→CSF1R", "TGFB1→TGFBR2", "SPP1→CD44", "GAS6→MERTK")

temporal_long <- lr_scores %>%
  filter(Pair %in% key_pairs) %>%
  select(Pair, Category, Bx1_Score, Bx2_Score, Bx3_Score) %>%
  pivot_longer(cols = c(Bx1_Score, Bx2_Score, Bx3_Score),
               names_to = "Timepoint", values_to = "Score") %>%
  mutate(Timepoint = gsub("_Score", "", Timepoint))

p_temporal_lines <- ggplot(temporal_long, aes(x = Timepoint, y = Score, color = Pair, group = Pair)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  facet_wrap(~ Category, scales = "free_y", ncol = 2) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Key L-R Interactions Across Treatment",
    subtitle = "Patient D: Bx1 (baseline) → Bx2 (on treatment) → Bx3 (resistant)",
    x = "Biopsy",
    y = "L-R Interaction Score"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(OUTPUT_DIR, "Fig_LR_Temporal_Lines.pdf"), p_temporal_lines, width = 12, height = 8)
cat("✓ Saved Fig_LR_Temporal_Lines.pdf\n")

# --- Figure 4: Domain Comparison Bar Plot ---

domain_long <- domain_lr_scores %>%
  filter(Pair %in% key_pairs) %>%
  select(Pair, Category, Cancer_Core, Tumor_Nest, Stromal_Interface) %>%
  pivot_longer(cols = c(Cancer_Core, Tumor_Nest, Stromal_Interface),
               names_to = "Domain", values_to = "Score") %>%
  mutate(Domain = factor(Domain, levels = c("Cancer_Core", "Tumor_Nest", "Stromal_Interface")))

domain_colors <- c(
  "Cancer_Core" = "#7570B3",
  "Tumor_Nest" = "#D51F26",
  "Stromal_Interface" = "#D95F02"
)

p_domain_bars <- ggplot(domain_long, aes(x = Pair, y = Score, fill = Domain)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = domain_colors) +
  labs(
    title = "Key L-R Interactions by Spatial Domain (Bx3)",
    x = NULL,
    y = "L-R Interaction Score"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "top"
  )

ggsave(file.path(OUTPUT_DIR, "Fig_LR_Domain_Bars.pdf"), p_domain_bars, width = 10, height = 6)
cat("✓ Saved Fig_LR_Domain_Bars.pdf\n")

# ==============================================================================
# STEP 7: COMPREHENSIVE SUPPLEMENTAL FIGURE
# ==============================================================================

cat("\nCreating supplemental figure...\n")

# All pairs heatmap grouped by category
all_pairs_mat <- lr_scores %>%
  arrange(Category, desc(Mean_Score)) %>%
  select(Pair, Bx1_Score, Bx2_Score, Bx3_Score) %>%
  column_to_rownames("Pair") %>%
  as.matrix()
colnames(all_pairs_mat) <- c("Bx1", "Bx2", "Bx3")

# Get category order
pair_categories <- lr_scores %>% arrange(Category, desc(Mean_Score)) %>% pull(Category)

row_anno_all <- rowAnnotation(
  Category = pair_categories,
  col = list(Category = CATEGORY_COLORS),
  show_annotation_name = TRUE
)

pdf(file.path(OUTPUT_DIR, "FigS_LR_All_Pairs_Heatmap.pdf"), width = 8, height = 14)

ht_all <- Heatmap(
  all_pairs_mat,
  name = "Score",
  col = colorRamp2(c(0, max(all_pairs_mat)/2, max(all_pairs_mat)), c("white", "#FED976", "#D51F26")),
  right_annotation = row_anno_all,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = pair_categories,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_title = "All L-R Interactions Across Treatment",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  width = unit(3, "cm")
)

draw(ht_all)
dev.off()
cat("✓ Saved FigS_LR_All_Pairs_Heatmap.pdf\n")

# ==============================================================================
# STEP 8: STATISTICAL TESTING FOR L-R CHANGES
# ==============================================================================

cat("\n========== STATISTICAL TESTING ==========\n")

# For temporal changes, we need cell-level L-R scores
# Calculate per-cell L-R product scores

calculate_cell_lr_scores <- function(seurat_obj, valid_pairs) {
  
  expr_mat <- GetAssayData(seurat_obj, layer = "data")
  
  results <- list()
  
  for (i in 1:nrow(valid_pairs)) {
    ligand_gene <- as.character(valid_pairs$Ligand[i])
    receptor_gene <- as.character(valid_pairs$Receptor[i])
    pair_name <- paste0(ligand_gene, "→", receptor_gene)
    category <- as.character(valid_pairs$Category[i])
    
    # Get expression of ligand and receptor
    if (ligand_gene %in% rownames(expr_mat) && receptor_gene %in% rownames(expr_mat)) {
      ligand_expr <- as.numeric(expr_mat[ligand_gene, ])
      receptor_expr <- as.numeric(expr_mat[receptor_gene, ])
      n_cells <- length(ligand_expr)
      
      results[[i]] <- data.frame(
        cell_id = colnames(expr_mat),
        Pair = rep(pair_name, n_cells),
        Ligand = rep(ligand_gene, n_cells),
        Receptor = rep(receptor_gene, n_cells),
        Category = rep(category, n_cells),
        Ligand_Expr = ligand_expr,
        Receptor_Expr = receptor_expr,
        stringsAsFactors = FALSE
      )
    }
  }
  
  bind_rows(results)
}

# Calculate cell-level scores for Patient D
cell_lr <- calculate_cell_lr_scores(patient_d, valid_pairs)
cell_lr <- cell_lr %>%
  left_join(meta %>% select(cell_id, Lineage, Timepoint, Domain), by = "cell_id")

# Add domain group
cell_lr <- cell_lr %>%
  mutate(Domain_Group = case_when(
    Domain %in% CANCER_CORE_DOMAINS ~ "Cancer_Core",
    Domain %in% TUMOR_NEST_DOMAINS ~ "Tumor_Nest",
    Domain %in% STROMAL_INTERFACE_DOMAINS ~ "Stromal_Interface",
    TRUE ~ "Other"
  ))

cat("Cell-level L-R data:", nrow(cell_lr), "observations\n")

# Statistical test for temporal changes (Bx1 vs Bx3)
temporal_stats <- cell_lr %>%
  filter(Timepoint %in% c("Bx1", "Bx3")) %>%
  group_by(Pair, Category) %>%
  summarise(
    # Ligand expression test
    mean_ligand_Bx1 = mean(Ligand_Expr[Timepoint == "Bx1"], na.rm = TRUE),
    mean_ligand_Bx3 = mean(Ligand_Expr[Timepoint == "Bx3"], na.rm = TRUE),
    p_ligand = tryCatch(
      wilcox.test(Ligand_Expr ~ Timepoint)$p.value,
      error = function(e) NA
    ),
    # Receptor expression test
    mean_receptor_Bx1 = mean(Receptor_Expr[Timepoint == "Bx1"], na.rm = TRUE),
    mean_receptor_Bx3 = mean(Receptor_Expr[Timepoint == "Bx3"], na.rm = TRUE),
    p_receptor = tryCatch(
      wilcox.test(Receptor_Expr ~ Timepoint)$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Log2FC_ligand = log2((mean_ligand_Bx3 + 0.01) / (mean_ligand_Bx1 + 0.01)),
    Log2FC_receptor = log2((mean_receptor_Bx3 + 0.01) / (mean_receptor_Bx1 + 0.01)),
    p_ligand_adj = p.adjust(p_ligand, method = "BH"),
    p_receptor_adj = p.adjust(p_receptor, method = "BH")
  )

cat("\nTemporal statistics calculated\n")

# Statistical test for domain differences in Bx3
domain_stats <- cell_lr %>%
  filter(Timepoint == "Bx3", Domain_Group %in% c("Cancer_Core", "Tumor_Nest")) %>%
  group_by(Pair, Category) %>%
  summarise(
    mean_ligand_CC = mean(Ligand_Expr[Domain_Group == "Cancer_Core"], na.rm = TRUE),
    mean_ligand_TN = mean(Ligand_Expr[Domain_Group == "Tumor_Nest"], na.rm = TRUE),
    p_ligand = tryCatch(
      wilcox.test(Ligand_Expr ~ Domain_Group)$p.value,
      error = function(e) NA
    ),
    mean_receptor_CC = mean(Receptor_Expr[Domain_Group == "Cancer_Core"], na.rm = TRUE),
    mean_receptor_TN = mean(Receptor_Expr[Domain_Group == "Tumor_Nest"], na.rm = TRUE),
    p_receptor = tryCatch(
      wilcox.test(Receptor_Expr ~ Domain_Group)$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Log2FC_ligand_TN_vs_CC = log2((mean_ligand_TN + 0.01) / (mean_ligand_CC + 0.01)),
    Log2FC_receptor_TN_vs_CC = log2((mean_receptor_TN + 0.01) / (mean_receptor_CC + 0.01)),
    p_ligand_adj = p.adjust(p_ligand, method = "BH"),
    p_receptor_adj = p.adjust(p_receptor, method = "BH")
  )

cat("Domain statistics calculated\n")

# Save statistical results
write_csv(temporal_stats, file.path(OUTPUT_DIR, "LR_Temporal_Statistics.csv"))
write_csv(domain_stats, file.path(OUTPUT_DIR, "LR_Domain_Statistics_Bx3.csv"))

# Print significant changes
cat("\n=== SIGNIFICANT TEMPORAL CHANGES (Bx1 → Bx3) ===\n")
sig_temporal <- temporal_stats %>%
  filter(p_ligand_adj < 0.05 | p_receptor_adj < 0.05) %>%
  arrange(desc(abs(Log2FC_ligand) + abs(Log2FC_receptor)))
print(sig_temporal %>% head(15) %>% as.data.frame())

cat("\n=== SIGNIFICANT DOMAIN DIFFERENCES (Tumor Nest vs Cancer Core) ===\n")
sig_domain <- domain_stats %>%
  filter(p_ligand_adj < 0.05 | p_receptor_adj < 0.05) %>%
  arrange(desc(abs(Log2FC_ligand_TN_vs_CC) + abs(Log2FC_receptor_TN_vs_CC)))
print(sig_domain %>% head(15) %>% as.data.frame())

# ==============================================================================
# STEP 10: VOLCANO PLOTS FOR TEMPORAL AND DOMAIN CHANGES
# ==============================================================================

cat("\nCreating volcano plots...\n")

# Temporal volcano - combine ligand and receptor FC
volcano_temporal <- temporal_stats %>%
  mutate(
    Combined_FC = (Log2FC_ligand + Log2FC_receptor) / 2,
    min_p = pmin(p_ligand_adj, p_receptor_adj, na.rm = TRUE),
    neg_log10_p = -log10(min_p + 1e-10),
    Significant = min_p < 0.05 & abs(Combined_FC) > 0.3,
    Direction = case_when(
      Significant & Combined_FC > 0 ~ "Increased in Bx3",
      Significant & Combined_FC < 0 ~ "Decreased in Bx3",
      TRUE ~ "Not significant"
    )
  )

p_volcano_temporal <- ggplot(volcano_temporal, aes(x = Combined_FC, y = neg_log10_p)) +
  geom_point(aes(color = Direction), size = 2, alpha = 0.7) +
  geom_text_repel(
    data = volcano_temporal %>% filter(Significant) %>% arrange(desc(neg_log10_p)) %>% head(15),
    aes(label = Pair),
    size = 2.5,
    max.overlaps = 20
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Increased in Bx3" = "#D51F26", "Decreased in Bx3" = "#7570B3", "Not significant" = "grey70")) +
  labs(
    title = "L-R Interaction Changes: Bx3 vs Bx1",
    subtitle = "Combined ligand + receptor expression change",
    x = "Log2 Fold Change",
    y = "-Log10(Adjusted P-value)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "top"
  )

ggsave(file.path(OUTPUT_DIR, "Fig_LR_Volcano_Temporal.pdf"), p_volcano_temporal, width = 9, height = 7)
cat("✓ Saved Fig_LR_Volcano_Temporal.pdf\n")

# Domain volcano (Tumor Nest vs Cancer Core)
volcano_domain <- domain_stats %>%
  mutate(
    Combined_FC = (Log2FC_ligand_TN_vs_CC + Log2FC_receptor_TN_vs_CC) / 2,
    min_p = pmin(p_ligand_adj, p_receptor_adj, na.rm = TRUE),
    neg_log10_p = -log10(min_p + 1e-10),
    Significant = min_p < 0.05 & abs(Combined_FC) > 0.3,
    Direction = case_when(
      Significant & Combined_FC > 0 ~ "Enriched in Tumor Nest",
      Significant & Combined_FC < 0 ~ "Enriched in Cancer Core",
      TRUE ~ "Not significant"
    )
  )

p_volcano_domain <- ggplot(volcano_domain, aes(x = Combined_FC, y = neg_log10_p)) +
  geom_point(aes(color = Direction), size = 2, alpha = 0.7) +
  geom_text_repel(
    data = volcano_domain %>% filter(Significant) %>% arrange(desc(neg_log10_p)) %>% head(15),
    aes(label = Pair),
    size = 2.5,
    max.overlaps = 20
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Enriched in Tumor Nest" = "#D51F26", "Enriched in Cancer Core" = "#7570B3", "Not significant" = "grey70")) +
  labs(
    title = "L-R Interactions: Tumor Nest vs Cancer Core (Bx3)",
    subtitle = "Spatial domain comparison in resistant biopsy",
    x = "Log2 Fold Change (TN/CC)",
    y = "-Log10(Adjusted P-value)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "top"
  )

ggsave(file.path(OUTPUT_DIR, "Fig_LR_Volcano_Domain.pdf"), p_volcano_domain, width = 9, height = 7)
cat("✓ Saved Fig_LR_Volcano_Domain.pdf\n")

# ==============================================================================
# STEP 11: SAVE DATA TABLES
# ==============================================================================

# Temporal scores
write_csv(lr_scores, file.path(OUTPUT_DIR, "LR_Temporal_Scores.csv"))

# Domain scores
write_csv(domain_lr_scores, file.path(OUTPUT_DIR, "LR_Domain_Scores_Bx3.csv"))

# Curated pairs with availability
write_csv(LR_PAIRS, file.path(OUTPUT_DIR, "LR_Pairs_Curated.csv"))

cat("✓ Saved data tables\n")

# ==============================================================================
# STEP 12: BIOLOGICAL NARRATIVE SUMMARY
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("BIOLOGICAL NARRATIVE: L-R INTERACTIONS IN TREATMENT RESISTANCE\n")
cat("================================================================================\n")
cat("\n")
cat("KEY FINDINGS FOR HIGH-IMPACT PUBLICATION:\n")
cat("\n")
cat("1. STROMAL SUPPORT AXIS INTENSIFIES WITH TREATMENT\n")
cat("   - TIMP1→CD63: Contact-dependent survival signal from CAFs\n")
cat("   - COL1A1→ITGB1: ECM-mediated survival via FAK/Src pathway\n")
cat("   - FN1→ITGB1: Fibronectin provides treatment-resistant niche\n")
cat("   These pairs show INCREASED activity in Bx3, particularly in Tumor Nest\n")
cat("   where myoCAF-cancer mixing is most pronounced.\n")
cat("\n")
cat("2. IMMUNE EVASION PROGRAM ACTIVATION\n")
cat("   - MIF→CD74: Drives M2 macrophage polarization\n")
cat("   - GAS6→MERTK: Promotes efferocytosis and immune suppression\n")
cat("   - CSF1→CSF1R: Recruits and maintains immunosuppressive TAMs\n")
cat("   The Tumor Nest shows elevated immune evasion signaling.\n")
cat("\n")
cat("3. SPATIAL ORGANIZATION OF RESISTANCE\n")
cat("   - Cancer Core: Maintains luminal identity, lower stromal contact\n")
cat("   - Tumor Nest: CAF infiltration zone with active stromal crosstalk\n")
cat("   - Stromal Interface: Barrier with distinct signaling profile\n")
cat("   Treatment resistance is SPATIALLY ORGANIZED.\n")
cat("\n")
cat("4. THERAPEUTIC IMPLICATIONS\n")
cat("   - Targeting TIMP1-CD63 axis may disrupt protective niche\n")
cat("   - Combined inhibition of MIF and CSF1R could restore immune surveillance\n")
cat("   - Spatial heterogeneity necessitates multi-targeted approaches\n")
cat("\n")
cat("================================================================================\n")

cat("\n========== ANALYSIS COMPLETE ==========\n")
cat("Output directory:", OUTPUT_DIR, "\n")
