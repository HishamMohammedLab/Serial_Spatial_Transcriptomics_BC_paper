# L-R Pair Scoring Methodology

## Document Purpose
This file defines the standardized methodology for calculating and classifying L-R (Ligand-Receptor) pair changes. **All analyses must follow this methodology for consistency.**

---

## 1. Expression Metric: MEAN

We use **MEAN expression** per cell population (not median) because:
- Single-cell data has high dropout (many zeros)
- Median is often zero, making fold-change calculations meaningless
- Mean captures the overall "signaling potential" of a population
- Even if only 20% of cells express a gene, those cells can still signal

---

## 2. Lineage-Specific Context: YES

L-R pairs are analyzed in their **biologically relevant context**:
- **Ligand** expression measured in the **SENDER** cell type
- **Receptor** expression measured in the **RECEIVER** cell type

### Context Mapping (CosMx → FELINE)

| CosMx Context | FELINE Sender | FELINE Receiver |
|---------------|---------------|-----------------|
| Fibroblast_to_Cancer | CAF | Cancer |
| Myeloid_to_Cancer | Myeloid | Cancer |
| Cancer_to_Myeloid | Cancer | Myeloid |
| Cancer_Autocrine | Cancer | Cancer |
| TCell_to_Cancer | TCell | Cancer |
| Cancer_to_TCell | Cancer | TCell |

**Note:** Cancer_Autocrine = Cancer→Cancer (same cell type for both)

---

## 3. L-R Score Calculation: HYBRID METHOD

### Step 1: Calculate Individual Component FCs
```
Ligand_FC = log2((Ligand_Mean_Late + 0.01) / (Ligand_Mean_Early + 0.01))
Receptor_FC = log2((Receptor_Mean_Late + 0.01) / (Receptor_Mean_Early + 0.01))
```

### Step 2: Calculate Combined Score FC
```
Score_Early = (Ligand_Mean_Early + 0.01) × (Receptor_Mean_Early + 0.01)
Score_Late = (Ligand_Mean_Late + 0.01) × (Receptor_Mean_Late + 0.01)
Combined_FC = log2(Score_Late / Score_Early)
```

### Step 3: Concordance Check
A pair is classified based on whether ligand and receptor move in the same direction:

| Ligand_FC | Receptor_FC | Classification |
|-----------|-------------|----------------|
| > 0 | > 0 | **Concordant UP** |
| < 0 | < 0 | **Concordant DOWN** |
| > 0 | < 0 | **Discordant** (ligand up, receptor down) |
| < 0 | > 0 | **Discordant** (ligand down, receptor up) |

### Step 4: Final Classification (Per Patient)

For a pair to be classified as **"UP"** in a patient:
1. Combined_FC > 0, AND
2. **Neither component strongly decreases**: Both Ligand_FC > -0.5 AND Receptor_FC > -0.5

For a pair to be classified as **"DOWN"** in a patient:
1. Combined_FC < 0, AND
2. **Neither component strongly increases**: Both Ligand_FC < 0.5 AND Receptor_FC < 0.5

Otherwise: **"Variable/Discordant"**

---

## 4. Consistency Threshold: ≥50%

A pair is **"Consistent UP"** if ≥50% of patients show UP direction.
A pair is **"Consistent DOWN"** if ≤50% of patients show UP direction (i.e., >50% show DOWN).
Otherwise: **"Variable"** (only when exactly split or high discordance)

---

## 5. Resistance-Specific Classification

A pair is **Resistance-Specific** if:
- ≥66.7% of NON-RESPONDERS show UP, AND
- <50% of RESPONDERS show UP

---

## 6. Pseudocount

A pseudocount of **0.01** is added to all expression values to avoid log(0) errors.

---

## Version
- Created: 2024-01-24
- Method: Hybrid (Combined Score + Concordance Filter)
