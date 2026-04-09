#!/usr/bin/env python3
"""
Figure 7: Patient D - Curated Horizontal Bar Plot
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Paths
DATA_PATH = Path("supplementary_input_data/LR_Scores")
OUTPUT_PATH = Path("Figure7_Patient_Heatmaps")

# Load data
print("Loading data...")
scores_df = pd.read_csv(DATA_PATH / "LR_Scores_AllPatients_AllTimepoints.csv")
scores_df['Patient'] = scores_df['Patient'].str.replace('Patient_', '')

# Filter for Patient D, Bx1 and Bx3
patient_d = scores_df[(scores_df['Patient'] == 'D') &
                      (scores_df['Timepoint'].isin(['Bx1', 'Bx3']))].copy()

# Colors - Dark blue and Gold
COLOR_UP = '#D4A017'    # Gold for increase
COLOR_DOWN = '#1B4F72'  # Dark blue for decrease

# Curated pairs for Patient D
CURATED_PAIRS = {
    'Fibroblast→Cancer': ['TIMP1→CD63', 'VTN→ITGA3', 'COL4A2→ITGA1'],
    'Myeloid→Cancer': ['FN1→ITGA3', 'S100A9→ALCAM', 'THBS1→CD36'],
    'Cancer Autocrine': ['MIF→CD74', 'WNT7B→WIF1', 'IGF2→IGF2R'],
    'TCell→Cancer': ['LTB→LTBR', 'TNF→LTBR', 'FASLG→FAS'],
    'Endothelial→Cancer': ['EFNB2→EPHB3', 'CXCL12→ITGB1', 'JAG1→NOTCH1'],
}

KEY_CONTEXTS = list(CURATED_PAIRS.keys())

# Build curated dataframe
print("\nBuilding curated data...")
curated_data = []

for context, pairs in CURATED_PAIRS.items():
    for pair in pairs:
        pair_data = patient_d[patient_d['LR_Pair'] == pair]
        if len(pair_data) == 0:
            print(f"  Warning: {pair} not found")
            continue
        bx1 = pair_data[pair_data['Timepoint'] == 'Bx1']['Score'].values
        bx3 = pair_data[pair_data['Timepoint'] == 'Bx3']['Score'].values
        if len(bx1) == 0 or len(bx3) == 0:
            print(f"  Warning: {pair} missing timepoint data")
            continue
        fc = np.log2((bx3[0] + 0.001) / (bx1[0] + 0.001))
        curated_data.append({'Context': context, 'LR_Pair': pair, 'Log2FC': fc})
        print(f"  {pair}: Log2FC = {fc:+.2f}")

curated_df = pd.DataFrame(curated_data)

# Create figure
print("\nCreating curated bar plot...")
fig, ax = plt.subplots(figsize=(8, 10))

# Calculate positions with gaps between contexts
y_positions = []
y_labels = []
y = 0
context_mids = {}

for context in KEY_CONTEXTS:
    ctx_data = curated_df[curated_df['Context'] == context].sort_values('Log2FC', ascending=True)
    start_y = y
    for _, row in ctx_data.iterrows():
        y_positions.append(y)
        y_labels.append(row['LR_Pair'])
        y += 1
    if len(ctx_data) > 0:
        context_mids[context] = (start_y + y - 1) / 2
    y += 1.5  # Gap between contexts

# Plot bars
for i, (_, row) in enumerate(curated_df.iterrows()):
    idx = y_labels.index(row['LR_Pair'])
    y_pos = y_positions[idx]
    color = COLOR_UP if row['Log2FC'] > 0 else COLOR_DOWN
    ax.barh(y_pos, row['Log2FC'], height=0.7, color=color, edgecolor='white', linewidth=0.5)

# Y-axis labels
ax.set_yticks(y_positions)
ax.set_yticklabels(y_labels, fontsize=10)

# Context labels on right
for context, mid_y in context_mids.items():
    ax.text(ax.get_xlim()[1] + 0.3, mid_y, context, va='center', ha='left',
            fontsize=10, fontweight='bold')

# Vertical line at 0
ax.axvline(0, color='black', linewidth=1)

# Labels
ax.set_xlabel('Log2 Fold Change (Bx3 / Bx1)', fontsize=11)
ax.set_title('Patient D: Selected L-R Pairs\n(Bx1 → Bx3)', fontsize=13, fontweight='bold')

# Adjust x limits for labels
xlim = ax.get_xlim()
ax.set_xlim(xlim[0] - 0.5, xlim[1] + 0.5)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(OUTPUT_PATH / 'PatientD_Curated_BarPlot.pdf', bbox_inches='tight', dpi=300)
plt.close()
print("\nSaved: PatientD_Curated_BarPlot.pdf")

# Save summary
curated_df.to_csv(OUTPUT_PATH / 'PatientD_Curated_Pairs_Summary.csv', index=False)
print("Saved: PatientD_Curated_Pairs_Summary.csv")

print("\n" + "="*60)
print("COMPLETE!")
print("="*60)
