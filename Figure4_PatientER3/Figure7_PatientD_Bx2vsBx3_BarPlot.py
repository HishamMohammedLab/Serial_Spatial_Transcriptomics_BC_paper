#!/usr/bin/env python3
"""
Figure 7: Patient D - Curated Horizontal Bar Plot (Bx2 vs Bx3)
GATA3 gain-of-function, intratumor heterogeneity
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Paths
DATA_PATH = Path("supplementary_input_data/LR_Scores")
OUTPUT_PATH = Path("Figure7_Patient_Heatmaps")
OUTPUT_PATH.mkdir(parents=True, exist_ok=True)

# Publication styling
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['savefig.facecolor'] = 'white'

# Colors - Dark blue and Gold
COLOR_UP = '#D4A017'    # Gold for increase
COLOR_DOWN = '#1B4F72'  # Dark blue for decrease

# Load data
print("=" * 60)
print("Loading Data")
print("=" * 60)

scores_df = pd.read_csv(DATA_PATH / "LR_Scores_AllPatients_AllTimepoints.csv")
scores_df['Patient'] = scores_df['Patient'].str.replace('Patient_', '')

# Filter for Patient D, Bx2 and Bx3
patient_d = scores_df[(scores_df['Patient'] == 'D') &
                      (scores_df['Timepoint'].isin(['Bx2', 'Bx3']))].copy()

print(f"Patient D data: {len(patient_d)} rows")
print(f"Timepoints: {sorted(patient_d['Timepoint'].unique())}")

# Curated pairs for Patient D (from user specification)
CURATED_PAIRS = {
    'Fibroblast→Cancer': ['TIMP1→CD63', 'VTN→ITGA3', 'COL4A2→ITGA1', 'COL4A1→ITGA3', 'COL4A2→ITGA3'],
    'Myeloid→Cancer': ['FN1→ITGA3', 'S100A9→ALCAM', 'FN1→ITGB1', 'FN1→ITGAV', 'THBS1→CD36'],
    'Cancer Autocrine': ['MIF→CD74', 'WNT7B→WIF1', 'WNT5B→WIF1', 'WNT3→WIF1', 'IGF2→IGF2R'],
    'TCell→Cancer': ['LTB→LTBR', 'TNF→LTBR', 'IFNG→IFNGR1', 'FASLG→FAS', 'IL2→IL2RB', 'IL2→CD53'],
    'Endothelial→Cancer': ['EFNB2→EPHB3', 'CXCL12→ITGB1', 'JAG1→NOTCH1', 'FGF2→FGFR1'],
}

CONTEXT_ORDER = list(CURATED_PAIRS.keys())

# Build curated dataframe
print("\n" + "=" * 60)
print("Building Curated Data (Bx2 → Bx3)")
print("=" * 60)

curated_data = []

for context in CONTEXT_ORDER:
    pairs = CURATED_PAIRS[context]
    print(f"\n{context}:")

    for pair in pairs:
        pair_data = patient_d[patient_d['LR_Pair'] == pair]

        if len(pair_data) == 0:
            print(f"  {pair}: NOT FOUND")
            curated_data.append({
                'Context': context,
                'LR_Pair': pair,
                'Bx2_Score': np.nan,
                'Bx3_Score': np.nan,
                'Log2FC': np.nan
            })
            continue

        bx2 = pair_data[pair_data['Timepoint'] == 'Bx2']['Score'].values
        bx3 = pair_data[pair_data['Timepoint'] == 'Bx3']['Score'].values

        bx2_score = bx2[0] if len(bx2) > 0 else np.nan
        bx3_score = bx3[0] if len(bx3) > 0 else np.nan

        if pd.notna(bx2_score) and pd.notna(bx3_score):
            fc = np.log2((bx3_score + 0.001) / (bx2_score + 0.001))
            print(f"  {pair}: Log2FC = {fc:+.2f}")
        else:
            fc = np.nan
            print(f"  {pair}: MISSING DATA")

        curated_data.append({
            'Context': context,
            'LR_Pair': pair,
            'Bx2_Score': bx2_score,
            'Bx3_Score': bx3_score,
            'Log2FC': fc
        })

curated_df = pd.DataFrame(curated_data)

# Create figure
print("\n" + "=" * 60)
print("Creating Bar Plot")
print("=" * 60)

fig, ax = plt.subplots(figsize=(10, 12))

# Calculate positions with gaps between contexts
y_positions = []
y_labels = []
context_boundaries = []
y = 0

for context in CONTEXT_ORDER:
    ctx_data = curated_df[curated_df['Context'] == context]

    if len(ctx_data) == 0:
        continue

    context_start = y

    # Sort by FC (highest to lowest)
    ctx_data_sorted = ctx_data.sort_values('Log2FC', ascending=False, na_position='last')

    for _, row in ctx_data_sorted.iterrows():
        y_positions.append(y)
        y_labels.append(row['LR_Pair'])
        y += 1

    context_boundaries.append((context_start, context, len(ctx_data)))
    y += 1  # Gap between contexts

# Reorder curated_df to match y_labels order
curated_df_ordered = curated_df.set_index('LR_Pair').loc[y_labels].reset_index()

# Plot bars
for i, (_, row) in enumerate(curated_df_ordered.iterrows()):
    y_pos = y_positions[i]
    fc = row['Log2FC']

    if np.isnan(fc):
        color = '#CCCCCC'  # Gray for missing
        fc_plot = 0
    elif fc > 0:
        color = COLOR_UP
        fc_plot = fc
    else:
        color = COLOR_DOWN
        fc_plot = fc

    ax.barh(y_pos, fc_plot, height=0.75, color=color, edgecolor='white', linewidth=0.5)

# Add context labels (no separator lines)
for start_idx, context, n_pairs in context_boundaries:
    # Context label on the right side
    mid_y = start_idx + n_pairs / 2 - 0.5

    # Get x position for label
    xlim = ax.get_xlim()
    if xlim[1] == 0:
        # If no data yet, use default
        x_label = 2.5
    else:
        x_label = max(abs(xlim[0]), abs(xlim[1])) * 1.1

# Y-axis labels
ax.set_yticks(y_positions)
ax.set_yticklabels(y_labels, fontsize=10)

# Invert y-axis so first item is at top
ax.invert_yaxis()

# Vertical line at 0
ax.axvline(0, color='black', linewidth=1.2)

# Set x limits symmetrically
fc_values = curated_df_ordered['Log2FC'].dropna().values
if len(fc_values) > 0:
    fc_max = max(abs(fc_values.min()), abs(fc_values.max()))
    fc_lim = fc_max * 1.3
else:
    fc_lim = 2
ax.set_xlim(-fc_lim, fc_lim * 1.5)  # Extra space on right for labels

# Add context labels on right
for start_idx, context, n_pairs in context_boundaries:
    mid_y = start_idx + n_pairs / 2 - 0.5
    ax.text(fc_lim * 1.15, mid_y, context, va='center', ha='left',
            fontsize=11, fontweight='bold')

# Labels and title
ax.set_xlabel('Log₂ Fold Change (Bx3 / Bx2)', fontsize=12, fontweight='bold')
ax.set_title('Patient D: Curated L-R Pairs\nGATA3 gain-of-function (Bx2 → Bx3)',
             fontsize=14, fontweight='bold', pad=15)

# Style
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Legend
legend_elements = [
    Patch(facecolor=COLOR_UP, edgecolor='white', label='Increased (Bx3 > Bx2)'),
    Patch(facecolor=COLOR_DOWN, edgecolor='white', label='Decreased (Bx3 < Bx2)'),
]
ax.legend(handles=legend_elements, loc='lower right', fontsize=10, framealpha=0.9)

plt.tight_layout()

# Save (PDF only)
plt.savefig(OUTPUT_PATH / 'PatientD_Bx2vsBx3_Curated_BarPlot.pdf', bbox_inches='tight', dpi=300)
plt.close()

print(f"\nSaved: PatientD_Bx2vsBx3_Curated_BarPlot.pdf")

# Save summary table
curated_df.to_csv(OUTPUT_PATH / 'PatientD_Bx2vsBx3_Curated_Summary.csv', index=False)
print(f"Saved: PatientD_Bx2vsBx3_Curated_Summary.csv")

# Print summary
print("\n" + "=" * 60)
print("Summary")
print("=" * 60)
print(f"Total pairs: {len(curated_df)}")
print(f"Pairs with data: {curated_df['Log2FC'].notna().sum()}")
print(f"Increased (Bx3 > Bx2): {(curated_df['Log2FC'] > 0).sum()}")
print(f"Decreased (Bx3 < Bx2): {(curated_df['Log2FC'] < 0).sum()}")

print("\n" + "=" * 60)
print("COMPLETE!")
print("=" * 60)
