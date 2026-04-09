#!/usr/bin/env python3
"""
Publication-ready M1/M2 figures
1. Main: HLA-DR vs CD163 scatter (clean, no quadrant gating)
2. Supplemental: Box/violin plots for all M1/M2 markers
3. CD68 levels visualization
4. Macrophage gating rationale
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from pathlib import Path
import gzip
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================

BX2_PATH = Path("supplementary_input_data/CosmX protein\\Bx2_0000459948/0000459948_01112024_exprMat_file.csv.gz")
BX4_PATH = Path("supplementary_input_data/CosmX protein\\Bx4_0000459956/0000459956_01112024_exprMat_file.csv.gz")

OUTPUT_PATH = Path("M1_M2_Analysis")
OUTPUT_PATH.mkdir(parents=True, exist_ok=True)

# Publication-quality settings
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2

# Colors
BX2_COLOR = '#4DBBD5'  # Cyan/teal
BX4_COLOR = '#E64B35'  # Red/coral

# Gating thresholds (75th percentile based)
CD68_THRESHOLD = 84.0
EPCAM_THRESHOLD = 70.0
CD3_THRESHOLD = 17.0

# =============================================================================
# LOAD DATA
# =============================================================================
print("Loading data...")

def load_cosmx_data(path, timepoint):
    with gzip.open(path, 'rt') as f:
        df = pd.read_csv(f)
    df['Timepoint'] = timepoint
    return df

df_bx2 = load_cosmx_data(BX2_PATH, 'Bx2')
df_bx4 = load_cosmx_data(BX4_PATH, 'Bx4')
df = pd.concat([df_bx2, df_bx4], ignore_index=True)

print(f"Total cells: {len(df):,}")

# =============================================================================
# IDENTIFY MACROPHAGES
# =============================================================================
print("\nIdentifying macrophages...")

# Hierarchical gating (excluding tumor and T cells first)
tumor_mask = df['EpCAM'] > EPCAM_THRESHOLD
tcell_mask = (df['CD3'] > CD3_THRESHOLD) & (~tumor_mask)
mac_mask = (df['CD68'] > CD68_THRESHOLD) & (~tumor_mask) & (~tcell_mask)

macs = df[mac_mask].copy()
macs_bx2 = macs[macs['Timepoint'] == 'Bx2']
macs_bx4 = macs[macs['Timepoint'] == 'Bx4']

print(f"Macrophages: {len(macs):,} (Bx2: {len(macs_bx2):,}, Bx4: {len(macs_bx4):,})")

# =============================================================================
# FIGURE 1: MAIN - HLA-DR vs CD163 (Clean scatter, no quadrants)
# =============================================================================
print("\n" + "=" * 70)
print("Creating main figure: HLA-DR vs CD163 scatter")
print("=" * 70)

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

# Same axis limits
x_max = 180
y_max = 120

for idx, (data, tp, color) in enumerate([(macs_bx2, 'Bx2', BX2_COLOR),
                                          (macs_bx4, 'Bx4', BX4_COLOR)]):
    ax = axes[idx]

    # Scatter plot
    ax.scatter(data['CD163'], data['HLA-DR'], c=color, s=3, alpha=0.4,
               edgecolors='none', rasterized=True)

    ax.set_xlim(0, x_max)
    ax.set_ylim(0, y_max)
    ax.set_xlabel('CD163', fontsize=12, fontweight='bold')
    ax.set_ylabel('HLA-DR', fontsize=12, fontweight='bold')
    ax.set_title(f'{tp} (n = {len(data):,})', fontsize=13, fontweight='bold')

    # Add median lines (subtle)
    med_x = data['CD163'].median()
    med_y = data['HLA-DR'].median()
    ax.axvline(x=med_x, color='gray', linestyle=':', linewidth=0.8, alpha=0.5)
    ax.axhline(y=med_y, color='gray', linestyle=':', linewidth=0.8, alpha=0.5)

    # Add median values as text
    ax.text(0.98, 0.02, f'Median CD163: {med_x:.1f}\nMedian HLA-DR: {med_y:.1f}',
            transform=ax.transAxes, fontsize=9, ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8, edgecolor='none'))

plt.tight_layout()
plt.savefig(OUTPUT_PATH / 'Fig_Main_HLADR_CD163_Scatter.pdf', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: Fig_Main_HLADR_CD163_Scatter.pdf")

# =============================================================================
# FIGURE 2: SUPPLEMENTAL - Box/Violin plots for M1/M2 markers
# =============================================================================
print("\nCreating supplemental: M1/M2 marker box plots")

M1_MARKERS = ['HLA-DR', 'CD40', 'IL-1b', 'IL-18', 'iNOS']
M2_MARKERS = ['CD163', 'IDO1']
ALL_MARKERS = M1_MARKERS + M2_MARKERS

fig, axes = plt.subplots(2, 4, figsize=(14, 7))
axes = axes.flatten()

for idx, marker in enumerate(ALL_MARKERS):
    ax = axes[idx]

    data_bx2 = macs_bx2[marker].values
    data_bx4 = macs_bx4[marker].values

    # Violin plot
    parts = ax.violinplot([data_bx2, data_bx4], positions=[0, 1],
                          showmeans=False, showmedians=False, showextrema=False)

    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor([BX2_COLOR, BX4_COLOR][i])
        pc.set_alpha(0.6)
        pc.set_edgecolor('black')
        pc.set_linewidth(0.5)

    # Box plot overlay
    bp = ax.boxplot([data_bx2, data_bx4], positions=[0, 1], widths=0.15,
                    patch_artist=True, showfliers=False)

    for i, (box, median) in enumerate(zip(bp['boxes'], bp['medians'])):
        box.set_facecolor([BX2_COLOR, BX4_COLOR][i])
        box.set_alpha(0.9)
        box.set_edgecolor('black')
        median.set_color('black')
        median.set_linewidth(1.5)

    for whisker in bp['whiskers']:
        whisker.set_color('black')
        whisker.set_linewidth(1)
    for cap in bp['caps']:
        cap.set_color('black')
        cap.set_linewidth(1)

    # Statistics
    stat, pval = stats.mannwhitneyu(data_bx2, data_bx4, alternative='two-sided')

    # Significance annotation
    y_max = max(np.percentile(data_bx2, 95), np.percentile(data_bx4, 95))
    if pval < 0.001:
        sig_text = '***'
    elif pval < 0.01:
        sig_text = '**'
    elif pval < 0.05:
        sig_text = '*'
    else:
        sig_text = 'ns'

    ax.plot([0, 1], [y_max * 1.1, y_max * 1.1], 'k-', linewidth=1)
    ax.text(0.5, y_max * 1.15, sig_text, ha='center', va='bottom', fontsize=11, fontweight='bold')

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Bx2', 'Bx4'], fontsize=11)
    ax.set_ylabel('Expression', fontsize=11)

    # Label M1/M2
    marker_type = 'M1' if marker in M1_MARKERS else 'M2'
    ax.set_title(f'{marker} ({marker_type})', fontsize=12, fontweight='bold')

    ax.set_ylim(bottom=0)

# Hide empty subplot
axes[7].axis('off')

# Add legend
legend_elements = [mpatches.Patch(facecolor=BX2_COLOR, edgecolor='black', label='Bx2'),
                   mpatches.Patch(facecolor=BX4_COLOR, edgecolor='black', label='Bx4')]
axes[7].legend(handles=legend_elements, loc='center', fontsize=12, frameon=False)
axes[7].text(0.5, 0.3, '* p < 0.05\n** p < 0.01\n*** p < 0.001',
             ha='center', va='center', fontsize=10, transform=axes[7].transAxes)

plt.tight_layout()
plt.savefig(OUTPUT_PATH / 'Fig_Supp_M1M2_Markers_BoxPlots.pdf', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: Fig_Supp_M1M2_Markers_BoxPlots.pdf")

# =============================================================================
# FIGURE 3: CD68 Levels
# =============================================================================
print("\nCreating CD68 levels figure")

fig, axes = plt.subplots(1, 3, figsize=(12, 4))

# Panel A: CD68 distribution in all cells
ax = axes[0]
ax.hist(df[df['Timepoint'] == 'Bx2']['CD68'], bins=50, alpha=0.6, color=BX2_COLOR,
        label='Bx2', density=True, edgecolor='none')
ax.hist(df[df['Timepoint'] == 'Bx4']['CD68'], bins=50, alpha=0.6, color=BX4_COLOR,
        label='Bx4', density=True, edgecolor='none')
ax.axvline(x=CD68_THRESHOLD, color='black', linestyle='--', linewidth=1.5,
           label=f'Threshold ({CD68_THRESHOLD})')
ax.set_xlabel('CD68 Expression', fontsize=11, fontweight='bold')
ax.set_ylabel('Density', fontsize=11, fontweight='bold')
ax.set_title('CD68 Distribution\n(All Cells)', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)
ax.set_xlim(0, 400)

# Panel B: CD68 in macrophages only
ax = axes[1]
ax.hist(macs_bx2['CD68'], bins=40, alpha=0.6, color=BX2_COLOR,
        label=f'Bx2 (n={len(macs_bx2):,})', density=True, edgecolor='none')
ax.hist(macs_bx4['CD68'], bins=40, alpha=0.6, color=BX4_COLOR,
        label=f'Bx4 (n={len(macs_bx4):,})', density=True, edgecolor='none')
ax.set_xlabel('CD68 Expression', fontsize=11, fontweight='bold')
ax.set_ylabel('Density', fontsize=11, fontweight='bold')
ax.set_title('CD68 Distribution\n(Macrophages)', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)

# Panel C: CD68 box plot comparison
ax = axes[2]

data_bx2 = macs_bx2['CD68'].values
data_bx4 = macs_bx4['CD68'].values

parts = ax.violinplot([data_bx2, data_bx4], positions=[0, 1],
                      showmeans=False, showmedians=False, showextrema=False)
for i, pc in enumerate(parts['bodies']):
    pc.set_facecolor([BX2_COLOR, BX4_COLOR][i])
    pc.set_alpha(0.6)
    pc.set_edgecolor('black')

bp = ax.boxplot([data_bx2, data_bx4], positions=[0, 1], widths=0.15,
                patch_artist=True, showfliers=False)
for i, box in enumerate(bp['boxes']):
    box.set_facecolor([BX2_COLOR, BX4_COLOR][i])
    box.set_alpha(0.9)
for median in bp['medians']:
    median.set_color('black')
    median.set_linewidth(1.5)

# Stats
stat, pval = stats.mannwhitneyu(data_bx2, data_bx4)
sig_text = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'
y_max = np.percentile(data_bx4, 95)
ax.plot([0, 1], [y_max * 1.1, y_max * 1.1], 'k-', linewidth=1)
ax.text(0.5, y_max * 1.15, sig_text, ha='center', fontsize=12, fontweight='bold')

ax.set_xticks([0, 1])
ax.set_xticklabels(['Bx2', 'Bx4'], fontsize=11)
ax.set_ylabel('CD68 Expression', fontsize=11, fontweight='bold')
ax.set_title('CD68 Levels\n(Macrophages)', fontsize=12, fontweight='bold')

# Add mean values
ax.text(0, np.median(data_bx2) - 30, f'Med: {np.median(data_bx2):.0f}', ha='center', fontsize=9)
ax.text(1, np.median(data_bx4) - 30, f'Med: {np.median(data_bx4):.0f}', ha='center', fontsize=9)

plt.tight_layout()
plt.savefig(OUTPUT_PATH / 'Fig_Supp_CD68_Levels.pdf', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: Fig_Supp_CD68_Levels.pdf")

# =============================================================================
# FIGURE 4: Macrophage Gating Strategy (Rationale)
# =============================================================================
print("\nCreating macrophage gating rationale figure")

fig, axes = plt.subplots(2, 3, figsize=(13, 8))

# Panel A: CD68 vs EpCAM (showing tumor exclusion)
ax = axes[0, 0]
sample = df.sample(n=min(20000, len(df)), random_state=42)
ax.scatter(sample['EpCAM'], sample['CD68'], c='gray', s=1, alpha=0.2, rasterized=True)

# Highlight regions
ax.axhline(y=CD68_THRESHOLD, color='green', linestyle='--', linewidth=1.5, label=f'CD68 > {CD68_THRESHOLD}')
ax.axvline(x=EPCAM_THRESHOLD, color='red', linestyle='--', linewidth=1.5, label=f'EpCAM > {EPCAM_THRESHOLD}')

# Shade macrophage region
ax.fill_between([0, EPCAM_THRESHOLD], CD68_THRESHOLD, 600, alpha=0.2, color='green')
ax.fill_between([EPCAM_THRESHOLD, 300], 0, 600, alpha=0.2, color='red')

ax.set_xlabel('EpCAM (Tumor)', fontsize=11, fontweight='bold')
ax.set_ylabel('CD68 (Macrophage)', fontsize=11, fontweight='bold')
ax.set_title('Step 1: Exclude Tumor Cells', fontsize=12, fontweight='bold')
ax.set_xlim(0, 250)
ax.set_ylim(0, 500)
ax.legend(fontsize=8, loc='upper right')

# Panel B: CD68 vs CD3 (showing T cell exclusion)
ax = axes[0, 1]
# Filter out tumor first
non_tumor = sample[sample['EpCAM'] <= EPCAM_THRESHOLD]
ax.scatter(non_tumor['CD3'], non_tumor['CD68'], c='gray', s=1, alpha=0.3, rasterized=True)

ax.axhline(y=CD68_THRESHOLD, color='green', linestyle='--', linewidth=1.5)
ax.axvline(x=CD3_THRESHOLD, color='blue', linestyle='--', linewidth=1.5, label=f'CD3 > {CD3_THRESHOLD}')

# Shade macrophage region
ax.fill_between([0, CD3_THRESHOLD], CD68_THRESHOLD, 600, alpha=0.2, color='green')
ax.fill_between([CD3_THRESHOLD, 100], 0, 600, alpha=0.2, color='blue')

ax.set_xlabel('CD3 (T cell)', fontsize=11, fontweight='bold')
ax.set_ylabel('CD68 (Macrophage)', fontsize=11, fontweight='bold')
ax.set_title('Step 2: Exclude T Cells', fontsize=12, fontweight='bold')
ax.set_xlim(0, 80)
ax.set_ylim(0, 500)
ax.legend(fontsize=8, loc='upper right')

# Panel C: Final macrophage population on CD68
ax = axes[0, 2]
ax.hist(df['CD68'], bins=50, alpha=0.4, color='gray', label='All cells', density=True)
ax.hist(macs['CD68'], bins=50, alpha=0.7, color='green', label='Macrophages', density=True)
ax.axvline(x=CD68_THRESHOLD, color='black', linestyle='--', linewidth=1.5)
ax.set_xlabel('CD68 Expression', fontsize=11, fontweight='bold')
ax.set_ylabel('Density', fontsize=11, fontweight='bold')
ax.set_title('Final Macrophage Population', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)
ax.set_xlim(0, 500)

# Panel D: Gating summary table
ax = axes[1, 0]
ax.axis('off')

table_text = """
MACROPHAGE GATING CRITERIA
─────────────────────────────────────

1. CD68 > 84 (75th percentile)
   → Identifies myeloid cells

2. EpCAM ≤ 70 (exclude if > 75th %ile)
   → Removes tumor/epithelial cells

3. CD3 ≤ 17 (exclude if > 75th %ile)
   → Removes T cells

─────────────────────────────────────

Hierarchical gating applied in order:
  Tumor → T cells → Macrophages

This ensures mutually exclusive populations.
"""
ax.text(0.05, 0.95, table_text, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#f0f0f0', edgecolor='gray'))

# Panel E: Cell counts by gating step
ax = axes[1, 1]

# Calculate at each step
total = len(df)
after_tumor = (~tumor_mask).sum()
after_tcell = ((~tumor_mask) & (~tcell_mask)).sum()
final_macs = mac_mask.sum()

categories = ['All Cells', 'After Tumor\nExclusion', 'After T cell\nExclusion', 'CD68+\nMacrophages']
counts = [total, after_tumor, after_tcell, final_macs]
colors_bar = ['gray', '#ffcccc', '#cce5ff', 'green']

bars = ax.bar(categories, counts, color=colors_bar, edgecolor='black', linewidth=1)
ax.set_ylabel('Cell Count', fontsize=11, fontweight='bold')
ax.set_title('Gating Progression', fontsize=12, fontweight='bold')

# Add percentages
for bar, count in zip(bars, counts):
    pct = count / total * 100
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1000,
            f'{count:,}\n({pct:.1f}%)', ha='center', va='bottom', fontsize=9)

ax.set_ylim(0, total * 1.2)

# Panel F: Final composition by timepoint
ax = axes[1, 2]

bx2_mac_pct = len(macs_bx2) / len(df_bx2) * 100
bx4_mac_pct = len(macs_bx4) / len(df_bx4) * 100

x = [0, 1]
bars = ax.bar(x, [bx2_mac_pct, bx4_mac_pct], color=[BX2_COLOR, BX4_COLOR],
              edgecolor='black', linewidth=1)
ax.set_xticks(x)
ax.set_xticklabels(['Bx2', 'Bx4'], fontsize=11)
ax.set_ylabel('% Macrophages', fontsize=11, fontweight='bold')
ax.set_title('Macrophage Proportion', fontsize=12, fontweight='bold')

# Add counts
ax.text(0, bx2_mac_pct + 0.5, f'n={len(macs_bx2):,}\n({bx2_mac_pct:.1f}%)', ha='center', fontsize=10)
ax.text(1, bx4_mac_pct + 0.5, f'n={len(macs_bx4):,}\n({bx4_mac_pct:.1f}%)', ha='center', fontsize=10)

plt.tight_layout()
plt.savefig(OUTPUT_PATH / 'Fig_Supp_Macrophage_Gating_Rationale.pdf', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: Fig_Supp_Macrophage_Gating_Rationale.pdf")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================
print("\n" + "=" * 70)
print("SUMMARY STATISTICS")
print("=" * 70)

print(f"""
MACROPHAGE IDENTIFICATION:
  Criteria: CD68 > {CD68_THRESHOLD} AND EpCAM ≤ {EPCAM_THRESHOLD} AND CD3 ≤ {CD3_THRESHOLD}

  Total cells: {len(df):,}
  Macrophages: {len(macs):,} ({len(macs)/len(df)*100:.1f}%)
    - Bx2: {len(macs_bx2):,} ({len(macs_bx2)/len(df_bx2)*100:.1f}% of Bx2 cells)
    - Bx4: {len(macs_bx4):,} ({len(macs_bx4)/len(df_bx4)*100:.1f}% of Bx4 cells)

CD68 LEVELS IN MACROPHAGES:
  Bx2: Median = {macs_bx2['CD68'].median():.1f}, Mean = {macs_bx2['CD68'].mean():.1f}
  Bx4: Median = {macs_bx4['CD68'].median():.1f}, Mean = {macs_bx4['CD68'].mean():.1f}

HLA-DR (M1) IN MACROPHAGES:
  Bx2: Median = {macs_bx2['HLA-DR'].median():.1f}
  Bx4: Median = {macs_bx4['HLA-DR'].median():.1f}

CD163 (M2) IN MACROPHAGES:
  Bx2: Median = {macs_bx2['CD163'].median():.1f}
  Bx4: Median = {macs_bx4['CD163'].median():.1f}

OUTPUT FILES (in {OUTPUT_PATH}):
  - Fig_Main_HLADR_CD163_Scatter.pdf (Main figure)
  - Fig_Supp_M1M2_Markers_BoxPlots.pdf (Supplemental)
  - Fig_Supp_CD68_Levels.pdf (Supplemental)
  - Fig_Supp_Macrophage_Gating_Rationale.pdf (Supplemental)
""")

# Save statistics
stats_summary = {
    'Metric': ['Total cells', 'Macrophages', 'Bx2 macrophages', 'Bx4 macrophages',
               'Bx2 CD68 median', 'Bx4 CD68 median', 'Bx2 HLA-DR median', 'Bx4 HLA-DR median',
               'Bx2 CD163 median', 'Bx4 CD163 median'],
    'Value': [len(df), len(macs), len(macs_bx2), len(macs_bx4),
              macs_bx2['CD68'].median(), macs_bx4['CD68'].median(),
              macs_bx2['HLA-DR'].median(), macs_bx4['HLA-DR'].median(),
              macs_bx2['CD163'].median(), macs_bx4['CD163'].median()]
}
pd.DataFrame(stats_summary).to_csv(OUTPUT_PATH / 'Statistics_Summary.csv', index=False)
print("Saved: Statistics_Summary.csv")
