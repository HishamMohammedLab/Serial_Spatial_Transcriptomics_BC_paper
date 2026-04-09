#!/usr/bin/env python3
"""
FINAL M1/M2 Analysis - Cleaned vs Original Comparison
- Uses GLOBAL thresholds (across all timepoints) for consistent classification
- Compares contaminated vs cleaned populations
- Shows dramatic differences in sample retention
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy import stats
from pathlib import Path
import gzip
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================

BX2_PATH = Path("supplementary_input_data/CosmX protein\\/Bx2_0000459948/0000459948_01112024_exprMat_file.csv.gz")
BX4_PATH = Path("supplementary_input_data/CosmX protein\\/Bx4_0000459956/0000459956_01112024_exprMat_file.csv.gz")

OUTPUT_PATH = Path("New Figures for paper/M1_M2_Analysis")
OUTPUT_PATH.mkdir(parents=True, exist_ok=True)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['figure.facecolor'] = 'white'

# =============================================================================
# LOAD DATA
# =============================================================================
print("=" * 70)
print("FINAL M1/M2 ANALYSIS - GLOBAL THRESHOLDS")
print("=" * 70)

print("\nLoading data...")

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
# GATING THRESHOLDS
# =============================================================================
print("\n" + "=" * 70)
print("Defining thresholds...")
print("=" * 70)

# CD68 gating
cd68_thresh = df['CD68'].quantile(0.75)
print(f"CD68 threshold (75th %ile): {cd68_thresh:.1f}")

# Contamination thresholds
epcam_thresh = df['EpCAM'].quantile(0.50)  # Median
cd3_thresh = df['CD3'].quantile(0.75)       # 75th percentile
print(f"EpCAM exclusion threshold (50th %ile): {epcam_thresh:.1f}")
print(f"CD3 exclusion threshold (75th %ile): {cd3_thresh:.1f}")

# =============================================================================
# CREATE POPULATIONS
# =============================================================================
print("\n" + "=" * 70)
print("Creating populations...")
print("=" * 70)

# Original CD68+ population
original_mask = df['CD68'] > cd68_thresh
original_pop = df[original_mask].copy()
print(f"\nOriginal CD68+: {len(original_pop):,}")

# Cleaned population (remove tumor and T cell contamination)
cleaned_mask = (df['CD68'] > cd68_thresh) & (df['EpCAM'] < epcam_thresh) & (df['CD3'] < cd3_thresh)
cleaned_pop = df[cleaned_mask].copy()
print(f"Cleaned CD68+: {len(cleaned_pop):,}")

# =============================================================================
# GLOBAL THRESHOLDS FOR M1/M2 CLASSIFICATION
# =============================================================================
print("\n" + "=" * 70)
print("Setting GLOBAL M1/M2 thresholds...")
print("=" * 70)

# Calculate thresholds on CLEANED population (better reference)
hladr_global = cleaned_pop['HLA-DR'].median()
cd163_global = cleaned_pop['CD163'].median()
print(f"HLA-DR median (global, cleaned): {hladr_global:.1f}")
print(f"CD163 median (global, cleaned): {cd163_global:.1f}")

# =============================================================================
# M1/M2 CLASSIFICATION FUNCTION
# =============================================================================

def classify_m1_m2(data, hladr_thresh, cd163_thresh):
    """Classify cells into M1, M2, Hybrid, Low using FIXED thresholds"""

    m1_mask = (data['HLA-DR'] > hladr_thresh) & (data['CD163'] <= cd163_thresh)
    m2_mask = (data['CD163'] > cd163_thresh) & (data['HLA-DR'] <= hladr_thresh)
    hybrid_mask = (data['HLA-DR'] > hladr_thresh) & (data['CD163'] > cd163_thresh)
    low_mask = (data['HLA-DR'] <= hladr_thresh) & (data['CD163'] <= cd163_thresh)

    return {
        'M1': m1_mask.sum(),
        'M2': m2_mask.sum(),
        'Hybrid': hybrid_mask.sum(),
        'Low': low_mask.sum(),
        'M1_pct': m1_mask.mean() * 100,
        'M2_pct': m2_mask.mean() * 100,
        'Hybrid_pct': hybrid_mask.mean() * 100,
        'Low_pct': low_mask.mean() * 100,
    }

# =============================================================================
# ANALYZE BOTH POPULATIONS
# =============================================================================
print("\n" + "=" * 70)
print("Running M1/M2 classification...")
print("=" * 70)

results = {}

for pop_name, pop in [('Original', original_pop), ('Cleaned', cleaned_pop)]:
    print(f"\n--- {pop_name} Population ---")
    results[pop_name] = {}

    for tp in ['Bx2', 'Bx4']:
        tp_data = pop[pop['Timepoint'] == tp]
        n = len(tp_data)

        if n < 10:
            print(f"  {tp}: n={n} (too few cells)")
            results[pop_name][tp] = {'n': n, 'M1_pct': np.nan, 'M2_pct': np.nan, 'M2_M1_ratio': np.nan}
            continue

        classification = classify_m1_m2(tp_data, hladr_global, cd163_global)
        classification['n'] = n
        classification['CD68_mean'] = tp_data['CD68'].mean()
        classification['EpCAM_mean'] = tp_data['EpCAM'].mean()
        classification['CD3_mean'] = tp_data['CD3'].mean()
        classification['M2_M1_ratio'] = classification['M2'] / classification['M1'] if classification['M1'] > 0 else np.nan

        results[pop_name][tp] = classification

        print(f"  {tp}: n={n:,}")
        print(f"    M1={classification['M1_pct']:.1f}%, M2={classification['M2_pct']:.1f}%, "
              f"Hybrid={classification['Hybrid_pct']:.1f}%, Low={classification['Low_pct']:.1f}%")
        print(f"    M2:M1 ratio = {classification['M2_M1_ratio']:.2f}")
        print(f"    CD68={classification['CD68_mean']:.1f}, EpCAM={classification['EpCAM_mean']:.1f}, CD3={classification['CD3_mean']:.1f}")

# =============================================================================
# COMPARISON SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("COMPARISON SUMMARY")
print("=" * 70)

print("\n            Original Population              Cleaned Population")
print("            ------------------              ------------------")
print("            Bx2          Bx4                Bx2          Bx4")
print("-" * 70)

# Sample sizes
orig_bx2_n = results['Original']['Bx2']['n']
orig_bx4_n = results['Original']['Bx4']['n']
clean_bx2_n = results['Cleaned']['Bx2']['n']
clean_bx4_n = results['Cleaned']['Bx4']['n']
print(f"N cells     {orig_bx2_n:>6,}     {orig_bx4_n:>6,}               {clean_bx2_n:>6,}     {clean_bx4_n:>6,}")

# Retention
bx2_retained = clean_bx2_n / orig_bx2_n * 100 if orig_bx2_n > 0 else 0
bx4_retained = clean_bx4_n / orig_bx4_n * 100 if orig_bx4_n > 0 else 0
print(f"% retained  {'100.0%':>10} {'100.0%':>10}              {bx2_retained:>6.1f}%    {bx4_retained:>6.1f}%")

# M1 percentages
print(f"\nM1 %        {results['Original']['Bx2']['M1_pct']:>10.1f} {results['Original']['Bx4']['M1_pct']:>10.1f}              {results['Cleaned']['Bx2']['M1_pct']:>6.1f}%    {results['Cleaned']['Bx4']['M1_pct']:>6.1f}%")
print(f"M2 %        {results['Original']['Bx2']['M2_pct']:>10.1f} {results['Original']['Bx4']['M2_pct']:>10.1f}              {results['Cleaned']['Bx2']['M2_pct']:>6.1f}%    {results['Cleaned']['Bx4']['M2_pct']:>6.1f}%")
print(f"Hybrid %    {results['Original']['Bx2']['Hybrid_pct']:>10.1f} {results['Original']['Bx4']['Hybrid_pct']:>10.1f}              {results['Cleaned']['Bx2']['Hybrid_pct']:>6.1f}%    {results['Cleaned']['Bx4']['Hybrid_pct']:>6.1f}%")

# M2:M1 ratios
orig_bx2_ratio = results['Original']['Bx2']['M2_M1_ratio']
orig_bx4_ratio = results['Original']['Bx4']['M2_M1_ratio']
clean_bx2_ratio = results['Cleaned']['Bx2']['M2_M1_ratio']
clean_bx4_ratio = results['Cleaned']['Bx4']['M2_M1_ratio']
print(f"\nM2:M1 ratio {orig_bx2_ratio:>10.2f} {orig_bx4_ratio:>10.2f}              {clean_bx2_ratio:>6.2f}      {clean_bx4_ratio:>6.2f}")

# Change
orig_change = (orig_bx4_ratio - orig_bx2_ratio) / orig_bx2_ratio * 100 if orig_bx2_ratio > 0 else np.nan
clean_change = (clean_bx4_ratio - clean_bx2_ratio) / clean_bx2_ratio * 100 if clean_bx2_ratio > 0 else np.nan
print(f"Change      {orig_change:>+10.0f}%                         {clean_change:>+6.0f}%")

# Contamination markers
print(f"\nEpCAM mean  {results['Original']['Bx2']['EpCAM_mean']:>10.1f} {results['Original']['Bx4']['EpCAM_mean']:>10.1f}              {results['Cleaned']['Bx2']['EpCAM_mean']:>6.1f}      {results['Cleaned']['Bx4']['EpCAM_mean']:>6.1f}")
print(f"CD3 mean    {results['Original']['Bx2']['CD3_mean']:>10.1f} {results['Original']['Bx4']['CD3_mean']:>10.1f}              {results['Cleaned']['Bx2']['CD3_mean']:>6.1f}      {results['Cleaned']['Bx4']['CD3_mean']:>6.1f}")

# =============================================================================
# VISUALIZATION 1: Side-by-side scatter plots
# =============================================================================
print("\n" + "=" * 70)
print("Creating visualizations...")
print("=" * 70)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Shared axis limits
x_max = max(original_pop['CD163'].quantile(0.99), 200)
y_max = max(original_pop['HLA-DR'].quantile(0.99), 100)

quadrant_colors = {
    'M1': '#E64B35',      # Red
    'M2': '#3C5488',      # Blue
    'Hybrid': '#9467BD',  # Purple
    'Low': '#CCCCCC'      # Gray
}

for col, (pop_name, pop) in enumerate([('Original', original_pop), ('Cleaned', cleaned_pop)]):
    for row, tp in enumerate(['Bx2', 'Bx4']):
        ax = axes[row, col]
        tp_data = pop[pop['Timepoint'] == tp]

        # Classify each cell
        m1_mask = (tp_data['HLA-DR'] > hladr_global) & (tp_data['CD163'] <= cd163_global)
        m2_mask = (tp_data['CD163'] > cd163_global) & (tp_data['HLA-DR'] <= hladr_global)
        hybrid_mask = (tp_data['HLA-DR'] > hladr_global) & (tp_data['CD163'] > cd163_global)
        low_mask = (tp_data['HLA-DR'] <= hladr_global) & (tp_data['CD163'] <= cd163_global)

        # Plot by category
        for mask, label, color in [(low_mask, 'Low', quadrant_colors['Low']),
                                    (hybrid_mask, 'Hybrid', quadrant_colors['Hybrid']),
                                    (m1_mask, 'M1', quadrant_colors['M1']),
                                    (m2_mask, 'M2', quadrant_colors['M2'])]:
            subset = tp_data[mask]
            ax.scatter(subset['CD163'], subset['HLA-DR'], c=color, s=2, alpha=0.4,
                      label=f'{label} ({mask.sum():,}, {mask.mean()*100:.1f}%)')

        # Threshold lines
        ax.axvline(x=cd163_global, color='black', linestyle='--', linewidth=1, alpha=0.7)
        ax.axhline(y=hladr_global, color='black', linestyle='--', linewidth=1, alpha=0.7)

        ax.set_xlim(0, x_max)
        ax.set_ylim(0, y_max)
        ax.set_xlabel('CD163 (M2 marker)', fontsize=10)
        ax.set_ylabel('HLA-DR (M1 marker)', fontsize=10)

        ratio = results[pop_name][tp]['M2_M1_ratio']
        ax.set_title(f'{pop_name} - {tp} (n={len(tp_data):,}, M2:M1={ratio:.2f})',
                    fontsize=11, fontweight='bold')
        ax.legend(loc='upper right', fontsize=7, markerscale=3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

plt.suptitle('M1/M2 Classification: Original vs Cleaned Populations (Global Thresholds)',
             fontsize=13, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig(OUTPUT_PATH / 'PatientA_M1M2_CleanedFinal_Scatter.pdf', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: PatientA_M1M2_CleanedFinal_Scatter.pdf")

# =============================================================================
# VISUALIZATION 2: Bar plot comparison
# =============================================================================
fig, axes = plt.subplots(1, 3, figsize=(14, 5))

# Panel A: Sample retention
ax = axes[0]
x = np.arange(2)
width = 0.35

orig_n = [results['Original']['Bx2']['n'], results['Original']['Bx4']['n']]
clean_n = [results['Cleaned']['Bx2']['n'], results['Cleaned']['Bx4']['n']]

bars1 = ax.bar(x - width/2, orig_n, width, label='Original', color='#888888', edgecolor='black')
bars2 = ax.bar(x + width/2, clean_n, width, label='Cleaned', color='#4DAF4A', edgecolor='black')

ax.set_xticks(x)
ax.set_xticklabels(['Bx2', 'Bx4'])
ax.set_ylabel('Cell Count', fontsize=11)
ax.set_title('Sample Size After Cleaning', fontsize=12, fontweight='bold')
ax.legend()

# Add retention percentages
for i, (o, c) in enumerate(zip(orig_n, clean_n)):
    pct = c / o * 100 if o > 0 else 0
    ax.text(i + width/2, c + max(orig_n)*0.02, f'{pct:.1f}%', ha='center', fontsize=9)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel B: M2:M1 ratio comparison
ax = axes[1]

ratios_orig = [results['Original']['Bx2']['M2_M1_ratio'], results['Original']['Bx4']['M2_M1_ratio']]
ratios_clean = [results['Cleaned']['Bx2']['M2_M1_ratio'], results['Cleaned']['Bx4']['M2_M1_ratio']]

bars1 = ax.bar(x - width/2, ratios_orig, width, label='Original', color='#888888', edgecolor='black')
bars2 = ax.bar(x + width/2, ratios_clean, width, label='Cleaned', color='#4DAF4A', edgecolor='black')

ax.set_xticks(x)
ax.set_xticklabels(['Bx2', 'Bx4'])
ax.set_ylabel('M2:M1 Ratio', fontsize=11)
ax.set_title('M2:M1 Ratio Comparison', fontsize=12, fontweight='bold')
ax.axhline(y=1.0, color='black', linestyle='--', alpha=0.5)
ax.legend()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel C: EpCAM contamination
ax = axes[2]

epcam_orig = [results['Original']['Bx2']['EpCAM_mean'], results['Original']['Bx4']['EpCAM_mean']]
epcam_clean = [results['Cleaned']['Bx2']['EpCAM_mean'], results['Cleaned']['Bx4']['EpCAM_mean']]

bars1 = ax.bar(x - width/2, epcam_orig, width, label='Original', color='#888888', edgecolor='black')
bars2 = ax.bar(x + width/2, epcam_clean, width, label='Cleaned', color='#4DAF4A', edgecolor='black')

ax.set_xticks(x)
ax.set_xticklabels(['Bx2', 'Bx4'])
ax.set_ylabel('EpCAM Mean Expression', fontsize=11)
ax.set_title('Tumor Contamination (EpCAM)', fontsize=12, fontweight='bold')
ax.legend()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(OUTPUT_PATH / 'PatientA_M1M2_CleanedFinal_BarPlots.pdf', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: PatientA_M1M2_CleanedFinal_BarPlots.pdf")

# =============================================================================
# VISUALIZATION 3: Stacked bar of M1/M2/Hybrid/Low composition
# =============================================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

categories = ['M1', 'M2', 'Hybrid', 'Low']
cat_colors = [quadrant_colors['M1'], quadrant_colors['M2'], quadrant_colors['Hybrid'], quadrant_colors['Low']]

for idx, (pop_name, pop) in enumerate([('Original', original_pop), ('Cleaned', cleaned_pop)]):
    ax = axes[idx]

    bx2_data = pop[pop['Timepoint'] == 'Bx2']
    bx4_data = pop[pop['Timepoint'] == 'Bx4']

    # Classify
    bx2_class = classify_m1_m2(bx2_data, hladr_global, cd163_global)
    bx4_class = classify_m1_m2(bx4_data, hladr_global, cd163_global)

    # Stacked bar
    bx2_vals = [bx2_class[f'{c}_pct'] for c in categories]
    bx4_vals = [bx4_class[f'{c}_pct'] for c in categories]

    x = np.arange(2)
    width = 0.5

    bottom_bx2 = 0
    bottom_bx4 = 0

    for i, cat in enumerate(categories):
        ax.bar(0, bx2_vals[i], width, bottom=bottom_bx2, color=cat_colors[i], label=cat if idx == 0 else None)
        ax.bar(1, bx4_vals[i], width, bottom=bottom_bx4, color=cat_colors[i])
        bottom_bx2 += bx2_vals[i]
        bottom_bx4 += bx4_vals[i]

    ax.set_xticks([0, 1])
    ax.set_xticklabels([f'Bx2\n(n={len(bx2_data):,})', f'Bx4\n(n={len(bx4_data):,})'])
    ax.set_ylabel('% of Population', fontsize=11)
    ax.set_title(f'{pop_name} Population', fontsize=12, fontweight='bold')
    ax.set_ylim(0, 100)

    if idx == 0:
        ax.legend(loc='upper right', fontsize=9)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add M2:M1 ratio annotation
    bx2_ratio = results[pop_name]['Bx2']['M2_M1_ratio']
    bx4_ratio = results[pop_name]['Bx4']['M2_M1_ratio']
    ax.text(0, 102, f'M2:M1={bx2_ratio:.2f}', ha='center', fontsize=9)
    ax.text(1, 102, f'M2:M1={bx4_ratio:.2f}', ha='center', fontsize=9)

plt.suptitle('M1/M2 Composition: Original vs Cleaned',
             fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(OUTPUT_PATH / 'PatientA_M1M2_CleanedFinal_Composition.pdf', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: PatientA_M1M2_CleanedFinal_Composition.pdf")

# =============================================================================
# FINAL SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("FINAL CONCLUSIONS")
print("=" * 70)

print(f"""
CONTAMINATION ASSESSMENT:

1. Bx2 MASSIVE CONTAMINATION:
   - Original CD68+ cells: {orig_bx2_n:,}
   - After cleaning: {clean_bx2_n:,} ({bx2_retained:.1f}% retained)
   - {100-bx2_retained:.1f}% of Bx2 "macrophages" were likely tumor/T cells
   - Mean EpCAM: {results['Original']['Bx2']['EpCAM_mean']:.1f} → {results['Cleaned']['Bx2']['EpCAM_mean']:.1f}

2. Bx4 MODERATE CONTAMINATION:
   - Original CD68+ cells: {orig_bx4_n:,}
   - After cleaning: {clean_bx4_n:,} ({bx4_retained:.1f}% retained)
   - {100-bx4_retained:.1f}% of Bx4 "macrophages" were likely tumor/T cells
   - Mean EpCAM: {results['Original']['Bx4']['EpCAM_mean']:.1f} → {results['Cleaned']['Bx4']['EpCAM_mean']:.1f}

M1/M2 ANALYSIS:

ORIGINAL (CONTAMINATED):
   Bx2: M2:M1 = {orig_bx2_ratio:.2f}
   Bx4: M2:M1 = {orig_bx4_ratio:.2f}
   Change: {orig_change:+.0f}%

CLEANED:
   Bx2: M2:M1 = {clean_bx2_ratio:.2f}
   Bx4: M2:M1 = {clean_bx4_ratio:.2f}
   Change: {clean_change:+.0f}%

INTERPRETATION:
- The M2 polarization shift persists in BOTH populations
- Original data shows {orig_change:+.0f}% increase in M2:M1 ratio
- Cleaned data shows {clean_change:+.0f}% increase in M2:M1 ratio
- However, Bx2 cleaned sample is VERY small (n={clean_bx2_n}) due to contamination
- Confidence in Bx2 data is limited

OUTPUT FILES:
- PatientA_M1M2_CleanedFinal_Scatter.pdf
- PatientA_M1M2_CleanedFinal_BarPlots.pdf
- PatientA_M1M2_CleanedFinal_Composition.pdf
""")

# Save summary table
summary_table = pd.DataFrame({
    'Population': ['Original', 'Original', 'Cleaned', 'Cleaned'],
    'Timepoint': ['Bx2', 'Bx4', 'Bx2', 'Bx4'],
    'N_cells': [orig_bx2_n, orig_bx4_n, clean_bx2_n, clean_bx4_n],
    'M1_pct': [results['Original']['Bx2']['M1_pct'], results['Original']['Bx4']['M1_pct'],
               results['Cleaned']['Bx2']['M1_pct'], results['Cleaned']['Bx4']['M1_pct']],
    'M2_pct': [results['Original']['Bx2']['M2_pct'], results['Original']['Bx4']['M2_pct'],
               results['Cleaned']['Bx2']['M2_pct'], results['Cleaned']['Bx4']['M2_pct']],
    'M2_M1_ratio': [orig_bx2_ratio, orig_bx4_ratio, clean_bx2_ratio, clean_bx4_ratio],
    'EpCAM_mean': [results['Original']['Bx2']['EpCAM_mean'], results['Original']['Bx4']['EpCAM_mean'],
                   results['Cleaned']['Bx2']['EpCAM_mean'], results['Cleaned']['Bx4']['EpCAM_mean']],
    'CD68_mean': [results['Original']['Bx2']['CD68_mean'], results['Original']['Bx4']['CD68_mean'],
                  results['Cleaned']['Bx2']['CD68_mean'], results['Cleaned']['Bx4']['CD68_mean']],
})
summary_table.to_csv(OUTPUT_PATH / 'PatientA_M1M2_CleanedFinal_Summary.csv', index=False)
print("Saved: PatientA_M1M2_CleanedFinal_Summary.csv")
