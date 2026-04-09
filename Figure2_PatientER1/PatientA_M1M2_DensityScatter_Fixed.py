#!/usr/bin/env python3
"""
Density scatter plots with fixed colorbar positioning (no overlap)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from mpl_toolkits.axes_grid1 import make_axes_locatable
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

# Publication settings
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2

# Gating thresholds
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

# Identify macrophages
tumor_mask = df['EpCAM'] > EPCAM_THRESHOLD
tcell_mask = (df['CD3'] > CD3_THRESHOLD) & (~tumor_mask)
mac_mask = (df['CD68'] > CD68_THRESHOLD) & (~tumor_mask) & (~tcell_mask)

macs = df[mac_mask].copy()
macs_bx2 = macs[macs['Timepoint'] == 'Bx2']
macs_bx4 = macs[macs['Timepoint'] == 'Bx4']

print(f"Macrophages: Bx2 = {len(macs_bx2):,}, Bx4 = {len(macs_bx4):,}")

# =============================================================================
# DENSITY SCATTER FUNCTION
# =============================================================================

def density_scatter(ax, x, y, cmap='viridis', s=4, alpha=0.8):
    """Create a density-colored scatter plot"""
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    if len(x) < 10:
        ax.scatter(x, y, c='steelblue', s=s, alpha=alpha)
        return None

    try:
        xy = np.vstack([x, y])
        density = gaussian_kde(xy)(xy)
        idx = density.argsort()
        x, y, density = x[idx], y[idx], density[idx]

        scatter = ax.scatter(x, y, c=density, cmap=cmap, s=s, alpha=alpha,
                            edgecolors='none', rasterized=True)
        return scatter
    except:
        ax.scatter(x, y, c='steelblue', s=s, alpha=alpha)
        return None

# =============================================================================
# FIGURE 1: MAIN - CD163 vs HLA-DR (Fixed colorbar)
# =============================================================================
print("\nCreating main figure: CD163 vs HLA-DR density scatter...")

# Axis limits based on Bx4
x_max = macs_bx4['CD163'].quantile(0.99)
y_max = macs_bx4['HLA-DR'].quantile(0.99)
x_max = np.ceil(x_max / 20) * 20
y_max = np.ceil(y_max / 20) * 20

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

scatter_obj = None
for idx, (data, tp) in enumerate([(macs_bx2, 'Bx2'), (macs_bx4, 'Bx4')]):
    ax = axes[idx]

    x = data['CD163'].values
    y = data['HLA-DR'].values

    scatter = density_scatter(ax, x, y, cmap='viridis', s=4, alpha=0.8)
    if scatter is not None:
        scatter_obj = scatter

    ax.set_xlim(0, x_max)
    ax.set_ylim(0, y_max)
    ax.set_xlabel('CD163', fontsize=13, fontweight='bold')
    ax.set_ylabel('HLA-DR', fontsize=13, fontweight='bold')
    ax.set_title(f'{tp} (n = {len(data):,})', fontsize=14, fontweight='bold')
    ax.grid(False)

# Add colorbar to the right, outside the plots
if scatter_obj is not None:
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(scatter_obj, cax=cbar_ax)
    cbar.set_label('Density', fontsize=12, fontweight='bold')
    cbar.ax.tick_params(labelsize=10)

plt.savefig(OUTPUT_PATH / 'Fig_Main_CD163_HLADR_Density.pdf', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: Fig_Main_CD163_HLADR_Density.pdf")

# =============================================================================
# SUPPLEMENTAL: Individual pair figures with fixed colorbar
# =============================================================================
print("\nCreating supplemental density scatter plots...")

PAIRS = [
    ('CD40', 'CD163'),
    ('IL-18', 'CD163'),
    ('IL-1b', 'CD163'),
    ('iNOS', 'CD163'),
    ('HLA-DR', 'IDO1'),
    ('CD40', 'IDO1'),
    ('IL-18', 'IDO1'),
    ('IL-1b', 'IDO1'),
    ('iNOS', 'IDO1'),
]

SUPP_PATH = OUTPUT_PATH / 'Supplemental_Pairs'
SUPP_PATH.mkdir(parents=True, exist_ok=True)

for m1, m2 in PAIRS:
    print(f"  {m1} vs {m2}...")

    # Axis limits based on Bx4
    x_max_pair = macs_bx4[m2].quantile(0.99)
    y_max_pair = macs_bx4[m1].quantile(0.99)
    x_max_pair = max(np.ceil(x_max_pair / 10) * 10, 20)
    y_max_pair = max(np.ceil(y_max_pair / 10) * 10, 20)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    scatter_obj = None
    for idx, (data, tp) in enumerate([(macs_bx2, 'Bx2'), (macs_bx4, 'Bx4')]):
        ax = axes[idx]

        x = data[m2].values
        y = data[m1].values

        scatter = density_scatter(ax, x, y, cmap='viridis', s=4, alpha=0.8)
        if scatter is not None:
            scatter_obj = scatter

        ax.set_xlim(0, x_max_pair)
        ax.set_ylim(0, y_max_pair)
        ax.set_xlabel(m2, fontsize=13, fontweight='bold')
        ax.set_ylabel(m1, fontsize=13, fontweight='bold')
        ax.set_title(f'{tp} (n = {len(data):,})', fontsize=14, fontweight='bold')
        ax.grid(False)

    # Colorbar to the right
    if scatter_obj is not None:
        fig.subplots_adjust(right=0.88)
        cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(scatter_obj, cax=cbar_ax)
        cbar.set_label('Density', fontsize=12, fontweight='bold')
        cbar.ax.tick_params(labelsize=10)
        
    plt.savefig(SUPP_PATH / f'Fig_Supp_{m1}_vs_{m2}_Density.pdf', dpi=300, bbox_inches='tight')
    plt.close()

print("\nAll individual pair figures saved.")

# =============================================================================
# COMBINED: CD163 pairs (2x2 with proper spacing)
# =============================================================================
print("\nCreating combined CD163 pairs figure...")

cd163_pairs = [('CD40', 'CD163'), ('IL-18', 'CD163'), ('IL-1b', 'CD163'), ('iNOS', 'CD163')]

fig = plt.figure(figsize=(14, 12))

# Create grid with space for colorbars
for pair_idx, (m1, m2) in enumerate(cd163_pairs):
    row = pair_idx // 2
    col = pair_idx % 2

    x_max_pair = macs_bx4[m2].quantile(0.99)
    y_max_pair = macs_bx4[m1].quantile(0.99)
    x_max_pair = max(np.ceil(x_max_pair / 10) * 10, 20)
    y_max_pair = max(np.ceil(y_max_pair / 10) * 10, 20)

    for tp_idx, (data, tp) in enumerate([(macs_bx2, 'Bx2'), (macs_bx4, 'Bx4')]):
        # Calculate subplot position
        subplot_idx = row * 4 + col * 2 + tp_idx + 1
        ax = fig.add_subplot(2, 4, subplot_idx)

        x = data[m2].values
        y = data[m1].values

        scatter = density_scatter(ax, x, y, cmap='viridis', s=3, alpha=0.7)

        ax.set_xlim(0, x_max_pair)
        ax.set_ylim(0, y_max_pair)
        ax.set_xlabel(m2, fontsize=10, fontweight='bold')
        ax.set_ylabel(m1, fontsize=10, fontweight='bold')
        ax.set_title(f'{tp} (n={len(data):,})', fontsize=11, fontweight='bold')
        ax.grid(False)

plt.suptitle('M1 Markers vs CD163 (M2)', fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(SUPP_PATH / 'Fig_Supp_CD163_AllM1_Combined.pdf', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: Fig_Supp_CD163_AllM1_Combined.pdf")

# =============================================================================
# COMBINED: IDO1 pairs (2x2 with proper spacing)
# =============================================================================
print("\nCreating combined IDO1 pairs figure...")

ido1_pairs = [('HLA-DR', 'IDO1'), ('CD40', 'IDO1'), ('IL-18', 'IDO1'), ('iNOS', 'IDO1')]

fig = plt.figure(figsize=(14, 12))

for pair_idx, (m1, m2) in enumerate(ido1_pairs):
    row = pair_idx // 2
    col = pair_idx % 2

    x_max_pair = macs_bx4[m2].quantile(0.99)
    y_max_pair = macs_bx4[m1].quantile(0.99)
    x_max_pair = max(np.ceil(x_max_pair / 10) * 10, 20)
    y_max_pair = max(np.ceil(y_max_pair / 10) * 10, 20)

    for tp_idx, (data, tp) in enumerate([(macs_bx2, 'Bx2'), (macs_bx4, 'Bx4')]):
        subplot_idx = row * 4 + col * 2 + tp_idx + 1
        ax = fig.add_subplot(2, 4, subplot_idx)

        x = data[m2].values
        y = data[m1].values

        scatter = density_scatter(ax, x, y, cmap='viridis', s=3, alpha=0.7)

        ax.set_xlim(0, x_max_pair)
        ax.set_ylim(0, y_max_pair)
        ax.set_xlabel(m2, fontsize=10, fontweight='bold')
        ax.set_ylabel(m1, fontsize=10, fontweight='bold')
        ax.set_title(f'{tp} (n={len(data):,})', fontsize=11, fontweight='bold')
        ax.grid(False)

plt.suptitle('M1 Markers vs IDO1 (M2)', fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(SUPP_PATH / 'Fig_Supp_IDO1_AllM1_Combined.pdf', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: Fig_Supp_IDO1_AllM1_Combined.pdf")

# =============================================================================
# IL-1b colored scatter (with fixed colorbar)
# =============================================================================
print("\nCreating IL-1b colored scatter with fixed colorbar...")

from matplotlib.colors import Normalize

il1b_min = macs['IL-1b'].quantile(0.01)
il1b_max = macs['IL-1b'].quantile(0.99)
norm = Normalize(vmin=il1b_min, vmax=il1b_max)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

scatter_obj = None
for idx, (data, tp) in enumerate([(macs_bx2, 'Bx2'), (macs_bx4, 'Bx4')]):
    ax = axes[idx]

    data_sorted = data.sort_values('IL-1b', ascending=True)

    scatter = ax.scatter(data_sorted['CD163'], data_sorted['HLA-DR'],
                        c=data_sorted['IL-1b'], cmap='inferno',
                        s=4, alpha=0.7, norm=norm,
                        edgecolors='none', rasterized=True)
    scatter_obj = scatter

    ax.set_xlim(0, x_max)
    ax.set_ylim(0, y_max)
    ax.set_xlabel('CD163', fontsize=13, fontweight='bold')
    ax.set_ylabel('HLA-DR', fontsize=13, fontweight='bold')
    ax.set_title(f'{tp} (n = {len(data):,})', fontsize=14, fontweight='bold')
    ax.grid(False)

# Colorbar to the right
fig.subplots_adjust(right=0.88)
cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
cbar = fig.colorbar(scatter_obj, cax=cbar_ax)
cbar.set_label('IL-1β', fontsize=12, fontweight='bold')
cbar.ax.tick_params(labelsize=10)

plt.savefig(OUTPUT_PATH / 'Fig_HLADR_CD163_ColorByIL1b.pdf', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: Fig_HLADR_CD163_ColorByIL1b.pdf")

# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("OUTPUT FILES")
print("=" * 70)

print(f"""
Main figure:
  - Fig_Main_CD163_HLADR_Density.pdf

IL-1b colored:
  - Fig_HLADR_CD163_ColorByIL1b.pdf

Supplemental (in Supplemental_Pairs/):
  Individual pairs:
    - Fig_Supp_CD40_vs_CD163_Density.pdf
    - Fig_Supp_IL-18_vs_CD163_Density.pdf
    - Fig_Supp_IL-1b_vs_CD163_Density.pdf
    - Fig_Supp_iNOS_vs_CD163_Density.pdf
    - Fig_Supp_HLA-DR_vs_IDO1_Density.pdf
    - Fig_Supp_CD40_vs_IDO1_Density.pdf
    - Fig_Supp_IL-18_vs_IDO1_Density.pdf
    - Fig_Supp_IL-1b_vs_IDO1_Density.pdf
    - Fig_Supp_iNOS_vs_IDO1_Density.pdf

  Combined:
    - Fig_Supp_CD163_AllM1_Combined.pdf
    - Fig_Supp_IDO1_AllM1_Combined.pdf

All colorbars are positioned to the right, outside the plot area.
""")
