#!/usr/bin/env python3
"""
Spatial Protein Plots using cell centroids with SHARED color scale
Patient A - Bx2 and Bx4 only
Markers: Ki-67, HLA-DR, CD163, EpCAM
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from pathlib import Path
import gzip
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_PATH = Path("supplementary_input_data/CosmX protein\\")

DATA_PATHS = {
    'Bx2': {
        'expr': BASE_PATH / "Bx2_0000459948/0000459948_01112024_exprMat_file.csv.gz",
        'poly': BASE_PATH / "Bx2_0000459948/0000459948_01112024-polygons.csv.gz"
    },
    'Bx4': {
        'expr': BASE_PATH / "Bx4_0000459956/0000459956_01112024_exprMat_file.csv.gz",
        'poly': BASE_PATH / "Bx4_0000459956/0000459956_01112024-polygons.csv.gz"
    }
}

# Target markers
TARGET_MARKERS = ['Ki-67', 'HLA-DR', 'CD163', 'EpCAM']

# Output directory
OUTPUT_PATH = Path("Spatial_Expression_Plots/Protein_Centroids_Bx2_Bx4")
OUTPUT_PATH.mkdir(parents=True, exist_ok=True)

# Publication settings
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 11
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

# =============================================================================
# LOAD ALL DATA
# =============================================================================

print("=" * 70)
print("Loading All Data...")
print("=" * 70)

all_data = {}

for timepoint, paths in DATA_PATHS.items():
    print(f"\n  Loading {timepoint}...")

    # Load expression
    with gzip.open(paths['expr'], 'rt') as f:
        expr_df = pd.read_csv(f)

    # Load polygons
    with gzip.open(paths['poly'], 'rt') as f:
        poly_df = pd.read_csv(f)

    # Calculate centroids
    centroids = poly_df.groupby(['fov', 'cellID']).agg({
        'x_global_px': 'mean',
        'y_global_px': 'mean'
    }).reset_index()

    # Create match key
    centroids['cell_key'] = centroids['fov'].astype(str) + '_' + centroids['cellID'].astype(str)
    expr_df['cell_key'] = expr_df['fov'].astype(str) + '_' + expr_df['cell_ID'].astype(str)

    # Merge expression with centroids
    merged = expr_df.merge(centroids[['cell_key', 'x_global_px', 'y_global_px']],
                          on='cell_key', how='inner')

    all_data[timepoint] = merged
    print(f"    {timepoint}: {len(merged):,} cells")

# =============================================================================
# CALCULATE SHARED COLOR SCALES
# =============================================================================

print("\n" + "=" * 70)
print("Calculating Shared Color Scales...")
print("=" * 70)

# Get global min/max for each marker across both timepoints
color_scales = {}

for marker in TARGET_MARKERS:
    all_values = []
    for timepoint in ['Bx2', 'Bx4']:
        if marker in all_data[timepoint].columns:
            all_values.extend(all_data[timepoint][marker].values)

    vmin = 0
    vmax = np.percentile(all_values, 98) if len(all_values) > 0 else 1
    if vmax == 0:
        vmax = 1

    color_scales[marker] = {'vmin': vmin, 'vmax': vmax}
    print(f"  {marker}: vmin={vmin:.1f}, vmax={vmax:.1f}")

# =============================================================================
# CREATE COMPARISON FIGURES (2 panels per marker)
# =============================================================================

print("\n" + "=" * 70)
print("Creating Comparison Figures...")
print("=" * 70)

for marker in TARGET_MARKERS:
    print(f"\n  Creating {marker} comparison...")

    fig, axes = plt.subplots(1, 2, figsize=(18, 8))

    vmin = color_scales[marker]['vmin']
    vmax = color_scales[marker]['vmax']
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.get_cmap('viridis')

    for idx, timepoint in enumerate(['Bx2', 'Bx4']):
        ax = axes[idx]
        data = all_data[timepoint]

        if marker not in data.columns:
            ax.text(0.5, 0.5, f'{marker} not available', ha='center', va='center',
                   transform=ax.transAxes, fontsize=14)
            ax.set_title(f'{timepoint}', fontsize=14, fontweight='bold')
            continue

        # Sort by expression (low values plotted first, high on top)
        data_sorted = data.sort_values(marker)

        scatter = ax.scatter(data_sorted['x_global_px'],
                            data_sorted['y_global_px'],
                            c=data_sorted[marker],
                            cmap=cmap,
                            norm=norm,
                            s=3,
                            alpha=0.8,
                            edgecolors='none',
                            rasterized=True)

        ax.set_aspect('equal')
        ax.set_xlabel('X (pixels)', fontsize=12)
        ax.set_ylabel('Y (pixels)', fontsize=12)

        # Calculate stats
        mean_expr = data[marker].mean()

        ax.set_title(f'{timepoint} (n = {len(data):,})\nMean = {mean_expr:.1f}',
                    fontsize=13, fontweight='bold')

    # Add shared colorbar
    fig.subplots_adjust(right=0.92)
    cbar_ax = fig.add_axes([0.94, 0.15, 0.015, 0.7])
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label(f'{marker} Expression', fontsize=12, fontweight='bold')

    plt.suptitle(f'{marker} Expression - Patient A\n(Bx2 → Bx4, Shared Color Scale)',
                fontsize=16, fontweight='bold', y=1.02)

    output_file = OUTPUT_PATH / f"PatientA_Bx2_Bx4_{marker.replace('-', '')}_Centroids.pdf"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {output_file.name}")

# =============================================================================
# CREATE INDIVIDUAL PLOTS
# =============================================================================

print("\n" + "=" * 70)
print("Creating Individual Plots...")
print("=" * 70)

for timepoint in ['Bx2', 'Bx4']:
    data = all_data[timepoint]

    for marker in TARGET_MARKERS:
        if marker not in data.columns:
            continue

        print(f"  Plotting {timepoint} - {marker}...")

        fig, ax = plt.subplots(figsize=(12, 10))

        vmin = color_scales[marker]['vmin']
        vmax = color_scales[marker]['vmax']
        norm = Normalize(vmin=vmin, vmax=vmax)
        cmap = cm.get_cmap('viridis')

        # Sort by expression
        data_sorted = data.sort_values(marker)

        scatter = ax.scatter(data_sorted['x_global_px'],
                            data_sorted['y_global_px'],
                            c=data_sorted[marker],
                            cmap=cmap,
                            norm=norm,
                            s=4,
                            alpha=0.8,
                            edgecolors='none',
                            rasterized=True)

        ax.set_aspect('equal')
        ax.set_xlabel('X (pixels)', fontsize=12)
        ax.set_ylabel('Y (pixels)', fontsize=12)

        mean_expr = data[marker].mean()
        ax.set_title(f'{marker} - Patient A {timepoint}\n(n = {len(data):,}, Mean = {mean_expr:.1f})',
                    fontsize=14, fontweight='bold')

        # Colorbar
        cbar = plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(f'{marker} Expression', fontsize=11)

        plt.tight_layout()

        output_file = OUTPUT_PATH / f"PatientA_{timepoint}_{marker.replace('-', '')}_Centroids.pdf"
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()

# =============================================================================
# CREATE 4-MARKER GRID FOR EACH TIMEPOINT
# =============================================================================

print("\n" + "=" * 70)
print("Creating 4-Marker Grids per Timepoint...")
print("=" * 70)

for timepoint in ['Bx2', 'Bx4']:
    print(f"  Creating grid for {timepoint}...")

    fig, axes = plt.subplots(2, 2, figsize=(20, 18))
    axes = axes.flatten()

    data = all_data[timepoint]

    for idx, marker in enumerate(TARGET_MARKERS):
        ax = axes[idx]

        if marker not in data.columns:
            ax.text(0.5, 0.5, f'{marker} not available', ha='center', va='center',
                   transform=ax.transAxes, fontsize=14)
            continue

        vmin = color_scales[marker]['vmin']
        vmax = color_scales[marker]['vmax']
        norm = Normalize(vmin=vmin, vmax=vmax)
        cmap = cm.get_cmap('viridis')

        data_sorted = data.sort_values(marker)

        scatter = ax.scatter(data_sorted['x_global_px'],
                            data_sorted['y_global_px'],
                            c=data_sorted[marker],
                            cmap=cmap,
                            norm=norm,
                            s=2,
                            alpha=0.7,
                            edgecolors='none',
                            rasterized=True)

        ax.set_aspect('equal')
        ax.set_xlabel('X (pixels)', fontsize=11)
        ax.set_ylabel('Y (pixels)', fontsize=11)

        mean_expr = data[marker].mean()
        ax.set_title(f'{marker} (Mean = {mean_expr:.1f})', fontsize=13, fontweight='bold')

        cbar = plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Expression', fontsize=10)

    plt.suptitle(f'Patient A {timepoint} - All Markers (n = {len(data):,} cells)\nShared Color Scales',
                fontsize=16, fontweight='bold', y=1.01)
    plt.tight_layout()
    
    output_file = OUTPUT_PATH / f"PatientA_{timepoint}_AllMarkers_Grid.pdf"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {output_file.name}")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"""
Output Directory: {OUTPUT_PATH}

Data Summary:
  Bx2: {len(all_data['Bx2']):,} cells
  Bx4: {len(all_data['Bx4']):,} cells

Shared Color Scales:
  Ki-67:  0 - {color_scales['Ki-67']['vmax']:.1f}
  HLA-DR: 0 - {color_scales['HLA-DR']['vmax']:.1f}
  CD163:  0 - {color_scales['CD163']['vmax']:.1f}
  EpCAM:  0 - {color_scales['EpCAM']['vmax']:.1f}

Files Created:

Comparison Plots (2 timepoints side-by-side):
  - PatientA_Bx2_Bx4_Ki67_Centroids.pdf
  - PatientA_Bx2_Bx4_HLADR_Centroids.pdf
  - PatientA_Bx2_Bx4_CD163_Centroids.pdf
  - PatientA_Bx2_Bx4_EpCAM_Centroids.pdf

Individual Plots (8 files):
  Bx2: Ki67, HLADR, CD163, EpCAM
  Bx4: Ki67, HLADR, CD163, EpCAM

4-Marker Grids:
  - PatientA_Bx2_AllMarkers_Grid.pdf
  - PatientA_Bx4_AllMarkers_Grid.pdf
""")
