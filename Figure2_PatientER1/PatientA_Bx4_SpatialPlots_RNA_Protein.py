#!/usr/bin/env python3
"""
Spatial Expression Plots for Patient A Bx4
- RNA: Selected markers (ESR1, TGFBR2, JAG1, etc.)
- Protein: All markers
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from pathlib import Path
import gzip
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================

# Paths
RNA_EXPR_PATH = Path("supplementary_input_data/PatientA_Bx4_RNA_Expression.csv")
RNA_POLYGON_PATH = Path("supplementary_input_data/Polygons/RNA/R1134_322078-polygons.csv")

PROTEIN_EXPR_PATH = Path("supplementary_input_data/CosmX protein\\Bx4_0000459956\\0000459956_01112024_exprMat_file.csv.gz")
PROTEIN_POLYGON_PATH = Path("supplementary_input_data/Polygons/Protein/0000459956_01112024-polygons.csv.gz")

# Output directories
OUTPUT_BASE = Path("Spatial_Expression_Plots")
RNA_OUTPUT = OUTPUT_BASE / "RNA_Spatial"
PROTEIN_OUTPUT = OUTPUT_BASE / "Protein_Spatial"

RNA_OUTPUT.mkdir(parents=True, exist_ok=True)
PROTEIN_OUTPUT.mkdir(parents=True, exist_ok=True)

# RNA target genes
RNA_GENES = ['ESR1', 'TGFBR2', 'JAG1', 'SNAI1', 'MKI67', 'MX1', 'STAT1',
             'HLA-A', 'HLA-B', 'SPP1', 'WNT5A', 'HGF', 'MET', 'MIF',
             'CD74', 'CDH1']

# Publication settings
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 10
plt.rcParams['figure.facecolor'] = 'white'

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def create_spatial_plot_centroid(df, x_col, y_col, value_col, title, output_path,
                                 cmap='viridis', figsize=(12, 10), point_size=1):
    """Create spatial scatter plot colored by expression"""
    fig, ax = plt.subplots(figsize=figsize)

    # Get values and normalize
    values = df[value_col].values
    vmin, vmax = np.percentile(values[values > 0], [5, 95]) if (values > 0).any() else (0, 1)
    if vmin == vmax:
        vmax = vmin + 0.1

    # Sort by expression (low on bottom, high on top)
    df_sorted = df.sort_values(value_col)

    scatter = ax.scatter(df_sorted[x_col], df_sorted[y_col],
                        c=df_sorted[value_col], cmap=cmap,
                        s=point_size, alpha=0.8,
                        vmin=0, vmax=vmax,
                        edgecolors='none', rasterized=True)

    ax.set_aspect('equal')
    ax.set_xlabel('X (pixels)', fontsize=11)
    ax.set_ylabel('Y (pixels)', fontsize=11)
    ax.set_title(f'{title}\nPatient A Bx4', fontsize=13, fontweight='bold')

    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Expression', fontsize=10)

    # Clean up
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

# =============================================================================
# LOAD RNA DATA
# =============================================================================
print("=" * 70)
print("Loading RNA Data...")
print("=" * 70)

rna_expr = pd.read_csv(RNA_EXPR_PATH)
print(f"RNA expression data: {len(rna_expr):,} cells")

# =============================================================================
# CREATE RNA SPATIAL PLOTS
# =============================================================================
print("\n" + "=" * 70)
print("Creating RNA Spatial Plots...")
print("=" * 70)

for gene in RNA_GENES:
    if gene in rna_expr.columns:
        print(f"  Plotting {gene}...")
        output_path = RNA_OUTPUT / f"PatientA_Bx4_RNA_{gene}_Spatial.pdf"

        create_spatial_plot_centroid(
            rna_expr,
            'x_global_px', 'y_global_px',
            gene,
            gene,
            output_path,
            cmap='magma',
            figsize=(12, 10),
            point_size=2
        )
    else:
        print(f"  WARNING: {gene} not found in RNA data")

print(f"\nRNA plots saved to: {RNA_OUTPUT}")

# =============================================================================
# LOAD PROTEIN DATA
# =============================================================================
print("\n" + "=" * 70)
print("Loading Protein Data...")
print("=" * 70)

# Load expression
with gzip.open(PROTEIN_EXPR_PATH, 'rt') as f:
    protein_expr = pd.read_csv(f)
print(f"Protein expression data: {len(protein_expr):,} cells")

# Load polygons to get centroids
with gzip.open(PROTEIN_POLYGON_PATH, 'rt') as f:
    protein_poly = pd.read_csv(f)
print(f"Protein polygon rows: {len(protein_poly):,}")

# Calculate centroids
protein_centroids = protein_poly.groupby(['fov', 'cellID']).agg({
    'x_global_px': 'mean',
    'y_global_px': 'mean'
}).reset_index()

# Create match key for merging
protein_centroids['match_key'] = protein_centroids['fov'].astype(str) + '_' + protein_centroids['cellID'].astype(str)
protein_expr['match_key'] = protein_expr['fov'].astype(str) + '_' + protein_expr['cell_ID'].astype(str)

# Merge
protein_data = protein_expr.merge(
    protein_centroids[['match_key', 'x_global_px', 'y_global_px']],
    on='match_key',
    how='inner'
)
print(f"Merged protein data: {len(protein_data):,} cells")

# =============================================================================
# GET PROTEIN MARKERS
# =============================================================================
# Exclude metadata columns
exclude_cols = ['fov', 'cell_ID', 'match_key', 'x_global_px', 'y_global_px',
                'Channel-CD45', 'Channel-DNA', 'Channel-G', 'Channel-Membrane',
                'Channel-PanCK', 'Ms IgG1', 'Rb IgG']
protein_markers = [c for c in protein_expr.columns if c not in exclude_cols]
print(f"\nProtein markers to plot: {len(protein_markers)}")

# =============================================================================
# CREATE PROTEIN SPATIAL PLOTS
# =============================================================================
print("\n" + "=" * 70)
print("Creating Protein Spatial Plots...")
print("=" * 70)

for marker in protein_markers:
    if marker in protein_data.columns:
        print(f"  Plotting {marker}...")

        # Clean marker name for filename
        safe_name = marker.replace('/', '_').replace(' ', '_')
        output_path = PROTEIN_OUTPUT / f"PatientA_Bx4_Protein_{safe_name}_Spatial.pdf"

        create_spatial_plot_centroid(
            protein_data,
            'x_global_px', 'y_global_px',
            marker,
            marker,
            output_path,
            cmap='viridis',
            figsize=(12, 10),
            point_size=3
        )

print(f"\nProtein plots saved to: {PROTEIN_OUTPUT}")

# =============================================================================
# CREATE COMBINED SUMMARY FIGURES
# =============================================================================
print("\n" + "=" * 70)
print("Creating Summary Figures...")
print("=" * 70)

# RNA Summary (4x4 grid)
print("Creating RNA summary grid...")
fig, axes = plt.subplots(4, 4, figsize=(20, 20))
axes = axes.flatten()

for idx, gene in enumerate(RNA_GENES):
    ax = axes[idx]
    if gene in rna_expr.columns:
        values = rna_expr[gene].values
        vmax = np.percentile(values[values > 0], 95) if (values > 0).any() else 1

        scatter = ax.scatter(rna_expr['x_global_px'], rna_expr['y_global_px'],
                            c=rna_expr[gene], cmap='magma',
                            s=0.5, alpha=0.7, vmin=0, vmax=vmax,
                            edgecolors='none', rasterized=True)
        ax.set_aspect('equal')
        ax.set_title(gene, fontsize=12, fontweight='bold')
        ax.axis('off')

plt.suptitle('RNA Expression - Patient A Bx4', fontsize=16, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig(OUTPUT_BASE / 'PatientA_Bx4_RNA_Summary_Grid.pdf', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: PatientA_Bx4_RNA_Summary_Grid.pdf")

# Protein Summary - Key markers (4x4 grid)
print("Creating Protein summary grid...")
key_protein_markers = ['EpCAM', 'CD68', 'CD3', 'CD8', 'HLA-DR', 'CD163',
                       'Bcl-2', 'EGFR', 'Beta-catenin', 'Ki-67', 'PD-L1',
                       'IDO1', 'Vimentin', 'Fibronectin', 'SMA', 'ICAM1']

fig, axes = plt.subplots(4, 4, figsize=(20, 20))
axes = axes.flatten()

for idx, marker in enumerate(key_protein_markers):
    ax = axes[idx]
    if marker in protein_data.columns:
        values = protein_data[marker].values
        vmax = np.percentile(values, 95)

        scatter = ax.scatter(protein_data['x_global_px'], protein_data['y_global_px'],
                            c=protein_data[marker], cmap='viridis',
                            s=0.5, alpha=0.7, vmin=0, vmax=vmax,
                            edgecolors='none', rasterized=True)
        ax.set_aspect('equal')
        ax.set_title(marker, fontsize=12, fontweight='bold')
        ax.axis('off')

plt.suptitle('Protein Expression - Patient A Bx4', fontsize=16, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig(OUTPUT_BASE / 'PatientA_Bx4_Protein_Summary_Grid.pdf', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: PatientA_Bx4_Protein_Summary_Grid.pdf")

# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
Output directory: {OUTPUT_BASE}

RNA Spatial Plots ({len(RNA_GENES)} genes):
  Directory: {RNA_OUTPUT}
  Genes: {', '.join(RNA_GENES)}

Protein Spatial Plots ({len(protein_markers)} markers):
  Directory: {PROTEIN_OUTPUT}

Summary Grids:
  - PatientA_Bx4_RNA_Summary_Grid.pdf
  - PatientA_Bx4_Protein_Summary_Grid.pdf
""")
