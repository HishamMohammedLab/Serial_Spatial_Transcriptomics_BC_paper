#!/usr/bin/env python3
"""
Patient A Spatial Expression by Major Lineage Groups
Creates composite scores for each lineage and visualizes spatially using cell polygons
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.collections import PolyCollection
import gzip
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Paths
DATA_PATH = Path("supplementary_input_data/LR_Scores")
OUTPUT_PATH = Path("Spatial_Lineage")
OUTPUT_PATH.mkdir(parents=True, exist_ok=True)

# Publication styling
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['savefig.facecolor'] = 'white'

# Pixel to micron conversion
PIXEL_SIZE = 0.18  # μm per pixel

# =============================================================================
# Define Major Lineage Groups
# =============================================================================
LINEAGE_GROUPS = {
    'T Cells': {
        'markers': ['CD3', 'CD4', 'CD8', 'CD45RA', 'FOXP3'],
        'cmap': 'Greens',
        'description': 'CD3+ T lymphocytes'
    },
    'B Cells': {
        'markers': ['CD19', 'CD20', 'CD27', 'IgD'],
        'cmap': 'Blues',
        'description': 'CD19/CD20+ B lymphocytes'
    },
    'Myeloid': {
        'markers': ['CD68', 'CD163', 'CD14', 'CD11b', 'CD11c', 'HLA-DR'],
        'cmap': 'Oranges',
        'description': 'Macrophages & Myeloid cells'
    },
    'NK Cells': {
        'markers': ['CD56', 'CD16'],
        'cmap': 'Purples',
        'description': 'Natural Killer cells'
    },
    'Tumor/Epithelial': {
        'markers': ['EpCAM', 'EGFR', 'Her2', 'pan-RAS', 'Ki-67'],
        'cmap': 'Reds',
        'description': 'Epithelial/Tumor markers'
    },
    'Stromal': {
        'markers': ['Vimentin', 'SMA', 'Fibronectin', 'CD31'],
        'cmap': 'YlOrBr',
        'description': 'Stromal & Endothelial'
    },
    'Immune Checkpoint': {
        'markers': ['PD-1', 'PD-L1', 'CTLA4', 'LAG3', 'Tim-3', 'VISTA'],
        'cmap': 'RdPu',
        'description': 'Checkpoint molecules'
    }
}

# =============================================================================
# Load Data
# =============================================================================
print("=" * 60)
print("Loading Data")
print("=" * 60)

def load_biopsy_data(biopsy_folder, biopsy_name):
    """Load all data for a biopsy."""
    print(f"\nLoading {biopsy_name}...")

    files = list(biopsy_folder.glob("*"))
    expr_file = [f for f in files if 'exprMat' in f.name][0]
    poly_file = [f for f in files if 'polygons' in f.name][0]

    # Load expression
    print(f"  Loading expression...")
    with gzip.open(expr_file, 'rt') as f:
        expr_df = pd.read_csv(f)
    print(f"    Rows: {len(expr_df):,}, FOVs: {expr_df['fov'].nunique()}")

    # Load polygons
    print(f"  Loading polygons...")
    with gzip.open(poly_file, 'rt') as f:
        poly_df = pd.read_csv(f)
    print(f"    Vertices: {len(poly_df):,}")

    # Convert polygon coordinates to microns
    poly_df['x_um'] = poly_df['x_global_px'] * PIXEL_SIZE
    poly_df['y_um'] = poly_df['y_global_px'] * PIXEL_SIZE

    # Create unique cell key: (fov, cellID)
    poly_df['cell_key'] = poly_df['fov'].astype(str) + '_' + poly_df['cellID'].astype(str)
    expr_df['cell_key'] = expr_df['fov'].astype(str) + '_' + expr_df['cell_ID'].astype(str)

    # Group polygons by cell_key
    print(f"  Grouping polygons by cell...")
    polygons = {}
    for cell_key, group in poly_df.groupby('cell_key'):
        vertices = group[['x_um', 'y_um']].values
        if len(vertices) >= 3:
            polygons[cell_key] = vertices

    # Create expression dict keyed by cell_key
    expr_dict = expr_df.set_index('cell_key').to_dict('index')

    # Check overlap
    poly_keys = set(polygons.keys())
    expr_keys = set(expr_dict.keys())
    overlap = poly_keys & expr_keys

    print(f"    Matched cells: {len(overlap):,}")

    return polygons, expr_dict, overlap, expr_df

# Load both biopsies
bx2_polys, bx2_expr, bx2_cells, bx2_df = load_biopsy_data(DATA_PATH / "Bx2_0000459948", "Bx2")
bx4_polys, bx4_expr, bx4_cells, bx4_df = load_biopsy_data(DATA_PATH / "Bx4_0000459956", "Bx4")

# =============================================================================
# Calculate Composite Scores
# =============================================================================
print("\n" + "=" * 60)
print("Calculating Lineage Composite Scores")
print("=" * 60)

def calculate_composite_scores_combined(bx2_df, bx4_df, bx2_cells, bx4_cells, markers):
    """Calculate composite scores normalized across both biopsies together."""
    # Get available markers
    available_bx2 = [m for m in markers if m in bx2_df.columns]
    available_bx4 = [m for m in markers if m in bx4_df.columns]
    available = list(set(available_bx2) & set(available_bx4))

    if not available:
        return None, None, []

    # Combine data from both biopsies for joint normalization
    bx2_marker_data = bx2_df.set_index('cell_key')[available].copy()
    bx4_marker_data = bx4_df.set_index('cell_key')[available].copy()

    combined = pd.concat([bx2_marker_data, bx4_marker_data])

    # Calculate z-scores using combined statistics
    combined_mean = combined.mean()
    combined_std = combined.std()

    bx2_z = (bx2_marker_data - combined_mean) / combined_std
    bx4_z = (bx4_marker_data - combined_mean) / combined_std

    # Calculate composite (mean of z-scores)
    bx2_composite = bx2_z.mean(axis=1).to_dict()
    bx4_composite = bx4_z.mean(axis=1).to_dict()

    return bx2_composite, bx4_composite, available

# Calculate composites for each lineage
composites = {}
for lineage, info in LINEAGE_GROUPS.items():
    print(f"\n{lineage}:")
    print(f"  Target markers: {info['markers']}")

    bx2_composite, bx4_composite, available = calculate_composite_scores_combined(
        bx2_df, bx4_df, bx2_cells, bx4_cells, info['markers']
    )

    if bx2_composite and bx4_composite:
        print(f"  Available markers: {available}")

        composites[lineage] = {
            'bx2': bx2_composite,
            'bx4': bx4_composite,
            'markers_used': available
        }

        # Summary stats
        bx2_vals = [bx2_composite[c] for c in bx2_cells if c in bx2_composite]
        bx4_vals = [bx4_composite[c] for c in bx4_cells if c in bx4_composite]
        print(f"  Bx2 mean: {np.mean(bx2_vals):.3f}, Bx4 mean: {np.mean(bx4_vals):.3f}")

# =============================================================================
# Create Spatial Plots
# =============================================================================
print("\n" + "=" * 60)
print("Creating Spatial Lineage Plots")
print("=" * 60)

def create_polygon_plot(ax, polygons, composite_dict, matched_cells, cmap, vmin, vmax):
    """Create plot with cell polygons colored by composite score."""

    poly_list = []
    colors = []

    for cell_key in matched_cells:
        if cell_key in polygons and cell_key in composite_dict:
            vertices = polygons[cell_key]
            score = composite_dict[cell_key]
            if not np.isnan(score):
                poly_list.append(vertices)
                colors.append(score)

    if len(poly_list) == 0:
        return None, 0

    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap_obj = plt.cm.get_cmap(cmap)

    collection = PolyCollection(
        poly_list,
        array=np.array(colors),
        cmap=cmap_obj,
        norm=norm,
        edgecolors='none',
        linewidths=0
    )

    ax.add_collection(collection)

    # Set limits
    all_x = np.concatenate([p[:, 0] for p in poly_list])
    all_y = np.concatenate([p[:, 1] for p in poly_list])

    ax.set_xlim(all_x.min() - 50, all_x.max() + 50)
    ax.set_ylim(all_y.min() - 50, all_y.max() + 50)

    return collection, len(poly_list)

# Create main figure: all lineages in a grid
n_lineages = len(LINEAGE_GROUPS)
fig, axes = plt.subplots(n_lineages, 2, figsize=(10, 3.5 * n_lineages))

for row, (lineage, info) in enumerate(LINEAGE_GROUPS.items()):
    print(f"\nPlotting {lineage}...")

    if lineage not in composites:
        print(f"  Skipped - no data")
        continue

    comp_data = composites[lineage]

    # Get scale from both biopsies
    bx2_vals = [comp_data['bx2'][c] for c in bx2_cells if c in comp_data['bx2'] and not np.isnan(comp_data['bx2'][c])]
    bx4_vals = [comp_data['bx4'][c] for c in bx4_cells if c in comp_data['bx4'] and not np.isnan(comp_data['bx4'][c])]
    all_vals = bx2_vals + bx4_vals

    vmin = np.percentile(all_vals, 2)
    vmax = np.percentile(all_vals, 98)
    print(f"  Scale: {vmin:.2f} - {vmax:.2f}")

    for col, (polys, comp, cells, name) in enumerate([
        (bx2_polys, comp_data['bx2'], bx2_cells, 'Bx2'),
        (bx4_polys, comp_data['bx4'], bx4_cells, 'Bx4')
    ]):
        ax = axes[row, col]

        result = create_polygon_plot(ax, polys, comp, cells, info['cmap'], vmin, vmax)

        if result[0]:
            n_cells = result[1]
            print(f"    {name}: {n_cells:,} cells")

        ax.set_aspect('equal')
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

        if row == 0:
            ax.set_title(name, fontsize=14, fontweight='bold', pad=10)
        if col == 0:
            ax.set_ylabel(f"{lineage}\n({len(comp_data['markers_used'])} markers)",
                         fontsize=10, fontweight='bold', labelpad=10)

        # Scale bar on first row, first column
        if row == 0 and col == 0:
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            scale_len = 500
            x0 = xlim[0] + (xlim[1]-xlim[0])*0.05
            y0 = ylim[0] + (ylim[1]-ylim[0])*0.05
            ax.plot([x0, x0+scale_len], [y0, y0], 'k-', lw=2)
            ax.text(x0+scale_len/2, y0-(ylim[1]-ylim[0])*0.03, '500 μm', ha='center', fontsize=8)

    # Colorbar
    cbar_ax = fig.add_axes([0.92, 0.88 - row*(0.85/n_lineages) - 0.08, 0.015, 0.85/n_lineages - 0.02])
    sm = plt.cm.ScalarMappable(cmap=info['cmap'],
                               norm=mcolors.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Z-score', fontsize=8)
    cbar.ax.tick_params(labelsize=7)

plt.suptitle('Patient A: Spatial Expression by Cell Lineage\n(Bx2 → Bx4)',
             fontsize=14, fontweight='bold', y=0.995)
plt.tight_layout(rect=[0, 0, 0.90, 0.98])
plt.savefig(OUTPUT_PATH / 'PatientA_Spatial_Lineage_All.pdf', dpi=150, bbox_inches='tight')
plt.close()
print("\nSaved: PatientA_Spatial_Lineage_All.pdf")

# =============================================================================
# Create Individual Lineage Plots (higher quality)
# =============================================================================
print("\n" + "=" * 60)
print("Creating Individual Lineage Plots")
print("=" * 60)

for lineage, info in LINEAGE_GROUPS.items():
    if lineage not in composites:
        continue

    print(f"\nCreating {lineage} plot...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    comp_data = composites[lineage]

    # Get scale
    bx2_vals = [comp_data['bx2'][c] for c in bx2_cells if c in comp_data['bx2'] and not np.isnan(comp_data['bx2'][c])]
    bx4_vals = [comp_data['bx4'][c] for c in bx4_cells if c in comp_data['bx4'] and not np.isnan(comp_data['bx4'][c])]
    all_vals = bx2_vals + bx4_vals

    vmin = np.percentile(all_vals, 2)
    vmax = np.percentile(all_vals, 98)

    for col, (polys, comp, cells, name) in enumerate([
        (bx2_polys, comp_data['bx2'], bx2_cells, 'Bx2'),
        (bx4_polys, comp_data['bx4'], bx4_cells, 'Bx4')
    ]):
        ax = axes[col]

        result = create_polygon_plot(ax, polys, comp, cells, info['cmap'], vmin, vmax)
        n_cells = result[1] if result[0] else 0

        ax.set_aspect('equal')
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_title(f'{name}\n({n_cells:,} cells)', fontsize=12, fontweight='bold')

        # Scale bar
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        x0 = xlim[0] + (xlim[1]-xlim[0])*0.05
        y0 = ylim[0] + (ylim[1]-ylim[0])*0.05
        ax.plot([x0, x0+500], [y0, y0], 'k-', lw=2)
        ax.text(x0+250, y0-(ylim[1]-ylim[0])*0.03, '500 μm', ha='center', fontsize=9)

    # Colorbar
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    sm = plt.cm.ScalarMappable(cmap=info['cmap'],
                               norm=mcolors.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Composite Z-score', fontsize=10)

    # Calculate mean change
    bx2_mean = np.mean(bx2_vals)
    bx4_mean = np.mean(bx4_vals)
    change = bx4_mean - bx2_mean
    arrow = "↑" if change > 0 else "↓"

    markers_str = ', '.join(comp_data['markers_used'])
    plt.suptitle(f"Patient A: {lineage} ({arrow} {abs(change):.2f} in Bx4)\n{info['description']}\nMarkers: {markers_str}",
                 fontsize=12, fontweight='bold', y=1.02)

    plt.tight_layout(rect=[0, 0, 0.88, 0.94])

    # Clean filename
    filename = lineage.replace(' ', '_').replace('/', '_')
    plt.savefig(OUTPUT_PATH / f'PatientA_Spatial_{filename}.pdf', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: PatientA_Spatial_{filename}.pdf")

# =============================================================================
# Summary Statistics
# =============================================================================
print("\n" + "=" * 60)
print("Lineage Expression Summary")
print("=" * 60)

summary_data = []
for lineage, info in LINEAGE_GROUPS.items():
    if lineage not in composites:
        continue

    comp_data = composites[lineage]

    bx2_vals = [comp_data['bx2'][c] for c in bx2_cells if c in comp_data['bx2'] and not np.isnan(comp_data['bx2'][c])]
    bx4_vals = [comp_data['bx4'][c] for c in bx4_cells if c in comp_data['bx4'] and not np.isnan(comp_data['bx4'][c])]

    bx2_mean = np.mean(bx2_vals)
    bx4_mean = np.mean(bx4_vals)
    change = bx4_mean - bx2_mean

    summary_data.append({
        'Lineage': lineage,
        'Bx2_Mean': bx2_mean,
        'Bx4_Mean': bx4_mean,
        'Change': change,
        'Direction': '↑' if change > 0 else '↓',
        'Markers': len(comp_data['markers_used'])
    })

    print(f"{lineage}: Bx2={bx2_mean:.3f}, Bx4={bx4_mean:.3f}, Δ={change:+.3f}")

# Save summary
summary_df = pd.DataFrame(summary_data)
summary_df.to_csv(OUTPUT_PATH / 'PatientA_Lineage_Summary.csv', index=False)
print(f"\nSaved: PatientA_Lineage_Summary.csv")

# =============================================================================
# Create Summary Bar Plot
# =============================================================================
print("\n" + "=" * 60)
print("Creating Summary Bar Plot")
print("=" * 60)

fig, ax = plt.subplots(figsize=(8, 5))

lineages = summary_df['Lineage'].tolist()
changes = summary_df['Change'].tolist()
colors_bar = ['#B2182B' if c > 0 else '#2166AC' for c in changes]

y_pos = np.arange(len(lineages))
bars = ax.barh(y_pos, changes, color=colors_bar, height=0.6, edgecolor='none')

ax.axvline(x=0, color='black', linestyle='-', linewidth=1)

ax.set_yticks(y_pos)
ax.set_yticklabels(lineages, fontsize=10)
ax.set_xlabel('Change in Composite Z-score (Bx2 → Bx4)', fontsize=10)
ax.set_title('Patient A: Lineage Expression Changes', fontsize=12, fontweight='bold')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add value labels
for i, (bar, val) in enumerate(zip(bars, changes)):
    x_pos = val + 0.02 if val >= 0 else val - 0.02
    ha = 'left' if val >= 0 else 'right'
    ax.text(x_pos, i, f'{val:+.3f}', va='center', ha=ha, fontsize=9)

ax.invert_yaxis()
plt.tight_layout()
plt.savefig(OUTPUT_PATH / 'PatientA_Lineage_Changes_Summary.pdf', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: PatientA_Lineage_Changes_Summary.pdf")

print("\n" + "=" * 60)
print("COMPLETE!")
print("=" * 60)
print(f"\nAll figures saved to: {OUTPUT_PATH}")
