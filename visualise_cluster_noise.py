#!/usr/bin/env python3

import argparse
import os
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys

def parse_args():
    parser = argparse.ArgumentParser(
        description="Visualize cluster noise by plotting expression profiles across tissues."
    )
    parser.add_argument(
        "--expression-file",
        required=True,
        help="Raw expression CSV (gene,tissue,median_tpm...)"
    )
    parser.add_argument(
        "--data-dir",
        default="results",
        help="Directory containing gene cluster files (default: results)"
    )
    parser.add_argument(
        "--clusters",
        nargs="+",
        required=True,
        help="Clusters to visualize. Can be names (e.g., 'k7 c2'), paths to gene list files, or 'all' to plot every cluster in the data directory."
    )
    parser.add_argument(
        "--stat",
        choices=["mean", "median"],
        default="median",
        help="Statistic to use from expression file (default: median)"
    )
    parser.add_argument(
        "--output-file",
        default="cluster_noise_visualization.png",
        help="Output image filename"
    )
    return parser.parse_args()

def load_expression_matrix(filepath, stat):
    """Loads and pivots the expression data to Gene x Tissue matrix."""
    if not os.path.exists(filepath):
        print(f"❌ Error: Expression file not found: {filepath}")
        return None
    
    try:
        df = pd.read_csv(filepath)
        value_col = f"{stat}_tpm"
        
        if value_col not in df.columns:
            # Fallback for already pivoted matrices
            if filepath.endswith('.csv'):
                # Check if first col is gene and others are tissues
                # If it's already pivoted, it might look like gtex_cluster_v2.py output
                print(f"⚠ Warning: '{value_col}' not found. Assuming file is an already pivoted matrix.")
                pivoted = df.set_index(df.columns[0])
                return pivoted.T
            return None
            
        # Pivot to (Genes x Tissues) then Transpose to (Tissues x Genes) for plotting logic
        matrix = df.pivot(index="gene", columns="tissue", values=value_col).fillna(0)
        return matrix.T
    except Exception as e:
        print(f"❌ Error loading expression data: {e}")
        return None

def get_genes_for_cluster(cluster_id, data_dir):
    """Resolves a cluster name or file path to a list of genes."""
    # Check if it's a direct file path first
    if os.path.exists(cluster_id):
        try:
            with open(cluster_id, 'r') as f:
                return [line.strip() for line in f if line.strip()], os.path.basename(cluster_id)
        except:
            pass

    # Try to parse 'k7 c2' or 'k7_c2' or 'k7 cluster 2'
    # Flexible match for k followed by digits and then cluster/c followed by digits
    m = re.search(r"k(\d+).*[cluster|c](\d+)", cluster_id, re.I)
    if m:
        k, c = m.group(1), m.group(2)
        # Match the filename pattern: cluster_k{k}_cluster_{c}_genes.txt
        target = os.path.join(data_dir, f"cluster_k{k}_cluster_{c}_genes.txt")
        if os.path.exists(target):
            with open(target, 'r') as f:
                return [line.strip() for line in f if line.strip()], f"k{k} Cluster {c}"
        else:
            # Try a slightly different variation if the first one fails
            target_v2 = os.path.join(data_dir, f"cluster_k{k}_c{c}_genes.txt")
            if os.path.exists(target_v2):
                with open(target_v2, 'r') as f:
                    return [line.strip() for line in f if line.strip()], f"k{k} Cluster {c}"
            else:
                print(f"⚠ Warning: Could not find cluster file for '{cluster_id}' in {data_dir}")
    
    return None, None

def plot_cluster_noise(df, cluster_genes, title, ax):
    """Plots the expression profile for a specific list of genes."""
    # Ensure genes are in columns
    valid_genes = [g for g in cluster_genes if g in df.columns]
    
    if not valid_genes:
        ax.text(0.5, 0.5, f'No gene overlap for {title}', ha='center', va='center')
        ax.set_title(title)
        return

    # Subset data
    subset = df[valid_genes]
    
    # Z-score standardization (relative to the gene's own distribution across tissues)
    # matches gtex_cluster_v2 standardization logic but applied per gene here
    subset_z = (subset - subset.mean()) / (subset.std() + 1e-8)
    
    # Melt for Seaborn
    # subset_z index is Tissue, columns are Genes
    plot_data = subset_z.reset_index().melt(id_vars=subset_z.index.name or 'index', 
                                            var_name='Gene', 
                                            value_name='Z-Score')
    tissue_col = plot_data.columns[0]

    # Plot individual gene lines
    sns.lineplot(data=plot_data, x=tissue_col, y='Z-Score', units='Gene', estimator=None, 
                 ax=ax, color='black', alpha=0.2, linewidth=0.8)
    
    # Highlight the mean trend
    sns.lineplot(data=plot_data, x=tissue_col, y='Z-Score', 
                 ax=ax, color='red', linewidth=2, label='Cluster Mean')

    ax.set_title(f"{title} (n={len(valid_genes)})")
    ax.set_ylabel("Z-Score")
    ax.set_xlabel("")
    ax.tick_params(axis='x', rotation=90, labelsize=7)
    ax.grid(True, alpha=0.15)

def main():
    args = parse_args()
    
    # 1. Load Expression Data
    df_t = load_expression_matrix(args.expression_file, args.stat)
    if df_t is None:
        sys.exit(1)

    # 2. Collect Cluster Data
    tasks = []
    requested_clusters = args.clusters
    
    # Handle "all" keyword
    if "all" in [c.lower() for c in requested_clusters]:
        import glob
        # Find all files matching the pattern
        pattern = os.path.join(args.data_dir, "cluster_k*_cluster_*_genes.txt")
        all_files = sorted(glob.glob(pattern))
        
        if not all_files:
            # Try alternate naming pattern
            pattern_v2 = os.path.join(args.data_dir, "cluster_k*_c*_genes.txt")
            all_files = sorted(glob.glob(pattern_v2))
            
        if all_files:
            # Custom sort to ensure k7 c10 comes after k7 c2
            def sort_key(path):
                fname = os.path.basename(path)
                m = re.search(r"k(\d+).*[cluster|c](\d+)", fname, re.I)
                if m:
                    return (int(m.group(1)), int(m.group(2)))
                return (0, 0)
            
            all_files.sort(key=sort_key)
            requested_clusters = all_files
        else:
            print(f"⚠ Warning: 'all' specified but no cluster files found in {args.data_dir}")

    for cid in requested_clusters:
        genes, name = get_genes_for_cluster(cid, args.data_dir)
        if genes:
            tasks.append((genes, name))
    
    if not tasks:
        print("❌ Error: No valid clusters or gene files found to visualize.")
        print(f"Tip: Ensure --data-dir ({args.data_dir}) is correct and contains files like 'cluster_k7_cluster_1_genes.txt'")
        sys.exit(1)

    # 3. Setup Multi-panel Plot
    n = len(tasks)
    fig, axes = plt.subplots(n, 1, figsize=(14, 4 * n + 2), sharex=True, squeeze=False)
    axes = axes.flatten()

    # 4. Generate Plots
    print(f"📊 Generating noise visualization for {n} cluster(s)...")
    for i, (genes, name) in enumerate(tasks):
        plot_cluster_noise(df_t, genes, name, axes[i])

    plt.tight_layout()
    plt.savefig(args.output_file, dpi=200, bbox_inches='tight')
    print(f"✅ Success! Plot saved to: {args.output_file}")

if __name__ == "__main__":
    main()