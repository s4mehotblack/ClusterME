#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
from sklearn.metrics import silhouette_score, silhouette_samples

# ─────────────────────────────────────────────
# CLI ARGUMENTS
# ─────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="GTEx gene expression clustering pipeline (Improved)"
)

parser.add_argument(
    "--expression-file",
    required=True,
    help="All-tissues expression CSV (gene,tissue,mean_tpm,median_tpm,n_samples)"
)

parser.add_argument(
    "--data-dir",
    default="results",
    help="Output directory (default: results)"
)

parser.add_argument(
    "--stat",
    choices=["mean", "median"],
    default="median",
    help="Statistic to use (default: median)"
)

parser.add_argument(
    "--no-log2",
    action="store_true",
    help="Disable log2(x+1) transform (default: enabled)"
)

parser.add_argument(
    "--z-score",
    action="store_true",
    help="Apply row-wise Z-score standardization (Recommended for tissue-specificity patterns)"
)

parser.add_argument(
    "--k",
    type=int,
    nargs="+",
    required=True,
    help="One or more k values to cut the dendrogram (e.g. --k 6 7 8)"
)

parser.add_argument(
    "--label-genes",
    action="store_true",
    help="Show gene labels on heatmaps"
)

parser.add_argument(
    "--max-label-genes",
    type=int,
    default=300,
    help="Only label genes if count <= this (default: 300)"
)

parser.add_argument(
    "--min-cluster-size",
    type=int,
    default=10,
    help="Minimum genes for zoomed cluster heatmap (default: 10)"
)

parser.add_argument(
    "--cluster-sorted-heatmap",
    action="store_true",
    help="Sort genes by cluster ID for clearer visual separation (recommended for k exploration)"
)

args = parser.parse_args()

os.makedirs(args.data_dir, exist_ok=True)

log2_enabled = not args.no_log2
norm_tag = "zscore" if args.z_score else "raw"
tag = f"{args.stat}_{'log2' if log2_enabled else 'nolog'}_{norm_tag}"


# ─────────────────────────────────────────────
# LOAD & PREPARE MATRIX
# ─────────────────────────────────────────────
if not os.path.exists(args.expression_file):
    print(f"❌ Error: Expression file not found: {args.expression_file}")
    sys.exit(1)

print("Loading expression data...")

try:
    expr = pd.read_csv(args.expression_file)
except Exception as e:
    print(f"❌ Error reading expression file: {e}")
    sys.exit(1)

value_col = f"{args.stat}_tpm"

# Pivot to (Genes x Tissues)
matrix = (
    expr
    .pivot(index="gene", columns="tissue", values=value_col)
    .fillna(0)
)

# 1. Log Transform
if log2_enabled:
    matrix = np.log2(matrix + 1)

# 2. Z-Score Standardization (Row-wise)
if args.z_score:
    print("Applying Z-score standardization (row-wise)...")
    # subtract mean, divide by std (add small epsilon to avoid div/0)
    matrix = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1) + 1e-8, axis=0)
    
    # Set visuals for Z-score (diverging)
    cmap = "vlag"
    center = 0
    # Robust scaling for visualization (clip outliers at 2nd/98th percentiles)
    vmin = np.percentile(matrix.values, 2)
    vmax = np.percentile(matrix.values, 98)
else:
    # Visuals for raw magnitudes (sequential)
    cmap = "viridis"
    center = None
    vmin = np.percentile(matrix.values, 5)
    vmax = np.percentile(matrix.values, 95)

matrix.to_csv(f"{args.data_dir}/expression_matrix_{tag}.csv")
print(f"Matrix shape: {matrix.shape[0]} genes × {matrix.shape[1]} tissues")


# ─────────────────────────────────────────────
# HIERARCHICAL CLUSTERING
# ─────────────────────────────────────────────
print("Computing hierarchical clustering (Ward)...")

# Ward method requires Euclidean distance. 
# Even if using Z-scores (correlation-like), Ward on Z-scores is valid.
dist_matrix = pdist(matrix.values)
row_linkage = linkage(dist_matrix, method="ward")


# ─────────────────────────────────────────────
# PROCESS EACH K
# ─────────────────────────────────────────────
metrics = []

for k in args.k:
    print(f"\nProcessing k = {k}")

    # 1. Cut Dendrogram to get labels
    cluster_ids = fcluster(
        row_linkage,
        k,
        criterion="maxclust"
    )

    # 2. Compute Silhouette Score (Scientific Validity Fix)
    # We calculate score based on the actual Hierarchical labels, not a separate KMeans model
    if k > 1:
        sil = silhouette_score(matrix, cluster_ids)
        sample_silhouette_values = silhouette_samples(matrix, cluster_ids)
        print(f"  Silhouette Score: {sil:.4f}")
    else:
        sil = 0
        sample_silhouette_values = np.zeros(matrix.shape[0])
        
    metrics.append({"k": k, "silhouette": sil})

    # Prepare Cluster Dataframe
    cluster_series = pd.Series(
        cluster_ids,
        index=matrix.index,
        name="cluster"
    )
    
    # Create DataFrame with gene, cluster, and silhouette score
    cluster_df = cluster_series.reset_index()
    cluster_df["silhouette_score"] = sample_silhouette_values
    
    # Save standard cluster file
    cluster_df[["gene", "cluster"]].to_csv(f"{args.data_dir}/gene_clusters_k{k}_{tag}.csv", index=False)
    
    # Save detailed silhouette scores
    cluster_df.to_csv(f"{args.data_dir}/gene_silhouette_scores_k{k}_{tag}.csv", index=False)

    # Save Gene Lists
    for cid, sub in cluster_df.groupby("cluster"):
        sub["gene"].to_csv(
            f"{args.data_dir}/cluster_k{k}_cluster_{cid}_genes.txt",
            index=False,
            header=False
        )

    # ─────────────────────────────────────────────
    # GLOBAL HEATMAP
    # ─────────────────────────────────────────────
    print("  Generating global heatmap...")
    
    try:
        # Generate colors for side bar
        palette = sns.color_palette("tab10", k) if k <= 10 else sns.color_palette("husl", k)
        row_colors = cluster_series.map(dict(zip(range(1, k + 1), palette)))

        show_labels = (args.label_genes and matrix.shape[0] <= args.max_label_genes)

        if args.cluster_sorted_heatmap:
            # Cluster-sorted mode: sort genes by cluster ID for clearer visual separation
            print("  Using cluster-sorted visualization (genes ordered by cluster)")
            ordered_genes = cluster_df.sort_values("cluster")["gene"].tolist()
            sorted_matrix = matrix.loc[ordered_genes]
            
            # Create figure with appropriate size
            fig, ax = plt.subplots(figsize=(16, max(12, matrix.shape[0] * 0.03)))
            
            # Plot heatmap with sorted genes
            g = sns.heatmap(
                sorted_matrix,
                cmap=cmap,
                center=center,
                vmin=vmin,
                vmax=vmax,
                yticklabels=show_labels,
                ax=ax
            )
            
            # Add row colors sidebar manually
            row_colors_df = pd.DataFrame({
                'cluster': cluster_series.loc[ordered_genes].values,
                'color': [row_colors[gene] for gene in ordered_genes]
            })
            
            # Create a color bar on the right
            colorbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.8])
            for i, (cluster_id, color) in enumerate(zip(row_colors_df['cluster'].unique(), palette)):
                colorbar_ax.add_patch(plt.Rectangle((0, i), 1, 1, color=color))
            colorbar_ax.set_yticks(np.arange(len(palette)) + 0.5)
            colorbar_ax.set_yticklabels(row_colors_df['cluster'].unique())
            colorbar_ax.set_title('Cluster')
            colorbar_ax.axis('off')
            
            g.set_xticklabels(g.get_xticklabels(), rotation=45, ha="right", fontsize=8)
            if show_labels:
                g.set_yticklabels(g.get_yticklabels(), fontsize=6)
            
            plt.suptitle(f"Global Expression Clustering (Cluster-Sorted)\nk={k} | {tag}", y=1.02)
            plt.savefig(f"{args.data_dir}/global_heatmap_k{k}_{tag}.png", dpi=200, bbox_inches='tight')
            plt.close()
        else:
            # Hierarchical mode: use dendrogram ordering (scientifically rigorous)
            print("  Using hierarchical visualization (dendrogram ordering)")
            g = sns.clustermap(
                matrix,
                row_linkage=row_linkage,
                col_cluster=True,        # Cluster tissues to see relationships
                row_colors=row_colors,
                cmap=cmap,
                center=center,
                vmin=vmin,
                vmax=vmax,
                yticklabels=show_labels,
                dendrogram_ratio=(0.15, 0.1),
                figsize=(16, max(12, matrix.shape[0] * 0.03))
            )

            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=8)
            if show_labels:
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=6)

            plt.suptitle(f"Global Expression Clustering (Hierarchical)\nk={k} | {tag}", y=1.02)
            plt.savefig(f"{args.data_dir}/global_heatmap_k{k}_{tag}.png", dpi=200, bbox_inches='tight')
            plt.close()
    except Exception as e:
        print(f"  ❌ Failed to generate global heatmap for k={k}: {e}")

    # ─────────────────────────────────────────────
    # ZOOMED CLUSTER HEATMAPS
    # ─────────────────────────────────────────────
    print("  Generating zoomed cluster heatmaps...")

    for cid in sorted(cluster_series.unique()):
        genes = cluster_series[cluster_series == cid].index.tolist()
        n_cgenes = len(genes)
        if n_cgenes < args.min_cluster_size:
            print(f"  Skipping cluster {cid} (size {n_cgenes} < {args.min_cluster_size})")
            continue

        try:
            sub_matrix = matrix.loc[genes]
            
            if sub_matrix.nunique().sum() <= sub_matrix.shape[1]:
                print(f"  ⚠ Skipping cluster {cid}: too little variance for hierarchical clustering.")
                continue

            # Zoomed height: more generous allocation (0.25 inches per gene)
            z_h = max(6, min(100, n_cgenes * 0.25))
            z_w = max(12, min(50, matrix.shape[1] * 0.3))
            z_y_font = max(4, min(10, 600 / n_cgenes))
            z_x_font = max(6, min(12, 600 / matrix.shape[1]))

            h = sns.clustermap(
                sub_matrix,
                row_cluster=True,
                col_cluster=True,
                cmap=cmap,
                center=center,
                vmin=vmin,
                vmax=vmax,
                yticklabels=True, # Always show genes in zoomed view
                figsize=(z_w, z_h),
                dendrogram_ratio=(0.15, 0.1)
            )
            
            h.ax_heatmap.set_xticklabels(h.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=z_x_font)
            h.ax_heatmap.set_yticklabels(h.ax_heatmap.get_yticklabels(), fontsize=z_y_font)

            plt.suptitle(f"Cluster {cid}\nk={k} | n={n_cgenes} | {tag}", y=1.02, fontsize=14)
            plt.savefig(f"{args.data_dir}/cluster_heatmap_k{k}_cluster_{cid}_{tag}.png", dpi=200, bbox_inches='tight')
            plt.close()
        except Exception as e:
            print(f"  ❌ Failed to generate heatmap for cluster {cid}: {e}")
            plt.close('all')


# ─────────────────────────────────────────────
# SAVE METRICS
# ─────────────────────────────────────────────
metrics_df = pd.DataFrame(metrics)
metrics_df.to_csv(f"{args.data_dir}/clustering_metrics_{tag}.csv", index=False)

print("\nProcessing Complete.")
print(f"Metrics:\n{metrics_df}")
print(f"Results saved to: {args.data_dir}")

