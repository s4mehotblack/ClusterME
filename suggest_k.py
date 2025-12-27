#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
from sklearn.metrics import silhouette_score
import os

def parse_args():
    parser = argparse.ArgumentParser(
        description="Suggest optimal number of clusters for GTEx-style expression data"
    )
    
    parser.add_argument(
        "--expression-file",
        required=True,
        help="All-tissues expression CSV (gene,tissue,mean_tpm,median_tpm,n_samples)"
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
    
    parser.add_argument("--kmin", type=int, default=2, help="Minimum k (default=2)")
    parser.add_argument("--kmax", type=int, default=15, help="Maximum k (default=15)")
    
    return parser.parse_args()

def calculate_sse(matrix_values, labels, k):
    """Calculate Sum of Squared Errors (distortion) for a given clustering."""
    sse = 0
    for i in range(1, k + 1):
        # Extract points in this cluster
        points = matrix_values[labels == i]
        if len(points) > 0:
            centroid = points.mean(axis=0)
            sse += ((points - centroid) ** 2).sum()
    return sse

def main():
    args = parse_args()

    # ─────────────────────────────────────────────
    # LOAD & PREPARE MATRIX (Identical to gtex_cluster_v2.py)
    # ─────────────────────────────────────────────
    if not os.path.exists(args.expression_file):
        print(f"❌ Error: Expression file not found: {args.expression_file}")
        return

    print(f"Loading expression data: {args.expression_file}")
    
    try:
        expr = pd.read_csv(args.expression_file)
    except Exception as e:
        print(f"❌ Error reading expression file: {e}")
        return

    value_col = f"{args.stat}_tpm"

    # Pivot to (Genes x Tissues)
    matrix = (
        expr
        .pivot(index="gene", columns="tissue", values=value_col)
        .fillna(0)
    )

    # 1. Log Transform
    if not args.no_log2:
        matrix = np.log2(matrix + 1)

    # 2. Z-Score Standardization (Row-wise)
    if args.z_score:
        print("Applying Z-score standardization (row-wise)...")
        # subtract mean, divide by std (add small epsilon to avoid div/0)
        matrix = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1) + 1e-8, axis=0)

    print(f"Matrix shape: {matrix.shape[0]} genes × {matrix.shape[1]} tissues")
    
    # ─────────────────────────────────────────────
    # HIERARCHICAL CLUSTERING
    # ─────────────────────────────────────────────
    print("Computing hierarchical clustering (Ward)...")
    
    # Use scipy to ensure exact match with pipeline
    dist_matrix = pdist(matrix.values)
    row_linkage = linkage(dist_matrix, method="ward")

    k_range = range(args.kmin, args.kmax + 1)
    silhouettes = []
    sse_values = []

    print("\nEvaluating k values...\n")

    for k in k_range:
        # Cut the linkage tree to get labels
        cluster_ids = fcluster(
            row_linkage,
            k,
            criterion="maxclust"
        )
        
        # Silhouette Score
        sil = silhouette_score(matrix, cluster_ids)
        silhouettes.append(sil)
        
        # SSE Calculation for Elbow Method
        sse = calculate_sse(matrix.values, cluster_ids, k)
        sse_values.append(sse)

        print(f"k={k:2d} | silhouette={sil:.3f} | SSE={sse:.1f}")

    # ─────────────────────────────────────────────
    # RECOMMENDATION LOGIC
    # ─────────────────────────────────────────────
    best_sil_k = k_range[np.argmax(silhouettes)]

    # Elbow heuristic: find point of maximum curvature (simplistic)
    # Using 2nd derivative approximation
    if len(sse_values) > 2:
        deltas = np.diff(sse_values)
        # Shift index by +1 because diff reduces length by 1, and +1 again for 1-based k-range start? 
        # Actually easier to map back to k_range
        # acceleration = diff(diff(sse))
        acceleration = np.diff(deltas)
        # max acceleration often indicates the elbow in an SSE curve (which is decreasing and convex)
        elbow_idx = np.argmax(acceleration) + 1 
        elbow_k = k_range[elbow_idx]
    else:
        elbow_k = k_range[0]

    # Conservative biological compromise
    recommended_k = int(round(np.median([best_sil_k, elbow_k])))

    print("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    print("📊 CLUSTER NUMBER SUGGESTION")
    print("━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    print(f"Best silhouette score:  k = {best_sil_k}")
    print(f"Elbow point estimate:  k = {elbow_k}")
    print(f"\n✔ Recommended starting k: {recommended_k}")
    print("━━━━━━━━━━━━━━━━━━━━━━━━━━━━")

if __name__ == "__main__":
    main()
