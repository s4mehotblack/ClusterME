#!/usr/bin/env python3

import argparse
import os
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(
        description="Visualise per-gene and per-cluster silhouette scores from GTEx clustering results"
    )
    parser.add_argument(
        "--data-dir",
        required=True,
        help="Directory containing gene_silhouette_scores_k*_*.csv files"
    )
    parser.add_argument(
        "--output-dir",
        help="Directory to save plots (defaults to data-dir)"
    )
    return parser.parse_args()

def plot_silhouette(df, k, tag, out_dir):
    """
    Generates a silhouette plot (per gene) and a summary bar chart (per cluster) for a specific k.
    """
    # Create a figure with 2 subplots
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(18, 8)

    # ---------------------------------------------------------
    # PLOT 1: The Silhouette Plot (Per Gene)
    # ---------------------------------------------------------
    # The silhouette coefficient ranges from -1 to 1
    ax1.set_xlim([-0.2, 1])
    
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    n_genes = len(df)
    ax1.set_ylim([0, n_genes + (k + 1) * 10])

    cluster_labels = sorted(df["cluster"].unique())
    
    y_lower = 10
    cluster_avg_scores = []
    
    # Match color palette used in clustering script
    if k <= 10:
        palette = sns.color_palette("tab10", k)
    else:
        palette = sns.color_palette("husl", k)

    for i, color in zip(cluster_labels, palette):
        # Aggregate the silhouette scores for samples belonging to cluster i, and sort them
        ith_cluster_silhouette_values = df[df["cluster"] == i]["silhouette_score"].values
        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i
        
        # Calculate statistics
        avg_score = np.mean(ith_cluster_silhouette_values)
        cluster_avg_scores.append({
            'cluster': i, 
            'avg_silhouette': avg_score, 
            'size': size_cluster_i,
            'color': color
        })

        # Fill the area for this cluster
        ax1.fill_betweenx(
            np.arange(y_lower, y_upper),
            0,
            ith_cluster_silhouette_values,
            facecolor=color,
            edgecolor=color,
            alpha=0.7,
        )

        # Label the silhouette plots with their cluster numbers at the middle
        ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    avg_score_total = df["silhouette_score"].mean()
    
    ax1.set_title(f"Silhouette Plot for k={k}\n(Each bar represents one gene)")
    ax1.set_xlabel("Silhouette Coefficient Values")
    ax1.set_ylabel("Cluster Label")

    ax1.set_yticks([])  # Clear the yaxis labels / ticks
    ax1.set_xticks([-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])

    # ---------------------------------------------------------
    # PLOT 2: Average Score per Cluster (Summary)
    # ---------------------------------------------------------
    avg_df = pd.DataFrame(cluster_avg_scores)
    
    sns.barplot(
        data=avg_df, 
        x='cluster', 
        y='avg_silhouette', 
        ax=ax2, 
        palette=palette
    )
    
    ax2.set_title(f"Average Silhouette Score per Cluster\n(Higher is better)")
    ax2.set_ylabel("Mean Silhouette Coefficient")
    ax2.set_xlabel("Cluster ID")
    
    # Adjust y-limit to show negative values if they exist, or at least 0-1
    y_min = min(0, avg_df['avg_silhouette'].min()) - 0.05
    ax2.set_ylim([y_min, 1.0])
    
    # Add value labels on bars
    for index, row in avg_df.iterrows():
        ax2.text(
            index, 
            row.avg_silhouette + 0.01, 
            f"{row.avg_silhouette:.2f}", 
            color='black', 
            ha="center"
        )
        # Add size label inside/below bar
        # ax2.text(index, 0.02, f"n={row['size']}", color='black', ha="center", fontsize=9)

    plt.suptitle(f"Silhouette Analysis (k={k}) | {tag}", fontsize=14, fontweight='bold', y=1.02)
    
    outfile = os.path.join(out_dir, f"silhouette_analysis_k{k}_{tag}.png")
    plt.tight_layout()
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Generated plot: {outfile}")

def main():
    args = parse_args()
    data_dir = args.data_dir
    output_dir = args.output_dir if args.output_dir else data_dir
    
    if not os.path.exists(data_dir):
        print(f"Error: Data directory '{data_dir}' does not exist.")
        return

    os.makedirs(output_dir, exist_ok=True)
    
    # Find files
    pattern = os.path.join(data_dir, "gene_silhouette_scores_k*_*.csv")
    files = glob.glob(pattern)
    
    if not files:
        print(f"No silhouette score files found in {data_dir}")
        print("Tip: Run gtex_cluster_v2.py first to generate clustering results.")
        return

    print(f"Found {len(files)} silhouette score files to process.")
    
    for f in sorted(files):
        try:
            # Extract info from filename
            basename = os.path.basename(f)
            # Expected format: gene_silhouette_scores_k{k}_{tag}.csv
            # Example: gene_silhouette_scores_k7_median_log2_zscore.csv
            
            parts = basename.replace("gene_silhouette_scores_", "").replace(".csv", "")
            
            # Extract k (everything before the first underscore)
            k_part = parts.split('_')[0] 
            
            if not k_part.startswith('k') or not k_part[1:].isdigit():
                print(f"Skipping file with unexpected format: {basename}")
                continue
                
            k = int(k_part.replace('k', ''))
            
            # Extract tag (everything after k_)
            tag = parts[len(k_part)+1:]
            
            print(f"Processing k={k}...")
            
            try:
                df = pd.read_csv(f)
                if df.empty:
                    print(f"  ⚠ Skipping {basename}: File is empty.")
                    continue
                if "silhouette_score" not in df.columns or "cluster" not in df.columns:
                    print(f"  ❌ Error in {basename}: Missing required columns ('silhouette_score', 'cluster').")
                    continue
                    
                plot_silhouette(df, k, tag, output_dir)
            except Exception as e:
                print(f"  ❌ Failed to plot {basename}: {e}")
                plt.close('all')
            
        except Exception as e:
            print(f"Error processing {f}: {e}")

    print("\nDone.")

if __name__ == "__main__":
    main()
