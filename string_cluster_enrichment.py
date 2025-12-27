#!/usr/bin/env python3

import requests
import pandas as pd
import argparse
import glob
import os
import re
import numpy as np

STRING_API_URL = "https://string-db.org/api/json/enrichment"
SPECIES = 9606
FDR_CUTOFF = 0.05

# Regex to capture k and cluster id
FILENAME_PATTERNS = [
    re.compile(r"k(?P<k>\d+).*c(?P<cluster>\d+)", re.IGNORECASE),
    re.compile(r"cluster[_-]?k(?P<k>\d+).*cluster[_-]?(?P<cluster>\d+)", re.IGNORECASE),
]

def parse_k_cluster(filename):
    for pat in FILENAME_PATTERNS:
        m = pat.search(filename)
        if m:
            return f"k{m.group('k')}", f"c{m.group('cluster')}"
    return None, None

def query_string_enrichment(genes):
    r = requests.post(
        STRING_API_URL,
        data={
            "identifiers": "%0d".join(genes),
            "species": SPECIES,
        },
        timeout=60
    )
    r.raise_for_status()
    return r.json()

def load_genes(path):
    with open(path) as f:
        return [g.strip() for g in f if g.strip()]

def main(cluster_dir, out_prefix):
    cluster_files = glob.glob(os.path.join(cluster_dir, "*genes.txt"))

    all_cluster_stats = []
    all_terms = []

    for f in sorted(cluster_files):
        fname = os.path.basename(f)
        k_label, cluster_id = parse_k_cluster(fname)

        if k_label is None:
            print(f"⚠ Skipping (cannot parse k/cluster): {fname}")
            continue

        genes = load_genes(f)
        if len(genes) < 3:
            continue

        print(f"Processing {fname} ({k_label}, {cluster_id}, {len(genes)} genes)")

        try:
            results = query_string_enrichment(genes)
        except Exception as e:
            print(f"❌ STRING API error for {fname}: {e}")
            continue

        sig = [r for r in results if r.get("fdr", 1) <= FDR_CUTOFF]

        # Per-term records
        for r in sig:
            all_terms.append({
                "k": k_label,
                "cluster": cluster_id,
                "category": r["category"],
                "term": r["term"],
                "description": r.get("description", ""),
                "fdr": r["fdr"],
                "genes_in_term": r["number_of_genes"],
            })

        # Per-cluster summary
        all_cluster_stats.append({
            "k": k_label,
            "cluster": cluster_id,
            "n_genes": len(genes),
            "n_significant_terms": len(sig),
            "best_fdr": min([r["fdr"] for r in sig], default=np.nan),
        })

    # Build dataframes
    df_clusters = pd.DataFrame(all_cluster_stats)
    df_terms = pd.DataFrame(all_terms)

    # Aggregate per-k metrics
    k_summary = (
        df_clusters
        .groupby("k")
        .agg(
            n_clusters=("cluster", "count"),
            n_clusters_enriched=("n_significant_terms", lambda x: (x > 0).sum()),
            total_terms=("n_significant_terms", "sum"),
            median_best_fdr=("best_fdr", "median"),
        )
        .reset_index()
    )

    k_summary["fraction_enriched"] = (
        k_summary["n_clusters_enriched"] / k_summary["n_clusters"]
    )

    # Write outputs
    df_clusters.to_csv(f"{out_prefix}_per_cluster.csv", index=False)
    df_terms.to_csv(f"{out_prefix}_terms.csv", index=False)
    k_summary.to_csv(f"{out_prefix}_k_summary.csv", index=False)

    print("\n✔ STRING enrichment complete")
    print(f"→ {out_prefix}_k_summary.csv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run STRING enrichment on cluster gene lists in a flat directory"
    )
    parser.add_argument(
        "--cluster-dir",
        required=True,
        help="Directory containing all *_genes.txt files for all k values"
    )
    parser.add_argument(
        "--out-prefix",
        default="string_multi_k",
        help="Output prefix"
    )
    args = parser.parse_args()

    main(args.cluster_dir, args.out_prefix)
