"""
STRING cluster enrichment reporting utility

Inputs:
- A directory containing gene cluster files, e.g.:
    cluster_k6_c3_genes.txt
    cluster_k8_c1_genes.txt

Outputs (written to same directory):
- string_per_cluster.csv
- string_terms.csv
- string_k_summary.csv
- string_category_summary.csv
- STRING_REPORT.md
"""

import requests
import pandas as pd
import argparse
import glob
import os
import re
import numpy as np
from collections import defaultdict

# ──────────────────────────────────────────
# CONFIG
# ──────────────────────────────────────────
STRING_API_URL = "https://string-db.org/api/json/enrichment"
SPECIES = 9606
FDR_CUTOFF = 0.05
MIN_GENES = 3

CATEGORY_LABELS = {
    "Process": "GO Biological Process",
    "Function": "GO Molecular Function",
    "Component": "GO Cellular Component",
    "KEGG": "KEGG Pathways",
    "Reactome": "Reactome Pathways",
    "WikiPathways": "WikiPathways",
    "UniProt Keywords": "UniProt Keywords",
    "PFAM": "Pfam domains",
    "INTERPRO": "InterPro domains",
    "SMART": "SMART domains",
    "PubMed": "PubMed publications",
}

FILENAME_PATTERNS = [
    re.compile(r"k(?P<k>\d+).*c(?P<cluster>\d+)", re.IGNORECASE),
    re.compile(r"cluster[_-]?k(?P<k>\d+).*cluster[_-]?(?P<cluster>\d+)", re.IGNORECASE),
]

# ──────────────────────────────────────────
# HELPERS
# ──────────────────────────────────────────
def parse_k_cluster(filename):
    for pat in FILENAME_PATTERNS:
        m = pat.search(filename)
        if m:
            return int(m.group("k")), int(m.group("cluster"))
    return None, None

def load_genes(path):
    with open(path) as f:
        return [g.strip() for g in f if g.strip()]

def query_string(genes):
    r = requests.post(
        STRING_API_URL,
        data={
            "identifiers": "%0d".join(genes),
            "species": SPECIES,
        },
        timeout=60,
    )
    r.raise_for_status()
    return r.json()

# ──────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────
def main(cluster_dir):
    if not os.path.exists(cluster_dir):
        print(f"❌ Error: Data directory not found: {cluster_dir}")
        return

    cluster_files = sorted(glob.glob(os.path.join(cluster_dir, "*genes.txt")))

    if not cluster_files:
        print(f"⚠ Warning: No gene cluster files (*genes.txt) found in {cluster_dir}")
        return

    cluster_rows = []
    term_rows = []

    print(f"Found {len(cluster_files)} cluster files")

    for path in cluster_files:
        fname = os.path.basename(path)
        k, cluster = parse_k_cluster(fname)
        if k is None:
            print(f"Skipping unrecognized filename: {fname}")
            continue

        genes = load_genes(path)
        if len(genes) < MIN_GENES:
            continue

        try:
            results = query_string(genes)
        except Exception as e:
            print(f"STRING query failed for {fname}: {e}")
            continue

        sig = [r for r in results if r.get("fdr", 1.0) <= FDR_CUTOFF]

        # Per-term rows
        for r in sig:
            term_rows.append({
                "k": k,
                "cluster": cluster,
                "n_genes": len(genes),
                "category": r.get("category"),
                "category_label": CATEGORY_LABELS.get(r.get("category"), r.get("category")),
                "term": r.get("term"),
                "description": r.get("description", ""),
                "fdr": r.get("fdr"),
            })

        # Per-cluster summary
        cluster_rows.append({
            "k": k,
            "cluster": cluster,
            "n_genes": len(genes),
            "n_terms": len(sig),
            "best_fdr": min([r["fdr"] for r in sig], default=np.nan),
        })

    df_clusters = pd.DataFrame(cluster_rows)
    df_terms = pd.DataFrame(term_rows)

    # ──────────────────────────────────────────
    # PER-k SUMMARY
    # ──────────────────────────────────────────
    k_summary = (
        df_clusters
        .groupby("k")
        .agg(
            n_clusters=("cluster", "count"),
            median_cluster_size=("n_genes", "median"),
            n_clusters_enriched=("n_terms", lambda x: (x > 0).sum()),
            fraction_enriched=("n_terms", lambda x: (x > 0).mean()),
            total_enriched_terms=("n_terms", "sum"),
            median_best_fdr=("best_fdr", "median"),
        )
        .reset_index()
        .sort_values("k")
    )

    # ──────────────────────────────────────────
    # CATEGORY BREAKDOWN
    # ──────────────────────────────────────────
    category_summary = (
        df_terms
        .groupby(["k", "category_label"])
        .agg(
            n_terms=("term", "count"),
            n_clusters=("cluster", "nunique"),
            best_fdr=("fdr", "min"),
        )
        .reset_index()
        .sort_values(["k", "n_terms"], ascending=[True, False])
    )

    # ──────────────────────────────────────────
    # WRITE OUTPUTS
    # ──────────────────────────────────────────
    df_clusters.to_csv(os.path.join(cluster_dir, "string_per_cluster.csv"), index=False)
    df_terms.to_csv(os.path.join(cluster_dir, "string_terms.csv"), index=False)
    k_summary.to_csv(os.path.join(cluster_dir, "string_k_summary.csv"), index=False)
    category_summary.to_csv(os.path.join(cluster_dir, "string_category_summary.csv"), index=False)

    # ──────────────────────────────────────────
    # HUMAN-READABLE REPORT
    # ──────────────────────────────────────────
    report_path = os.path.join(cluster_dir, "STRING_REPORT.md")

    with open(report_path, "w") as f:
        f.write("# STRING Enrichment Report\n\n")

        f.write(
            "This report summarizes functional enrichment of gene clusters using the STRING API.\n\n"
            "Important interpretation notes:\n"
            "- Lower k values favor enrichment due to larger cluster size.\n"
            "- Increasing k may *increase the number of enriched clusters* as distinct biology separates.\n"
            "- Loss of enrichment at high k reflects statistical underpowering, not absence of biology.\n\n"
        )

        f.write("## Per-k enrichment overview\n\n")
        f.write(k_summary.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Enrichment source breakdown\n\n")
        f.write(
            "The table below shows how enrichment is distributed across annotation sources.\n"
            "This helps distinguish pathway-driven, functional, structural, or literature-driven signal.\n\n"
        )
        f.write(category_summary.to_markdown(index=False))
        f.write("\n\n")

        f.write("## How to use this report\n\n")
        f.write(
            "- Look for k values where enrichment spreads across *multiple clusters*.\n"
            "- Prefer k where GO / pathway categories dominate over domain-only or PubMed-only signal.\n"
            "- Use heatmaps alongside this report to confirm tissue coherence.\n"
        )

    print("\n✔ STRING enrichment analysis complete")
    print(f"→ Results written to: {cluster_dir}")
    print(f"→ Report: {report_path}")

# ──────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STRING enrichment for gene clusters")
    parser.add_argument("--data-dir", required=True, help="Directory containing gene cluster files")
    args = parser.parse_args()

    main(args.data_dir)
