#!/usr/bin/env python3

import requests
import pandas as pd
import argparse
import glob
import os
import re
import time
import numpy as np

# ──────────────────────────────────────────
# ENRICHR API ENDPOINTS (OFFICIAL)
# ──────────────────────────────────────────
ADD_LIST_URL = "https://maayanlab.cloud/Enrichr/addList"
ENRICH_URL = "https://maayanlab.cloud/Enrichr/enrich"

# Recommended libraries (balanced sensitivity & interpretability)
LIBRARIES = [
    "GO_Biological_Process_2025",
    "GO_Molecular_Function_2025",
    "GO_Cellular_Component_2025",
    "KEGG_2025_Human",
    "Reactome_2022",
    "WikiPathways_2023_Human",
    "MSigDB_Hallmark_2020",
    "Human_Phenotype_Ontology",
    "Pfam_2019",
    "InterPro_2019",
    "SMART_Domains_2018",
    "OMIM_Disease",
]

# Analysis parameters
FDR_CUTOFF = 0.05
MIN_GENES = 3
REQUEST_DELAY = 0.5  # seconds between requests (polite)

CLUSTER_RE = re.compile(r"k(?P<k>\d+).*cluster_(?P<c>\d+)", re.I)

# ──────────────────────────────────────────
# HELPERS
# ──────────────────────────────────────────
def parse_k_cluster(fname):
    m = CLUSTER_RE.search(fname)
    if not m:
        return None, None
    return int(m.group("k")), int(m.group("c"))

def load_genes(path):
    with open(path) as f:
        genes = [l.strip() for l in f if l.strip()]
    return sorted(set(genes))

def submit_gene_list(genes):
    gene_text = "\n".join(genes)

    r = requests.post(
        ADD_LIST_URL,
        files={"list": ("genes.txt", gene_text)},
        timeout=60,
    )

    if not r.ok:
        raise RuntimeError(f"Enrichr addList failed: HTTP {r.status_code}")

    data = r.json()
    if "userListId" not in data:
        raise RuntimeError("No userListId returned by Enrichr")

    return data["userListId"]

def fetch_enrichment(user_list_id, library):
    r = requests.get(
        ENRICH_URL,
        params={
            "userListId": user_list_id,
            "backgroundType": library,
        },
        timeout=60,
    )

    if not r.ok:
        return []

    return r.json().get(library, [])

# ──────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────
def main(cluster_dir):
    if not os.path.exists(cluster_dir):
        print(f"❌ Error: Data directory not found: {cluster_dir}")
        return

    files = sorted(glob.glob(os.path.join(cluster_dir, "*genes.txt")))
    if not files:
        print(f"⚠ Warning: No cluster gene files (*genes.txt) found in {cluster_dir}")
        return

    print(f"Found {len(files)} cluster files")

    term_rows = []
    cluster_rows = []

    for path in files:
        fname = os.path.basename(path)
        k, cluster = parse_k_cluster(fname)
        if k is None:
            continue

        genes = load_genes(path)
        if len(genes) < MIN_GENES:
            continue

        print(f"→ k={k} cluster={cluster} ({len(genes)} genes)")

        try:
            user_list_id = submit_gene_list(genes)
        except Exception as e:
            print(f"❌ Enrichr submission failed: {e}")
            continue

        total_terms = 0
        best_fdr = np.nan

        for lib in LIBRARIES:
            results = fetch_enrichment(user_list_id, lib)

            for r in results:
                # Defensive parsing: Enrichr rows vary in length
                if len(r) < 7:
                    continue

                term = r[1]
                p_val = r[2]
                z_score = r[3]
                combined = r[4]
                adj_p = r[6]

                try:
                    adj_p = float(adj_p)
                except (TypeError, ValueError):
                    continue

                if adj_p <= FDR_CUTOFF:
                    term_rows.append({
                        "k": k,
                        "cluster": cluster,
                        "library": lib,
                        "term": term,
                        "p_value": p_val,
                        "adjusted_p": adj_p,
                        "z_score": z_score,
                        "combined_score": combined,
                        "n_genes_cluster": len(genes),
                    })
                    total_terms += 1
                    best_fdr = min(best_fdr, adj_p) if not np.isnan(best_fdr) else adj_p

        cluster_rows.append({
            "k": k,
            "cluster": cluster,
            "n_genes": len(genes),
            "n_enriched_terms": total_terms,
            "best_fdr": best_fdr,
        })

        time.sleep(REQUEST_DELAY)

    # ──────────────────────────────────────────
    # OUTPUT
    # ──────────────────────────────────────────
    df_terms = pd.DataFrame(term_rows)
    df_clusters = pd.DataFrame(cluster_rows)

    k_summary = (
        df_clusters
        .groupby("k")
        .agg(
            n_clusters=("cluster", "count"),
            median_cluster_size=("n_genes", "median"),
            n_clusters_enriched=("n_enriched_terms", lambda x: (x > 0).sum()),
            fraction_enriched=("n_enriched_terms", lambda x: (x > 0).mean()),
            total_enriched_terms=("n_enriched_terms", "sum"),
            median_best_fdr=("best_fdr", "median"),
        )
        .reset_index()
        .sort_values("k")
    )

    df_terms.to_csv(os.path.join(cluster_dir, "enrichr_terms.csv"), index=False)
    df_clusters.to_csv(os.path.join(cluster_dir, "enrichr_per_cluster.csv"), index=False)
    k_summary.to_csv(os.path.join(cluster_dir, "enrichr_k_summary.csv"), index=False)

    print("\n✔ Enrichr analysis complete")
    print(f"Results written to: {cluster_dir}")

# ──────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run Enrichr enrichment for GTEx gene clusters"
    )
    parser.add_argument(
        "--data-dir",
        required=True,
        help="Directory containing cluster_*_genes.txt files",
    )
    args = parser.parse_args()
    main(args.data_dir)
