#!/usr/bin/env python3

import pandas as pd
import argparse
import os
from datetime import datetime

# ──────────────────────────────────────────
# HELPERS
# ──────────────────────────────────────────
def safe_read(path):
    if os.path.exists(path):
        return pd.read_csv(path)
    return None

def fmt(x, nd=3):
    try:
        return f"{float(x):.{nd}g}"
    except Exception:
        return "NA"

def detect_category_column(df):
    """Try to infer enrichment source column name"""
    for col in ["category", "category_label", "source", "library", "annotation", "term_type"]:
        if col in df.columns:
            return col
    return None

# ──────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────
def main(data_dir, out_name):
    if not os.path.exists(data_dir):
        print(f"❌ Error: Data directory not found: {data_dir}")
        return

    # Load STRING outputs
    string_k = safe_read(os.path.join(data_dir, "string_k_summary.csv"))
    string_cat = safe_read(os.path.join(data_dir, "string_category_summary.csv"))

    # Load Enrichr outputs
    enrichr_k = safe_read(os.path.join(data_dir, "enrichr_k_summary.csv"))
    enrichr_terms = safe_read(os.path.join(data_dir, "enrichr_terms.csv"))
    
    if all(x is None for x in [string_k, string_cat, enrichr_k, enrichr_terms]):
        print(f"⚠ Warning: No enrichment data found in {data_dir}. Report will be empty.")

    out_path = os.path.join(data_dir, out_name)

    try:
        with open(out_path, "w") as f:
            f.write("# Combined STRING + Enrichr Enrichment Report\n\n")
            f.write(f"Generated: {datetime.now().isoformat(timespec='seconds')}\n\n")

            # ──────────────────────────────────────────
            # OVERVIEW
            # ──────────────────────────────────────────
            f.write("## Overview\n\n")
            f.write(
                "This report summarizes functional enrichment results from:\n"
                "- **STRING** (network-aware enrichment)\n"
                "- **Enrichr** (gene set–based enrichment)\n\n"
                "Results are reported across multiple clustering resolutions (k values).\n"
                "No single k is assumed to be optimal; trends and consistency are emphasized.\n\n"
            )

            # ──────────────────────────────────────────
            # STRING SUMMARY
            # ──────────────────────────────────────────
            if string_k is not None:
                f.write("## STRING enrichment summary\n\n")

                for _, row in string_k.iterrows():
                    f.write(
                        f"### k = {int(row['k'])}\n\n"
                        f"- Clusters: {int(row['n_clusters'])}\n"
                        f"- Median cluster size: {fmt(row['median_cluster_size'],1)} genes\n"
                        f"- Clusters enriched: {int(row['n_clusters_enriched'])} "
                        f"({fmt(row['fraction_enriched'],2)} fraction)\n"
                        f"- Total enriched terms: {int(row['total_enriched_terms'])}\n"
                        f"- Median best FDR: {fmt(row['median_best_fdr'])}\n\n"
                    )

                if string_cat is not None:
                    cat_col = detect_category_column(string_cat)

                    if cat_col is None:
                        f.write(
                            "⚠ STRING category breakdown unavailable "
                            "(no recognizable category column found).\n\n"
                        )
                    else:
                        f.write("### STRING enrichment sources\n\n")
                        for k, sub in string_cat.groupby("k"):
                            f.write(f"**k = {k}**\n\n")
                            for _, r in sub.sort_values("n_terms", ascending=False).iterrows():
                                f.write(
                                    f"- {r[cat_col]}: "
                                    f"{int(r['n_terms'])} terms "
                                    f"(best FDR {fmt(r['best_fdr'])})\n"
                                )
                            f.write("\n")

            # ──────────────────────────────────────────
            # ENRICHR SUMMARY
            # ──────────────────────────────────────────
            if enrichr_k is not None:
                f.write("## Enrichr enrichment summary\n\n")

                for _, row in enrichr_k.iterrows():
                    f.write(
                        f"### k = {int(row['k'])}\n\n"
                        f"- Clusters: {int(row['n_clusters'])}\n"
                        f"- Median cluster size: {fmt(row['median_cluster_size'],1)} genes\n"
                        f"- Clusters enriched: {int(row['n_clusters_enriched'])} "
                        f"({fmt(row['fraction_enriched'],2)} fraction)\n"
                        f"- Total enriched terms: {int(row['total_enriched_terms'])}\n"
                        f"- Median best FDR: {fmt(row['median_best_fdr'])}\n\n"
                    )

                if enrichr_terms is not None:
                    f.write("### Enrichr enrichment sources\n\n")
                    for k, sub in enrichr_terms.groupby("k"):
                        f.write(f"**k = {k}**\n\n")
                        for lib, n in (
                            sub.groupby("library").size().sort_values(ascending=False).items()
                        ):
                            f.write(f"- {lib}: {n} terms\n")
                        f.write("\n")

            # ──────────────────────────────────────────
            # INTERPRETATION
            # ──────────────────────────────────────────
            f.write("## Interpretation guidance\n\n")
            f.write(
                "- Lower k values tend to show broader, more generic enrichment.\n"
                "- Intermediate k values often reveal pathway-level or tissue-specific biology.\n"
                "- Higher k values produce more specific but sparser enrichment signals.\n\n"
                "Concordance between STRING and Enrichr strengthens biological confidence.\n"
                "Divergence may reflect literature bias, domain-level signal, or annotation gaps.\n\n"
            )

            f.write("---\n\nEnd of report.\n")

        print(f"\n✔ Combined report written to:\n{out_path}\n")

    except Exception as e:
        print(f"❌ Error writing report to {out_path}: {e}")

# ──────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a combined STRING + Enrichr enrichment report"
    )
    parser.add_argument(
        "--data-dir",
        required=True,
        help="Directory containing STRING and Enrichr output CSVs",
    )
    parser.add_argument(
        "--output-file",
        default="COMBINED_ENRICHMENT_REPORT.md",
        help="Output markdown report filename",
    )
    args = parser.parse_args()
    main(args.data_dir, args.output_file)
