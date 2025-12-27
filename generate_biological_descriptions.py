#!/usr/bin/env python3

"""
Generate human-readable biological summaries for GTEx expression clusters.

Inputs:
- data_dir containing:
    - cluster_kX_cluster_Y_genes.txt
    - string_terms.csv
    - enrichr_terms.csv
- expression file (gene, tissue, mean_tpm / median_tpm)

Output:
- CLUSTER_BIOLOGICAL_DESCRIPTIONS.md (written inside data_dir)
"""

import os
import re
import pandas as pd
import requests
from collections import defaultdict

# -----------------------------
# Parameters
# -----------------------------

FDR_CUTOFF = 0.05
ENRICHR_COMBINED_SCORE_MIN = 10
TOP_TERMS_PER_SOURCE = 5
TOP_TISSUES = 5
TOP_GENES = 20

OLS_API = "https://www.ebi.ac.uk/ols/api/ontologies/{onto}/terms"
ONTOLOGY_CACHE = {}

# -----------------------------
# Utilities
# -----------------------------

def resolve_ontology_term(term_id):
    """
    Resolve EFO / HP / DOID IDs to human-readable labels using OLS API.
    Uses a global cache to avoid redundant network calls.
    """
    if term_id in ONTOLOGY_CACHE:
        return ONTOLOGY_CACHE[term_id]

    if term_id.startswith("EFO:"):
        onto = "efo"
    elif term_id.startswith("HP:"):
        onto = "hp"
    elif term_id.startswith("DOID:"):
        onto = "doid"
    else:
        return None

    url = OLS_API.format(onto=onto)
    params = {"obo_id": term_id}

    try:
        r = requests.get(url, params=params, timeout=10)
        r.raise_for_status()
        data = r.json()

        terms = data.get("_embedded", {}).get("terms", [])
        if terms:
            label = terms[0].get("label")
            ONTOLOGY_CACHE[term_id] = label
            return label
    except Exception:
        pass

    return None


def parse_run_config(data_dir):
    """Parse run_config.cfg into a dictionary."""
    config = {}
    config_path = os.path.join(data_dir, "run_config.cfg")
    if not os.path.exists(config_path):
        return None
    
    pattern = re.compile(r'^(\w+)="?([^"]*)"?')
    with open(config_path, 'r') as f:
        for line in f:
            m = pattern.match(line.strip())
            if m:
                config[m.group(1)] = m.group(2)
    return config


def infer_tag(data_dir):
    """Try to infer the naming tag (e.g. median_log2_zscore) from filenames."""
    for fname in os.listdir(data_dir):
        if fname.startswith("expression_matrix_") and fname.endswith(".csv"):
            return fname.replace("expression_matrix_", "").replace(".csv", "")
    return "results"


def format_p(val):
    """Format p-values for readability."""
    if val < 0.001:
        return f"{val:.2e}"
    else:
        return f"{val:.3f}"


def generate_html_report(data_dir, report_data, config, tag, report_title):
    """Generate a comprehensive HTML report with hyperlinks."""
    html_path = os.path.join(data_dir, "CLUSTER_BIOLOGICAL_REPORT.html")
    
    css = """
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif; line-height: 1.6; color: #24292e; max-width: 1000px; margin: 0 auto; padding: 20px; background-color: #f6f8fa; }
    .card { background: white; border: 1px solid #e1e4e8; border-radius: 6px; padding: 20px; margin-bottom: 20px; box-shadow: 0 1px 3px rgba(27,31,35,0.12); }
    h1, h2, h3 { border-bottom: 1px solid #eaecef; padding-bottom: 0.3em; margin-top: 24px; }
    .tag { display: inline-block; padding: 2px 8px; font-size: 12px; font-weight: 600; line-height: 1; border-radius: 10px; background-color: #e1e4e8; margin-right: 5px; }
    .gene-list { font-family: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, "Liberation Mono", monospace; font-size: 12px; background: #f6f8fa; padding: 10px; border-radius: 6px; border: 1px solid #e1e4e8; overflow-x: auto; }
    .links { margin-top: 10px; font-size: 14px; }
    .links a { margin-right: 15px; color: #0366d6; text-decoration: none; }
    .links a:hover { text-decoration: underline; }
    .stats { color: #586069; font-size: 14px; }
    .config-table { width: 100%; border-collapse: collapse; font-size: 14px; }
    .config-table th, .config-table td { text-align: left; padding: 8px; border-bottom: 1px solid #e1e4e8; }
    .config-table th { background: #f6f8fa; width: 30%; }
    """

    html = [f"<!DOCTYPE html><html><head><meta charset='utf-8'><title>{report_title}</title><style>{css}</style></head><body>"]
    
    html.append(f"<h1>{report_title}</h1>")
    
    # Run Configuration
    if config:
        html.append("<div class='card'><h2>Run Configuration</h2><table class='config-table'>")
        for key in ['RUN_NAME', 'K_VALUES', 'STAT_METHOD', 'Z_FLAG', 'SORT_FLAG']:
            if key in config:
                html.append(f"<tr><th>{key}</th><td>{config[key] if config[key] else 'None'}</td></tr>")
        html.append("</table>")
        
        html.append("<h3>Report Thresholds</h3><table class='config-table'>")
        html.append(f"<tr><th>FDR Cutoff</th><td>{FDR_CUTOFF}</td></tr>")
        html.append(f"<tr><th>Enrichr Combined Score Min</th><td>{ENRICHR_COMBINED_SCORE_MIN}</td></tr>")
        html.append(f"<tr><th>Max Terms per Category</th><td>{TOP_TERMS_PER_SOURCE}</td></tr>")
        html.append(f"<tr><th>Max Tissues Shown</th><td>{TOP_TISSUES}</td></tr>")
        html.append(f"<tr><th>Max Genes Shown</th><td>{TOP_GENES}</td></tr>")
        html.append("</table></div>")

    # Table of Contents
    html.append("<div class='card'><h2>Table of Contents</h2><ul>")
    for k_val in sorted(report_data.keys()):
        html.append(f"<li><strong><a href='#res_{k_val}'>Clustering k = {k_val}</a></strong><ul>")
        html.append(f"<li><a href='#viz_{k_val}'>Resolution Visualizations</a></li>")
        for cluster_id in sorted(report_data[k_val].keys()):
            html.append(f"<li><a href='#cluster_{k_val}_{cluster_id}'>Cluster {cluster_id}</a></li>")
        html.append("</ul></li>")
    html.append("</ul>")
    
    # Global Files for this tag
    html.append("<div class='links'><strong>Global Files:</strong> ")
    html.append(f"<a href='expression_matrix_{tag}.csv' target='_blank'>Expression Matrix</a>")
    html.append(f"<a href='clustering_metrics_{tag}.csv' target='_blank'>Clustering Metrics</a>")
    html.append(f"<a href='string_terms.csv' target='_blank'>STRING Terms</a>")
    html.append(f"<a href='enrichr_terms.csv' target='_blank'>Enrichr Terms</a>")
    html.append("</div></div>")

    # Clusters
    for k_val in sorted(report_data.keys()):
        html.append(f"<h2 id='res_{k_val}'>Clustering with k = {k_val}</h2>")
        
        # Link to k-specific visualization
        html.append(f"<div class='card' id='viz_{k_val}'><h3>Resolution Visualizations (k={k_val})</h3><div class='links'>")
        html.append(f"<a href='global_heatmap_k{k_val}_{tag}.png' target='_blank'>Global Heatmap</a>")
        html.append(f"<a href='silhouette_analysis_k{k_val}_{tag}.png' target='_blank'>Silhouette Analysis Plot</a>")
        html.append(f"<a href='gene_clusters_k{k_val}_{tag}.csv' target='_blank'>Cluster Assignments (CSV)</a>")
        html.append(f"<a href='gene_silhouette_scores_k{k_val}_{tag}.csv' target='_blank'>Gene Silhouette Scores (CSV)</a>")
        html.append("</div></div>")

        for cluster_id, data in sorted(report_data[k_val].items()):
            html.append(f"<div class='card' id='cluster_{k_val}_{cluster_id}'>")
            html.append(f"<h3>Cluster {cluster_id} <span class='stats'>(k={k_val})</span></h3>")
            
            # Cluster Links
            html.append("<div class='links'>")
            html.append(f"<a href='cluster_k{k_val}_cluster_{cluster_id}_genes.txt' target='_blank'>Gene List (TXT)</a>")
            html.append(f"<a href='cluster_heatmap_k{k_val}_cluster_{cluster_id}_{tag}.png' target='_blank'>Cluster Heatmap (PNG)</a>")
            html.append("</div>")

            # Description
            desc_html = data['description'].replace('\n', '<br>')
            desc_html = re.sub(r'\*\*(.*?)\*\*', r'<strong>\1</strong>', desc_html)
            html.append(f"<p>{desc_html}</p>")
            
            # Top Genes
            html.append(f"<p><strong>Top Genes (by silhouette):</strong></p>")
            
            linked_genes = []
            for g, s in data['genes_with_scores']:
                link = f"<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene={g}' target='_blank'>{g}</a>"
                if s is not None:
                    linked_genes.append(f"{link} (sil: {s:.2f})")
                else:
                    linked_genes.append(link)
            
            html.append(f"<div class='gene-list'>{', '.join(linked_genes)} ...</div>")
            
            html.append("</div>")

    html.append("</body></html>")
    
    with open(html_path, 'w') as f:
        f.write("\n".join(html))
    print(f"✓ HTML report written to {html_path}")


def load_silhouette_scores(data_dir):
    """
    Load silhouette scores from gene_silhouette_scores_k*_*.csv files.
    Returns: dict[(k, gene)] -> score
    """
    scores = {}
    pattern = re.compile(r"gene_silhouette_scores_k(\d+)_.*\.csv")
    
    for fname in os.listdir(data_dir):
        m = pattern.match(fname)
        if not m:
            continue
            
        k = int(m.group(1))
        try:
            df = pd.read_csv(os.path.join(data_dir, fname))
            if "gene" in df.columns and "silhouette_score" in df.columns:
                for _, row in df.iterrows():
                    scores[(k, row["gene"])] = row["silhouette_score"]
        except Exception as e:
            print(f"Warning: Failed to load silhouette scores from {fname}: {e}")
            
    return scores


def load_gene_clusters(data_dir):
    clusters = defaultdict(dict)
    pattern = re.compile(r"cluster_k(\d+)_cluster_(\d+)_genes\.txt")

    for fname in os.listdir(data_dir):
        m = pattern.match(fname)
        if not m:
            continue

        k = int(m.group(1))
        cluster = int(m.group(2))

        with open(os.path.join(data_dir, fname)) as f:
            genes = [g.strip() for g in f if g.strip()]

        clusters[k][cluster] = genes

    return clusters


def summarize_tissues(expr_df, genes):
    sub = expr_df[expr_df["gene"].isin(genes)]
    if sub.empty:
        return []

    tissue_series = (
        sub.groupby("tissue", as_index=True)["value"]
        .mean()
        .sort_values(ascending=False)
    )

    return list(tissue_series.head(TOP_TISSUES).index)


def summarize_string_enrichment(df, k, cluster):
    sub = df[(df["k"] == k) & (df["cluster"] == cluster)]
    if sub.empty or "fdr" not in sub.columns:
        return {}, {}

    summaries = defaultdict(list)
    total_counts = defaultdict(int)

    for _, r in sub.iterrows():
        if r["fdr"] <= FDR_CUTOFF:
            total_counts[r["category"]] += 1
            # Capture description if available
            desc = r.get("description", "")
            summaries[r["category"]].append((r["term"], r["fdr"], desc))

    for cat in summaries:
        summaries[cat] = sorted(summaries[cat], key=lambda x: x[1])[:TOP_TERMS_PER_SOURCE]

    return summaries, total_counts


def summarize_enrichr_enrichment(df, k, cluster):
    sub = df[(df["k"] == k) & (df["cluster"] == cluster)]
    if sub.empty:
        return {}, {}

    if not {"adjusted_p", "combined_score"}.issubset(sub.columns):
        return {}, {}

    summaries = defaultdict(list)
    total_counts = defaultdict(int)

    for _, r in sub.iterrows():
        if r["adjusted_p"] <= FDR_CUTOFF and r["combined_score"] >= ENRICHR_COMBINED_SCORE_MIN:
            total_counts[r["library"]] += 1
            summaries[r["library"]].append((r["term"], r["adjusted_p"], r["combined_score"]))

    for lib in summaries:
        # Sort by Combined Score (index 2) descending
        summaries[lib] = sorted(summaries[lib], key=lambda x: x[2], reverse=True)[:TOP_TERMS_PER_SOURCE]

    return summaries, total_counts


def speculative_description(genes, tissues, string_terms, enrichr_terms, avg_sil=None, string_counts=None, enrichr_counts=None):
    lines = []

    lines.append(f"This cluster contains **{len(genes)} genes**.")
    
    if avg_sil is not None:
        lines.append(f"Average Silhouette Score: **{avg_sil:.3f}**.")

    if tissues:
        lines.append(
            f"Highest expression is observed in: **{', '.join(tissues)}**."
        )
        if len(tissues) >= TOP_TISSUES:
            lines.append(f"*(Showing top {TOP_TISSUES} tissues. Refer to full expression matrix for all values.)*")

    if string_terms:
        lines.append("\n**STRING network enrichment:**")
        lines.append(f"*(Thresholds: FDR < {FDR_CUTOFF}. Showing top {TOP_TERMS_PER_SOURCE} terms per category.)*")
        for cat, terms in string_terms.items():
            # terms is list of (term_id, fdr, description)
            total = string_counts.get(cat, 0) if string_counts else len(terms)
            formatted_terms = []
            for tid, p, desc in terms:
                # Use provided description if available
                label = desc if desc else resolve_ontology_term(tid)
                
                if label:
                    # Avoid redundant Label (Label) if desc == tid
                    if label.lower() == tid.lower():
                        term_str = label
                    else:
                        term_str = f"{label} ({tid})"
                else:
                    term_str = tid
                
                formatted_terms.append(f"{term_str} (FDR: {format_p(p)})")
            
            term_line = f"- {cat}: " + "; ".join(formatted_terms)
            if total > TOP_TERMS_PER_SOURCE:
                term_line += f" ... *(Total significant: {total}. See string_terms.csv for full list.)*"
            else:
                term_line += f" *(Total significant: {total})*"
            lines.append(term_line)

    if enrichr_terms:
        lines.append("\n**Enrichr annotation enrichment:**")
        lines.append(f"*(Thresholds: Adj. P < {FDR_CUTOFF}, Combined Score > {ENRICHR_COMBINED_SCORE_MIN}. Showing top {TOP_TERMS_PER_SOURCE} terms per library.)*")
        for lib, terms in enrichr_terms.items():
            # terms is list of (term_name, adj_p, combined_score)
            total = enrichr_counts.get(lib, 0) if enrichr_counts else len(terms)
            formatted_terms = [f"{t} (Adj. P: {format_p(p)}, Combined Score: {c:.2f})" for t, p, c in terms]
            
            term_line = f"- {lib}: " + "; ".join(formatted_terms)
            if total > TOP_TERMS_PER_SOURCE:
                term_line += f" ... *(Total significant: {total}. See enrichr_terms.csv for full list.)*"
            else:
                term_line += f" *(Total significant: {total})*"
            lines.append(term_line)

    return "\n".join(lines)

# -----------------------------
# Main
# -----------------------------

def main(data_dir, expression_file):
    if not os.path.exists(data_dir):
        print(f"❌ Error: Data directory not found: {data_dir}")
        return
        
    if not os.path.exists(expression_file):
        print(f"❌ Error: Expression file not found: {expression_file}")
        return

    clusters = load_gene_clusters(data_dir)
    if not clusters:
        print(f"⚠ Warning: No gene cluster files (cluster_k*_cluster_*_genes.txt) found in {data_dir}")
        print("  Skipping biological description generation.")
        return

    sil_scores = load_silhouette_scores(data_dir)
    config = parse_run_config(data_dir)
    if config is None:
        config = {}
    
    run_name = config.get('RUN_NAME', 'GTEx Analysis')
    report_title = f"{run_name} Gene-Tissue Expression and Enrichment Report"
    tag = infer_tag(data_dir)

    string_path = os.path.join(data_dir, "string_terms.csv")
    enrichr_path = os.path.join(data_dir, "enrichr_terms.csv")
    
    string_df = pd.read_csv(string_path) if os.path.exists(string_path) else pd.DataFrame(columns=["k", "cluster", "term", "fdr", "category"])
    enrichr_df = pd.read_csv(enrichr_path) if os.path.exists(enrichr_path) else pd.DataFrame(columns=["k", "cluster", "term", "adjusted_p", "combined_score", "library"])
    
    try:
        expr_df = pd.read_csv(expression_file)
    except Exception as e:
        print(f"❌ Error reading expression file: {e}")
        return

    if "median_tpm" in expr_df.columns:
        expr_df["value"] = expr_df["median_tpm"]
    elif "mean_tpm" in expr_df.columns:
        expr_df["value"] = expr_df["mean_tpm"]
    else:
        print("❌ Error: Expression file must contain 'mean_tpm' or 'median_tpm' columns.")
        return

    report = []
    report.append(f"# {report_title}\n")
    
    html_report_data = defaultdict(dict)

    try:
        for k in sorted(clusters):
            report.append(f"\n## Clustering with k = {k}\n")

            for cluster in sorted(clusters[k]):
                genes = clusters[k][cluster]
                
                # Get silhouette stats
                cluster_sils = [sil_scores.get((k, g), None) for g in genes]
                valid_sils = [s for s in cluster_sils if s is not None]
                avg_sil = sum(valid_sils) / len(valid_sils) if valid_sils else None
                
                # Sort genes by silhouette score if available
                genes_with_scores = []
                for g in genes:
                    s = sil_scores.get((k, g), None)
                    genes_with_scores.append((g, s))
                
                # Sort descending by score (None at end)
                genes_with_scores.sort(key=lambda x: x[1] if x[1] is not None else -1, reverse=True)
                
                # Format gene list string
                gene_strs = []
                for g, s in genes_with_scores[:TOP_GENES]:
                    if s is not None:
                        gene_strs.append(f"{g} (sil: {s:.2f})")
                    else:
                        gene_strs.append(g)

                tissues = summarize_tissues(expr_df, genes)
                string_terms, string_counts = summarize_string_enrichment(string_df, k, cluster)
                enrichr_terms, enrichr_counts = summarize_enrichr_enrichment(enrichr_df, k, cluster)

                report.append(f"\n### Cluster {cluster} (k{k})\n")
                
                # Gene list line
                gene_list_line = f"**Genes ({len(genes)}):** {', '.join(gene_strs)}"
                if len(genes) > TOP_GENES:
                    gene_list_line += f" … *(Showing top {TOP_GENES}. See cluster_k{k}_cluster_{cluster}_genes.txt for full list.)*"
                report.append(gene_list_line + "\n")

                description = speculative_description(
                    genes, tissues, string_terms, enrichr_terms, avg_sil,
                    string_counts=string_counts, enrichr_counts=enrichr_counts
                )
                report.append(description + "\n")
                
                # Collect for HTML
                html_report_data[k][cluster] = {
                    'genes_with_scores': genes_with_scores[:TOP_GENES],
                    'description': description,
                    'n_genes': len(genes)
                }

        out_path = os.path.join(data_dir, "CLUSTER_BIOLOGICAL_DESCRIPTIONS.md")
        with open(out_path, "w") as f:
            f.write("\n".join(report))
        print(f"✓ Markdown report written to {out_path}")
        
        # Generate HTML Report
        generate_html_report(data_dir, html_report_data, config, tag, report_title)

    except Exception as e:
        import traceback
        print(f"❌ Error generating report: {e}")
        traceback.print_exc()



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate biological interpretations for GTEx expression clusters"
    )
    parser.add_argument("--data-dir", required=True, help="Directory with clustering & enrichment outputs")
    parser.add_argument(
        "--expression-file",
        required=True,
        help="Gene–tissue expression CSV (mean_tpm or median_tpm)"
    )

    args = parser.parse_args()
    main(args.data_dir, args.expression_file)
