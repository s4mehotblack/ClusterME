# GTEx Gene Expression Analysis Pipeline

A comprehensive suite of tools for extracting, clustering, and analyzing gene expression patterns using GTEx data. This pipeline takes a list of genes and produces biological insights through tissue-specificity clustering, statistical validation, and functional enrichment analysis.

## 🚀 Quick Start

The easiest way to run the pipeline is via the interactive wizard, which handles file paths, creates directories, and manages the workflow for you.

```bash
./run_pipeline_wizard.sh
```

Alternatively, you can run individual steps manually as described in the [Pipeline Steps](#-pipeline-steps) section.

---

## 📋 Prerequisites

### Python Environment
- **Python 3.7+**
- Install required packages using the provided requirements file:
  ```bash
  pip install -r requirements.txt
  ```

### Data Sources
1.  **Gene List**: You must provide your own list of genes of interest (one gene symbol per line in a `.txt` file).
2.  **GTEx Data**: Download the following from the [GTEx Portal (Adult Data)](https://gtexportal.org/home/downloads/adult-gtex). Select **"Bulk Tissue Expression"** and **"Metadata"** data types to obtain:
    *   **Gene TPMs**: `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz` (or v10 equivalent).
    *   **Sample Attributes**: `GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt` (or v10 equivalent).

---

## 🛠️ Pipeline Steps

### 1. Data Extraction (`extract_all_chunked.py`)
Extracts RNA-seq expression data (TPM) for your genes from the massive GTEx dataset. 
*   **Optimization**: Extracted matrices are saved to `analysis_runs/` for reuse across different analysis projects.

### 2. Clustering (`gtex_cluster_v2.py`)
Groups genes based on their tissue expression profiles using Ward's hierarchical clustering.
*   Calculates **Silhouette Scores** for every gene to quantify cluster membership quality.

### 3. Visualization (`visualise_silhouette_scores.py`)
Generates plots to assess cluster cohesion and identify "core" vs. ambiguous cluster members.

### 4. Enrichment Analysis
*   **STRING (`string_cluster_enrichment_report.py`)**: Queries protein interaction networks and functional pathways.
*   **Enrichr (`enrichr_cluster_enrichment_report.py`)**: Checks for over-representation in GO (2025), KEGG (2025), Human Phenotype Ontology, and more. Logs all four statistical scores (p-value, adjusted p-value, Z-score, and Combined Score).

### 5. Reporting (`generate_biological_descriptions.py`)
Synthesizes all results into:
*   **Markdown Report**: `CLUSTER_BIOLOGICAL_DESCRIPTIONS.md`
*   **HTML Dashboard**: `CLUSTER_BIOLOGICAL_REPORT.html` (includes hyperlinks to all PNG/CSV/TXT artifacts).

---

## 🔧 Support & Diagnostic Scripts

Beyond the core pipeline, the following scripts assist in parameter selection and quality control:

### Cluster Number Suggestion (`suggest_k.py`)
Helps determine the optimal number of clusters ($k$) by testing a range of values and calculating Silhouette and SSE (elbow) metrics.
```bash
python suggest_k.py --expression-file your_data.csv --kmin 2 --kmax 15 --z-score
```

### Cluster Noise Visualization (`visualise_cluster_noise.py`)
Generates line plots of gene expression profiles across tissues to qualitatively assess how "noisy" or "tight" a cluster is.
```bash
# Plot all clusters in a results directory
./visualise_cluster_noise.py --expression-file your_data.csv --data-dir results_k7 --clusters all
```

---

## 📂 Output Files (in `--data-dir`)

| File Type | Filename Pattern | Description |
|-----------|------------------|-------------|
| **Clusters** | `gene_clusters_k{k}_*.csv` | Table assigning each gene to a cluster. |
| **Scores** | `gene_silhouette_scores_k{k}_*.csv` | Silhouette score for every gene. |
| **Heatmaps** | `global_heatmap_k{k}_*.png` | Visualizations of expression patterns. |
| **Quality Plots**| `silhouette_analysis_k{k}_*.png` | Cluster cohesion visualization. |
| **Reports** | `CLUSTER_BIOLOGICAL_REPORT.html` | **Interactive HTML results dashboard.** |

---

## 💡 Best Practices
1.  **Use Z-scores**: Always use the `--z-score` flag to identify relative tissue specificity rather than just expression magnitude.
2.  **Run `suggest_k.py`**: Before committing to a final analysis, use the suggestion tool to find the most natural clustering resolution for your gene set.
3.  **Inspect the HTML Report**: Use the hyperlinks in the HTML report to quickly toggle between heatmaps and gene lists for deeper investigation.
