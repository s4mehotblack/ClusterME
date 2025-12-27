#!/usr/bin/env python3

import argparse
import pandas as pd
import plotly.graph_objects as go
from pathlib import Path


def load_clusters(path, gene_col, cluster_col):
    df = pd.read_csv(path)

    if gene_col not in df.columns or cluster_col not in df.columns:
        raise ValueError(
            f"{path} must contain columns: {gene_col}, {cluster_col}"
        )

    return df[[gene_col, cluster_col]].rename(
        columns={gene_col: "gene", cluster_col: "cluster"}
    )


def map_ancestry(low_df, high_df):
    merged = pd.merge(
        low_df.rename(columns={"cluster": "low_cluster"}),
        high_df.rename(columns={"cluster": "high_cluster"}),
        on="gene",
        how="inner"
    )

    overlap = (
        merged
        .groupby(["high_cluster", "low_cluster"])
        .size()
        .reset_index(name="n_shared_genes")
    )

    high_sizes = high_df.groupby("cluster").size().rename("high_cluster_size")
    low_sizes = low_df.groupby("cluster").size().rename("low_cluster_size")

    overlap = overlap.merge(
        high_sizes, left_on="high_cluster", right_index=True
    )
    overlap = overlap.merge(
        low_sizes, left_on="low_cluster", right_index=True
    )

    overlap["fraction_of_high_cluster"] = (
        overlap["n_shared_genes"] / overlap["high_cluster_size"]
    )
    overlap["fraction_of_low_cluster"] = (
        overlap["n_shared_genes"] / overlap["low_cluster_size"]
    )

    overlap = overlap.sort_values(
        ["high_cluster", "fraction_of_high_cluster"],
        ascending=[True, False]
    )

    return overlap


def extract_primary_origins(overlap):
    return (
        overlap
        .sort_values("fraction_of_high_cluster", ascending=False)
        .groupby("high_cluster")
        .first()
        .reset_index()
    )


def make_sankey(overlap, low_label, high_label, output_html):
    overlap = overlap.copy()

    overlap["low_label"] = overlap["low_cluster"].astype(str).apply(
        lambda x: f"{low_label} {x}"
    )
    overlap["high_label"] = overlap["high_cluster"].astype(str).apply(
        lambda x: f"{high_label} {x}"
    )

    labels = pd.Index(
        overlap["low_label"].tolist() + overlap["high_label"].tolist()
    ).unique()

    label_to_index = {label: i for i, label in enumerate(labels)}

    sources = overlap["low_label"].map(label_to_index)
    targets = overlap["high_label"].map(label_to_index)
    values = overlap["n_shared_genes"]

    fig = go.Figure(
        go.Sankey(
            node=dict(
                pad=15,
                thickness=15,
                line=dict(width=0.5),
                label=labels.tolist(),
            ),
            link=dict(
                source=sources,
                target=targets,
                value=values,
            ),
        )
    )

    fig.update_layout(
        title=f"Cluster Ancestry: {low_label} → {high_label}",
        font_size=10
    )

    fig.write_html(output_html)


def main():
    parser = argparse.ArgumentParser(
        description="Map ancestry between two clustering resolutions"
    )

    parser.add_argument("--low", required=True, help="Lower-resolution cluster CSV")
    parser.add_argument("--high", required=True, help="Higher-resolution cluster CSV")

    parser.add_argument("--gene-col", default="gene", help="Gene column name")
    parser.add_argument("--cluster-col", default="cluster", help="Cluster column name")

    parser.add_argument("--low-label", default="k_low", help="Label for low k")
    parser.add_argument("--high-label", default="k_high", help="Label for high k")

    parser.add_argument("--out-prefix", default="cluster_mapping")

    args = parser.parse_args()

    low_df = load_clusters(args.low, args.gene_col, args.cluster_col)
    high_df = load_clusters(args.high, args.gene_col, args.cluster_col)

    overlap = map_ancestry(low_df, high_df)
    primary = extract_primary_origins(overlap)

    overlap_out = f"{args.out_prefix}_full_mapping.csv"
    primary_out = f"{args.out_prefix}_primary_origins.csv"
    sankey_out = f"{args.out_prefix}_sankey.html"

    overlap.to_csv(overlap_out, index=False)
    primary.to_csv(primary_out, index=False)

    make_sankey(
        overlap,
        args.low_label,
        args.high_label,
        sankey_out
    )

    print("Outputs written:")
    print(f"  {overlap_out}")
    print(f"  {primary_out}")
    print(f"  {sankey_out}")


if __name__ == "__main__":
    main()
