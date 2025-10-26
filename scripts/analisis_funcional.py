#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analisis_funcional.py — Functional enrichment analysis of genes
Author: [Your Name]
Date: [Date]

Description:
  This script performs a functional enrichment analysis of genes
  listed in a comma-separated file "genes_input.txt" (do not modify it).
  It fetches functional annotations from MyGene, performs enrichment
  using Enrichr API, and generates bar plots for top terms.
  Results are saved in the "results" folder, overwriting previous outputs.
"""

import argparse
import os
import sys
import requests
import pandas as pd
import mygene
import matplotlib.pyplot as plt
import seaborn as sns

# -------------------- Utilities --------------------
def read_gene_file(path: str) -> list[str]:
    """Read comma-separated genes from a file and return a clean list."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")
    with open(path, "r") as f:
        content = f.read()
    genes = [g.strip().upper() for g in content.split(",") if g.strip()]
    return genes


def ensure_results_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


# -------------------- Functional Annotations --------------------
def fetch_annotations(genes: list[str]) -> pd.DataFrame:
    """Retrieve gene annotations from MyGene."""
    MITO_GENES = ["ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
                  "CO1", "CO2", "CO3", "CYB", "ATP6", "ATP8"]
    corrected_genes = ["MT-" + g if g in MITO_GENES else g for g in genes]

    if corrected_genes != genes:
        print("Detected mitochondrial genes. Adjusted query names:")
        for orig, corr in zip(genes, corrected_genes):
            print(f"  {orig} → {corr}")

    mg = mygene.MyGeneInfo()
    annotations = mg.querymany(
        corrected_genes,
        scopes='symbol',
        fields=['name', 'go.BP.term', 'go.MF.term', 'go.CC.term', 'pathway.kegg.name'],
        species='human'
    )
    df = pd.DataFrame(annotations)
    return df


# -------------------- Enrichment with Enrichr --------------------
def enrichr_analysis(genes: list[str], libraries: list[str]) -> dict:
    """Submit genes to Enrichr API and fetch top enriched terms."""
    add_list_url = "https://maayanlab.cloud/Enrichr/addList"
    enrich_url = "https://maayanlab.cloud/Enrichr/enrich"

    gene_str = "\n".join(genes)
    payload = {'list': (None, gene_str), 'description': (None, 'Functional analysis')}
    response = requests.post(add_list_url, files=payload)
    if response.status_code != 200:
        raise RuntimeError("Error submitting gene list to Enrichr.")

    user_list_id = response.json()['userListId']
    enrichment_results = {}

    for lib in libraries:
        params = {'userListId': user_list_id, 'backgroundType': lib}
        enrich_response = requests.get(enrich_url, params=params)
        if enrich_response.status_code != 200:
            print(f"Error retrieving enrichment results for {lib}")
            continue
        results = enrich_response.json()
        if lib not in results:
            continue
        enrichment_results[lib] = results[lib][:5]  # Keep top 5 terms
    return enrichment_results


# -------------------- Plotting --------------------
def plot_enrichment(enrichment_results: dict, out_dir: str) -> None:
    """Create barplots for enrichment results (overwriting old ones)."""
    sns.set(style="whitegrid")

    for lib, terms_data in enrichment_results.items():
        if not terms_data:
            continue

        terms = [entry[1] for entry in terms_data]
        scores = [entry[4] for entry in terms_data]

        plt.figure(figsize=(8, 5))
        sns.barplot(x=scores, y=terms, palette="viridis")
        plt.xlabel("Combined Score")
        plt.ylabel("Term")
        plt.title(f"Top 5 Enriched Terms: {lib}")
        plt.tight_layout()

        # Fixed filename → overwrite each run
        plot_filename = f"top5_{lib.replace(' ', '_')}.png"
        plot_path = os.path.join(out_dir, plot_filename)
        plt.savefig(plot_path, dpi=150)
        plt.close()

        print(f"Saved barplot for {lib} to '{plot_path}'")


# -------------------- CLI --------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Functional enrichment analysis CLI.")
    parser.add_argument("--input", type=str, default="../data/genes_input.txt",
                        help="Input file with genes (comma-separated).")
    parser.add_argument("--results", type=str, default="../results",
                        help="Folder to save results and plots.")
    return parser.parse_args()


# -------------------- Main --------------------
def main():
    args = parse_args()
    ensure_results_dir(args.results)

    # Step 1: Load genes
    genes = read_gene_file(args.input)
    print("Genes loaded from file:")
    print(genes)
    print("=" * 60)

    # Step 2: Fetch annotations
    print("Fetching gene annotations from MyGene...")
    annotations_df = fetch_annotations(genes)
    annotations_csv_path = os.path.join(args.results, "gene_annotations.csv")
    annotations_df.to_csv(annotations_csv_path, index=False)
    print(f"Functional annotations saved to '{annotations_csv_path}'")
    print("=" * 60)

    # Step 3: Enrichment analysis
    print("Performing enrichment analysis using Enrichr...\n")
    libraries = ['KEGG_2021_Human', 'GO_Biological_Process_2023']
    enrichment_results = enrichr_analysis(genes, libraries)

    print("Top enrichment results:")
    for lib, results in enrichment_results.items():
        print(f"\n--- {lib} ---")
        for entry in results:
            print(f"{entry[1]} (p={entry[2]:.3e}, combined score={entry[4]:.2f})")
    print("=" * 60)

    # Step 4: Plot barplots (overwrites existing files)
    plot_enrichment(enrichment_results, args.results)
    print("Analysis completed successfully.")


if __name__ == "__main__":
    main()
