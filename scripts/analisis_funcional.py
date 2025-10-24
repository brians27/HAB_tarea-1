#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
analisis_funcional.py

Script completo per analisi funzionale dei geni:
- Conversione simboli → Entrez/Ensembl
- Arricchimento funzionale Enrichr
- Interazioni STRING
- Analisi GOATOOLS opzionale
"""

import os
import sys
import requests
import pandas as pd
from typing import List, Optional

# -----------------------------
# MyGene per conversione ID
# -----------------------------
try:
    import mygene
except ImportError:
    print("Installa 'mygene' con pip install mygene")
    sys.exit(1)

# -----------------------------
# gseapy per Enrichr
# -----------------------------
try:
    from gseapy import enrichr
except ImportError:
    print("Installa 'gseapy' con pip install gseapy")
    sys.exit(1)

# -----------------------------
# GOATOOLS opzionale
# -----------------------------
try:
    from goatools.obo_parser import GODag
    from goatools.go_enrichment_ns import GOEnrichmentStudy
    from goatools.associations import read_ncbi_gene2go
    GOATOOLS_AVAILABLE = True
except ImportError:
    GOATOOLS_AVAILABLE = False

# -----------------------------
# Conversione simboli → Entrez/Ensembl
# -----------------------------
class IDConverter:
    def __init__(self):
        self.mg = mygene.MyGeneInfo()

    def convert(self, symbols: List[str]) -> pd.DataFrame:
        res = self.mg.querymany(symbols, scopes='symbol', fields='entrezgene,ensembl.gene,name', species='human')
        rows = []
        for r in res:
            rows.append({
                'query': r.get('query'),
                'symbol': r.get('symbol'),
                'entrez': r.get('entrezgene'),
                'ensembl': r.get('ensembl', {}).get('gene') if isinstance(r.get('ensembl'), dict) else None,
                'name': r.get('name'),
                'found': not r.get('notfound', False)
            })
        return pd.DataFrame(rows)

# -----------------------------
# STRING API
# -----------------------------
class STRINGdb:
    STRING_API = "https://version-11-5.string-db.org/api/json"

    def __init__(self):
        self.session = requests.Session()

    def get_interactions(self, identifier: str, species: int = 9606, required_score: int = 700) -> List[dict]:
        url = f"{self.STRING_API}/network"
        params = {'identifier': identifier, 'species': species, 'required_score': required_score}
        try:
            r = self.session.get(url, params=params, timeout=15)
            if r.status_code != 200:
                return []
            return r.json()
        except Exception:
            return []

# -----------------------------
# Enrichr
# -----------------------------
def run_enrichr(genes: List[str], outdir: str) -> pd.DataFrame:
    libs = [
        'GO_Biological_Process_2023',
        'GO_Cellular_Component_2023',
        'GO_Molecular_Function_2023',
        'KEGG_2021_Human',
        'Reactome_2022'
    ]
    all_res = []
    for lib in libs:
        try:
            enr = enrichr(
                gene_list=genes,
                gene_sets=lib,
                organism='Human',
                description='enrichr_analysis',
                outdir=outdir,
                cutoff=0.05,
                no_plot=True
            )
            df = enr.results
            if df is not None and not df.empty:
                df['library'] = lib
                all_res.append(df)
        except Exception as e:
            print(f'Enrichr fallito per {lib}: {e}')

    if all_res:
        df_all = pd.concat(all_res, ignore_index=True)
    else:
        # crea file vuoto se non ci sono risultati
        df_all = pd.DataFrame(columns=[
            'Term','Overlap','P-value','Adjusted P-value','Old P-value','Old Adjusted',
            'Z-score','Combined Score','Genes','library'
        ])
    df_all.to_csv(os.path.join(outdir, 'enrichr_combined.tsv'), sep='\t', index=False)
    print(f"Arricchimento Enrichr salvato in: {os.path.join(outdir, 'enrichr_combined.tsv')}")
    return df_all

# -----------------------------
# GOATOOLS opzionale
# -----------------------------
def run_goatools(study_genes: List[str], gene2go_path: str, obo_path: str, outdir: str) -> Optional[pd.DataFrame]:
    if not GOATOOLS_AVAILABLE:
        print('GOATOOLS non disponibile. Saltando analisi GO locale.')
        return None

    godag = GODag(obo_path)
    associations = read_ncbi_gene2go(gene2go_path, taxids=[9606])
    background = set(associations.keys())

    if any(not g.isdigit() for g in study_genes):
        print('GOATOOLS richiede Entrez Gene IDs per study_genes.')
        return None

    study_set = set(int(g) for g in study_genes)

    goeaobj = GOEnrichmentStudy(
        list(background),
        associations,
        godag,
        propagate_counts=True,
        alpha=0.05,
        methods=['fdr_bh']
    )

    goea_results = goeaobj.run_study(list(study_set))

    rows = []
    for r in goea_results:
        rows.append({
            'GO': r.GO,
            'name': r.name,
            'p_uncorrected': r.p_uncorrected,
            'p_fdr_bh': r.p_fdr_bh,
            'study_count': r.study_count,
            'study_items': ';'.join(map(str, r.study_items))
        })

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(outdir, 'goatools_enrichment.tsv'), sep='\t', index=False)
    print(f"Analisi GOATOOLS salvata in: {os.path.join(outdir, 'goatools_enrichment.tsv')}")
    return df

# -----------------------------
# Main
# -----------------------------
def main():
    import argparse

    parser = argparse.ArgumentParser(description="Analisi funzionale dei geni")
    parser.add_argument("--input", required=True, help="File di geni (uno per linea)")
    parser.add_argument("--out", required=True, help="Cartella di output")
    parser.add_argument("--gene2go", required=False, help="File gene2go (per GOATOOLS)")
    parser.add_argument("--obo", required=False, help="File go-basic.obo (per GOATOOLS)")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # Legge geni
    with open(args.input) as f:
        genes = [line.strip() for line in f if line.strip()]

    # Conversione simboli
    print("Converting symbols to Entrez/Ensembl...")
    converter = IDConverter()
    df_ids = converter.convert(genes)
    df_ids.to_csv(os.path.join(args.out, 'id_mapping.tsv'), sep='\t', index=False)
    print(f"ID mapping salvato in: {os.path.join(args.out, 'id_mapping.tsv')}")

    # STRING interazioni
    stringdb = STRINGdb()
    all_interactions = []
    print("Consultando interazioni STRING...")
    for gene in genes:
        interactions = stringdb.get_interactions(gene)
        for inter in interactions:
            inter['gene'] = gene
        all_interactions.extend(interactions)
    # Salva interazioni in file
    if all_interactions:
        pd.DataFrame(all_interactions).to_csv(os.path.join(args.out, 'string_interactions.tsv'), sep='\t', index=False)
        print(f"Interazioni STRING salvate in: {os.path.join(args.out, 'string_interactions.tsv')}")

    # Enrichr
    print("Eseguendo arricchimento Enrichr...")
    run_enrichr(genes, args.out)

    # GOATOOLS opzionale
    if args.gene2go and args.obo:
        print("Eseguendo analisi GOATOOLS locale...")
        entrez_genes = [str(g) for g in df_ids['entrez'] if g is not None]
        run_goatools(entrez_genes, args.gene2go, args.obo, args.out)

    print("Analisi completa! Risultati in:", args.out)

if __name__ == "__main__":
    main()
