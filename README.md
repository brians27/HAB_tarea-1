
# Task 1: Functional Analysis

This project implements a **functional analysis of the genes COX4I2, ND1, and ATP6** using Python.  
The analysis retrieves functional annotations through **MyGene.info** and performs an **Over-Representation Analysis (ORA)** using the **Enrichr API**.  
The results include both a table of functional annotations and bar plots of the most enriched biological processes and pathways.

---

## Repository Structure

```
/functional-analysis/
├── data/
│ └── genes_input.txt                       # Provided genes for the analysis (do not modify)
├── scripts/
│ └── analisis_funcional.py                 # Main functional analysis script
├── results/                                # Automatically generated results
│ ├── gene_annotations.csv                  # Functional annotations (GO + KEGG)
│ ├── top5_KEGG_2021_Human.png              # Bar plot for KEGG enrichment
│ └── top5_GO_Biological_Process_2023.png   # Bar plot for GO Biological Process enrichment
├── README.md                               # Project description and instructions
└── requirements.txt                        # Python dependencies
```

## Command-Line Usage (CLI)

Run the script from the terminal as follows:

```bash
cd scripts/
python analisis_funcional.py --input ../data/genes_input.txt --results ../results
```

## Dependencies

Install all necessary dependencies present in `requirements.txt`:

```bash
pip install -r requirements.txt
```
