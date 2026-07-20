# Multi-omics Cancer Analysis Project

Differential gene expression (DGE) analysis of TCGA RNA-seq data, comparing tumor vs. normal
samples within two cancer types (**BRCA** – breast cancer, **LUAD** – lung adenocarcinoma) and
directly comparing the two cancer types against each other.

## Overview

This project analyzes bulk RNA-seq gene counts from The Cancer Genome Atlas (TCGA) using
[`pydeseq2`](https://pydeseq2.readthedocs.io/), the Python implementation of DESeq2, together
with `pandas`, `scikit-learn`, `seaborn`, and `gseapy`.

The analysis pipeline covers:

1. **Data loading & preprocessing** – matching RNA-seq counts to sample metadata, converting
   counts to integers, filtering low-expression genes.
2. **Differential expression analysis (per cancer type)** – Tumor vs. Normal for BRCA and LUAD,
   run separately with `pydeseq2` (`DeseqDataSet` + `DeseqStats`).
3. **Result export** – DE result tables (`log2FoldChange`, `padj`, etc.) as CSV.
4. **Volcano plots** – visualizing significant genes per cancer type.
5. **Principal Component Analysis (PCA)** – on variance-stabilized/normalized expression across
   all samples, including a scree plot and 95% confidence ellipses per group.
6. **Heatmap** – expression of the top differentially expressed genes across all samples,
   annotated by project and tissue type.
7. **BRCA vs. LUAD comparison**
   - Overlap of significant genes found in the two separate Tumor-vs-Normal analyses.
   - A direct DESeq2 comparison of BRCA tumors vs. LUAD tumors.
8. **Gene enrichment analysis** – over-representation analysis (ORA) of significant genes against
   GO Biological Process and KEGG gene sets using `gseapy`/Enrichr.

## Data

| File | Description |
|---|---|
| `rna_seq_counts.csv` | Raw RNA-seq gene counts. Rows = samples, columns = genes. |
| `rna_seq_metadata.csv` | Sample metadata. Key columns: `File Name` (sample ID, matches the row index of the counts file), `Project ID` (`TCGA-BRCA` / `TCGA-LUAD`), `condition`, `Tissue Type` (`Tumor` / `Normal`). |

Data originates from [TCGA](https://www.cancer.gov/ccg/research/genome-sequencing/tcga) and is
not redistributed in this repository unless explicitly included — see the notebook for expected
file locations.

> **Note:** TCGA counts are provided as (rounded) floats rather than raw integers; the notebook
> rounds them to integers before running DESeq2, as required by `pydeseq2`.

## Requirements

```bash
pip install pandas numpy matplotlib seaborn scikit-learn scipy pydeseq2 gseapy
```

`gseapy`'s enrichment step queries the [Enrichr](https://maayanlab.cloud/Enrichr/) web service and
therefore requires an internet connection at runtime; all other steps run fully offline.

## Usage

1. Place `rna_seq_counts.csv` and `rna_seq_metadata.csv` in the project directory (or update the
   file paths at the top of the notebook).
2. Open and run `TCGA_analysis.ipynb` top to bottom in Jupyter.
3. Result tables and plots are generated inline; DE result tables are additionally exported as
   CSV:
   - `deseq2_results_BRCA.csv`
   - `deseq2_results_LUAD.csv`
   - `deseq2_results_BRCA_vs_LUAD.csv`
   - `enrichment_results_BRCA.csv`
   - `enrichment_results_LUAD.csv`

## Methods summary

- **Design formula:** `~Tissue_Type` (per-project Tumor vs. Normal), `~Project_ID` (direct
  BRCA-vs-LUAD comparison on tumor samples only).
- **Significance thresholds:** `padj < 0.05` and `|log2FoldChange| > 1`, used consistently for
  volcano plots, heatmap gene selection, and enrichment input.
- **Low-count filtering:** genes with a total count below 10 across the relevant sample subset
  are removed before fitting.
- **PCA:** top 500 most variable genes (log2 of size-factor normalized counts).
- **Gene Enrichment Analysis**

## Repository structure

```
.
├── rna_seq_counts.csv          # RNA-seq gene counts (samples x genes)
├── rna_seq_metadata.csv        # Sample metadata
├── TCGA_analysis.ipynb         # Main analysis notebook
└── README.md
```
