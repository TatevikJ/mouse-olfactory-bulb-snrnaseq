# mouse-olfactory-bulb-snrnaseq

Single-nucleus RNA-seq of the mouse olfactory bulb in the ASO model of Parkinson's disease
(mice overexpressing human alpha-synuclein) and wild-type controls.

- 12 mice: 6 ASO and 6 WT
- Parse Biosciences split-pool libraries, 302,758 nuclei before filtering
- Reference: GRCm39 (Ensembl v113) with the human SNCA transgene added (`hSNCA-transgene`)
- Cell types: Allen Brain Cell Atlas whole mouse brain taxonomy (CCN20230722), assigned with MapMyCells

This repository holds the analysis scripts only. Count matrices and results are not included.

## Pipeline

All scripts are in `mouse_ob_GRCm39_v113_hSNCA/src/` and are run from `mouse_ob_GRCm39_v113_hSNCA/`
(for example `Rscript src/normalize_pre.R`). Paths inside the scripts are relative to that folder.
Outputs are written to `outputs/`.

| Step | Script | What it does |
|---|---|---|
| 1 | `merge_count_matrices.py` | Merge the 12 per-sample Parse count matrices into one h5ad with raw counts; keep genes detected in at least 3 nuclei |
| 2 | MapMyCells (web tool, no script) | Annotate the merged h5ad with the Allen taxonomy (hierarchical mapping) |
| 3 | `convert_h5ad_to_rds.R` | Build a Seurat object from the h5ad and attach the MapMyCells labels |
| 4 | `normalize_pre.R` | Normalize, find variable genes, scale, PCA, UMAP |
| 5 | `qc_pre.R` | QC metrics and plots before filtering (mitochondrial fraction, scDblFinder doublet score, mapping confidence) |
| 6 | `qc_filtering_v1.R`, `qc_filtering_v2.R` | Basic QC filter: nCount_RNA < 400,000, percent.mt < 0.5, doublet score < 0.2 (v1) or < 0.3 (v2) |
| 7 | `qc_post_v1.R`, `qc_post_v2.R` | QC plots and tables after the basic filter |
| 8 | `bootstrap_diagnostics_post_v1_v2.R` | Diagnostic: mapping confidence against nFeature_RNA and nCount_RNA per subclass |
| 9 | `cellcycle_diagnostics_post_v1_v2.R` | Diagnostic: cell-cycle scores of nuclei with very high gene counts |
| 10 | `qc_filtering_v3.R`, `qc_filtering_v4.R` | Remove nuclei mapped to supertypes of brain regions outside the olfactory bulb (v3 = v1 + region filter, v4 = v2 + region filter) |
| 11 | `qc_post_v3.R`, `qc_post_v4.R` | QC plots and tables after the region filter |
| 12 | `normalize_post_v3.R`, `normalize_post_v4.R` | Renormalize the filtered object, Leiden clustering, FindAllMarkers, and write a protein-coding-only object |
| 13 | `qc_post_v3_normalized.R`, `qc_post_v4_normalized.R` | UMAPs, transgene expression and composition tables on the final objects |
| 14 | `markers_0149_vs_0148.R` | Markers separating the two OB Trdn Gaba supertypes (0149 vs 0148) |
| 15 | `deg_post_v3_normalized.R`, `deg_post_v4_normalized.R` | ASO vs WT differential expression per class, subclass and supertype |
| 16 | `gsea_post_v3_normalized.R`, `gsea_post_v4_normalized.R` | Gene set enrichment analysis (GO biological process) on the ranked DEG tables |
| 17 | `ora_post_v3_normalized.R`, `ora_post_v4_normalized.R` | Over-representation analysis (GO biological process) of up and down DEGs |

Steps 8 and 9 are diagnostics and are not needed for the later steps.

## QC versions

| Version | Doublet score cutoff | Brain-region filter | Nuclei kept |
|---|---|---|---|
| v1 | < 0.2 | no | 221,974 |
| v2 | < 0.3 | no | 226,833 |
| v3 | < 0.2 | yes | 220,078 |
| v4 | < 0.3 | yes | 224,395 |

Differential expression and pathway analysis are run on v3 and v4.

## Differential expression

`deg_post_v*_normalized.R` compares ASO and WT within each cell type on protein-coding genes, with three methods:

- **Pseudobulk (edgeR):** counts summed per mouse, quasi-likelihood test.
  DEG: FDR < 0.05 and |log2FC| > 0.5.
- **Single-nucleus (MAST):** hurdle model on nuclei, at most 2,000 per condition.
  DEG: Bonferroni < 0.01 and |coef| > 0.5.
- **Single-nucleus (Wilcoxon):** Seurat FindMarkers on nuclei, all genes reported.
  DEG: no cutoff applied.

## Software

R 4.3 with Seurat 5, scDblFinder, scran, scuttle, edgeR, MAST, clusterProfiler, enrichplot, org.Mm.eg.db, EnhancedVolcano
and anndata. Python 3.11 with anndata, pandas, numpy and scipy for the merge step.
