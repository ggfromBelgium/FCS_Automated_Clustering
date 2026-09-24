# FCS Automated Clustering — FlowSOM & UMAP Pipeline

Automated R pipeline for flow cytometry data analysis: interactive gating,
asinh transformation, FlowSOM clustering with robust cluster-number estimation
(Dynamic Tree Cut, silhouette, gap statistic, bootstrap stability), UMAP
visualisation, AI-assisted cluster annotation (DeepSeek), and comparative
analysis across conditions with statistical testing.

> **Note:** The R source code contains French comments and French console
> messages (the original development language). This documentation is provided
> in English for international accessibility. This does not affect execution.

## Main features

| Step | Description |
|---|---|
| 1–2 | FCS loading, configurable asinh transformation |
| 3–4 | Interactive gating (target populations, FSC.A/FSC.H doublets, FSC.A/SSC.A) |
| 5–6 | Channel exclusion, subsampling (global / balanced by condition or file) |
| 7 | **Batch-effect diagnostic** (UMAP by condition/file) + optional **CytoNorm** normalisation |
| 8–9 | FlowSOM (10×10 SOM) + robust *k* selection (bootstrap Jaccard ≥ 0.75, consensus) |
| 10–11 | Meta-clustering, AI cluster annotation (DeepSeek) with contradiction audit |
| 12 | Visualisations: heatmap, UMAP, node→metacluster plot, interactive filtering/merging |
| 13 | Condition-level analysis: proportions, heatmaps, stacked bars, Welch tests (BH-FDR) |

## Usage

```r
source("FCS_Automated_Clustering.R")
```

The pipeline is **interactive**: it prompts for the FCS file, lets the user
draw gates manually (`locator()`), and guides parameter selection through the
console. Run it in RStudio or an interactive R session (not via batch `Rscript`).

### DeepSeek API key (optional automatic annotation)

The key is read in this order: `~/.deepseek_api_key` → `DEEPSEEK_API_KEY`
environment variable → interactive prompt. **Never commit an API key**. 
Without a key, the pipeline falls back to manual cluster
annotation.

## Dependencies

```r
# Bioconductor
BiocManager::install(c("flowCore", "FlowSOM", "CytoNorm", "PeacoQC"))
# CRAN
install.packages(c("uwot", "dplyr", "ComplexHeatmap", "Cairo", "ggpubr",
                   "cluster", "viridis", "ConsensusClusterPlus", "ggplot2",
                   "reshape2", "ggalluvial", "data.table", "httr", "jsonlite",
                   "scales", "RColorBrewer", "digest", "keyring",
                   "dynamicTreeCut", "fpc", "renv"))
```

Each run produces a timestamped `renv.lock` and a `session_info.txt` in the
output folder, ensuring full version traceability.

## Data

FCS files are not included in this repository (file size + consent). They are
available on [Zenodo: https://doi.org/10.5281/zenodo.22893790]. The pipeline works with any
concatenated FCS file containing `FileID`/`FolderID` columns.

## Outputs

For each sample, a folder `outputs/<sample>/` is created containing: PDFs
(histograms, gates, UMAP, heatmap, alluvial plot, CytoNorm densities), CSVs
(cluster medians, proportions, statistical tests, analysis metrics), and
annotation reports.

## License

Code under the [MIT License](LICENSE). Associated data under [CC-BY 4.0]
(deposited on Zenodo).

## Citation

If you use this code, please cite: [to be completed after publication],
and this GitHub repository (see `CITATION.cff`).
BibTeX entries for the underlying R packages are provided in [`CITATIONS.bib`](CITATIONS.bib).

## Disclaimer

This pipeline is an analysis-aid tool. AI-generated annotations must be
**validated by an experienced cytometrist** before any biological
interpretation.
