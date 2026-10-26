# Visium — Tutorial

**Read the tutorial:** https://nightwalkie.github.io/visium-guided-clustering-tutorial/

An English, Seurat-style guided workflow from unnormalized Visium counts to QC,
normalization, variable features, PCA, graph-based clustering, UMAP, and tissue maps.

## Contents

- **Visium v1 and v2/CytAssist:** worked examples with figures generated from public
  10x Genomics mouse brain datasets.
- **Visium HD bins:** a complete sketch-clustering and full-data projection template.
- **Visium HD segmented cells:** a separate cell-level clustering template.

The HD templates are explicitly distinguished from executed results. Example QC
and clustering parameters are teaching choices, not universal quality thresholds.

## Files

- `index.html` — the complete self-contained tutorial, including its fonts and figures.
- `Visium_standard_workflow.R` — standard-resolution count-to-UMAP code.
- `Visium_HD_bins_workflow.R` — HD bin import, sketch clustering, and projection.
- `Visium_HD_cells_workflow.R` — segmentation-derived cell analysis.
- `.nojekyll` — preserves the existing static HTML presentation on GitHub Pages.

The HTML file is published directly, without conversion to Markdown or a replacement
theme. The three R files remain beside it so their download links work.

## Attribution and reproducibility

The article includes references, dataset provenance, run parameters, and session
information. Public example data are from **10x Genomics (CC BY 4.0)**. Figures were
generated for these notes. Styling follows the Seurat PBMC3k tutorial's open-source
Flatly/pkgdown conventions; third-party style and font notices are embedded in HTML.

These are independent learning notes, not an official Seurat or 10x Genomics release.
No private sequencing data or workstation credentials are included.
