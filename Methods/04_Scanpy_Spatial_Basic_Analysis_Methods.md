# 04 — Scanpy Spatial Basic Analysis Methods

**Tutorial source:** https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html  
**Author:** Giovanni Palla  
**Datasets:** 10x Visium Mouse Brain + MERFISH Mouse Cortex  
**Key focus:** Foundational spatial scRNA-seq workflow in Scanpy

---

## Overview

This tutorial establishes the foundational workflow for spatial transcriptomics analysis using Scanpy. It covers two datasets and two platforms — Visium (spot-based, ~55 µm resolution) and MERFISH (cell-based, subcellular resolution) — demonstrating how the same analytical framework applies across different spatial technologies.

The key addition over standard scRNA-seq analysis is the use of spatial coordinates to overlay expression data and cluster annotations directly on the tissue image, enabling spatially-informed interpretation.

---

## Part A — Visium Mouse Brain

### 1. Data Loading

```python
adata = sq.datasets.visium_hne_adata()
```

The dataset was loaded from Squidpy's dataset module. See [02_visium_hne_methods.md](02_visium_hne_methods.md) for background on the Visium platform and data structure.

The critical spatial information is stored in:
- `adata.obsm['spatial']` — (x, y) pixel coordinates of each spot
- `adata.uns['spatial']` — image data, scale factors, and spot diameter metadata

### 2. Quality Control

```python
adata.var['mt'] = adata.var_names.str.startswith('Mt-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
```

Mitochondrial genes were identified using the `Mt-` prefix (mouse gene nomenclature uses capitalised first letter followed by lowercase, unlike human `MT-`). QC metrics were visualised as histograms:

| Metric | Threshold applied | Rationale |
|---|---|---|
| `total_counts` | Visual inspection | Spots with very low counts are likely off-tissue |
| `n_genes_by_counts` | `sc.pp.filter_genes(min_cells=10)` | Remove genes detected in too few spots |
| `pct_counts_mt` | Visual inspection | High MT% may indicate tissue damage |

For this pre-processed tutorial dataset, no additional cell filtering was required beyond gene filtering.

### 3. Preprocessing

```python
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=2000)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added='leiden_clusters', resolution=0.5)
```

**Resolution parameter:** A resolution of 0.5 produces fewer, broader clusters than the default (1.0), which is appropriate for brain region-level annotation where over-clustering would fragment anatomically coherent regions.

**HVG flavour:** The `'seurat'` flavour (original Seurat method) was used here rather than `'seurat_v3'`. It selects genes based on mean expression and dispersion (coefficient of variation), which is a slightly simpler approach than the v3 method but performs similarly on well-characterised datasets.

### 4. Spatial Visualisation

The key function for spatial transcriptomics visualisation in Scanpy is `sc.pl.spatial()` (or equivalently `sq.pl.spatial_scatter()`):

```python
sq.pl.spatial_scatter(adata, color='cluster')
```

This function:
1. Reads the tissue image from `adata.uns['spatial']`
2. Plots the image as background
3. Overlays each spot as a circle at its spatial coordinate, coloured by the specified variable
4. Applies the appropriate scale factors to map pixel coordinates to the rendered image

**What spatial visualisation reveals that UMAP cannot:** Clusters that appear as a single UMAP blob may occupy completely different anatomical regions in the tissue. Conversely, cells that are spatially adjacent may belong to transcriptionally distinct clusters (e.g., the sharp boundary between cortex and white matter). The spatial view is essential for biological interpretation.

### 5. Marker Gene Identification

```python
sc.tl.rank_genes_groups(adata, 'cluster', method='wilcoxon')
```

Marker genes were identified using the Wilcoxon rank-sum test comparing each cluster against all others. The top marker genes were then visualised spatially:

```python
sq.pl.spatial_scatter(adata, color=marker_genes)
```

Plotting marker gene expression on the tissue reveals the spatial extent and boundaries of each transcriptionally defined region — for example, hippocampal markers should show expression concentrated in the hippocampal formation.

---

## Part B — MERFISH Mouse Cortex

### 1. Platform background

MERFISH (Multiplexed Error-Robust Fluorescence In Situ Hybridization) is a highly multiplexed single-molecule fluorescence in situ hybridisation (smFISH) method developed by Zhuang et al. (2016). It assigns transcript identities to individual fluorescent spots using combinatorial barcoding with error-robust encoding, allowing simultaneous detection of hundreds to thousands of genes at single-cell and subcellular resolution.

Key differences from Visium:

| Property | Visium | MERFISH |
|---|---|---|
| Unit of observation | Spot (~55 µm, multiple cells) | Individual cell |
| Resolution | ~55 µm | Subcellular |
| Gene coverage | Whole transcriptome | Targeted panel (~500 genes) |
| Throughput | ~5,000 spots per slide | ~50,000–1,000,000 cells |
| Data type | Count matrix per spot | Count matrix per cell with (x,y) coordinates |

### 2. Data Loading

```python
adata_merfish = sq.datasets.merfish()
```

The MERFISH mouse cortex dataset contains ~3,000 cells with pre-annotated cell type labels (`adata.obs['Cell_class']`).

### 3. Spatial Visualisation

```python
sq.pl.spatial_scatter(
    adata_merfish,
    color='Cell_class',
    shape=None,
    size=2,
)
```

`shape=None` plots cells as points (no polygon boundaries). `size=2` uses small points appropriate for the high cell density of the cortex. The resulting spatial map shows the cortical layer organisation — excitatory neuron subtypes form distinct horizontal bands corresponding to cortical layers L2/3, L4, L5, and L6, while inhibitory neurons and non-neuronal cells are distributed throughout.

### 4. Spatial Statistics — Neighborhood Enrichment

```python
sq.gr.spatial_neighbors(adata_merfish, coord_type='generic', delaunay=True)
sq.gr.nhood_enrichment(adata_merfish, cluster_key='Cell_class')
```

Delaunay triangulation was used as the spatial graph method, appropriate for irregularly distributed single cells. The neighborhood enrichment matrix reveals which cell type pairs are spatially adjacent — in the cortex, adjacent layer-specific excitatory neuron subtypes (e.g., L4 and L5) are expected to show positive enrichment, while distant layer types (e.g., L2/3 and L6) show lower enrichment.

---

## Key Parameters Summary

| Section | Step | Parameter | Value | Rationale |
|---|---|---|---|---|
| Visium | Gene filter | `min_cells` | 10 | Remove rarely detected genes |
| Visium | Clustering | `resolution` | 0.5 | Brain region-level granularity |
| Visium | HVG | `n_top_genes` | 2,000 | Standard for Visium whole-transcriptome data |
| Visium | HVG flavour | `seurat` | — | Simple dispersion-based HVG selection |
| MERFISH | Spatial graph | `delaunay` | True | Appropriate for irregular single-cell layout |
| MERFISH | Spatial graph | `coord_type` | `'generic'` | Non-grid cell layout |

---

## Comparison of All Four Tutorials

| Tutorial | Platform | Unit | Gene coverage | Spatial resolution | Image available |
|---|---|---|---|---|---|
| 01 Xenium | Xenium | Cell | Targeted (480) | Subcellular | Yes (multiscale) |
| 02 Visium H&E | Visium | Spot | Whole transcriptome | ~55 µm | Yes (H&E) |
| 03 Visium Fluo | Visium | Spot | Whole transcriptome | ~55 µm | Yes (fluorescence) |
| 04 Scanpy Spatial | Visium + MERFISH | Spot / Cell | Both | ~55 µm / subcellular | Yes (H&E) |

---

## References

- Wolf, F.A. et al. (2018). Scanpy: large-scale single-cell gene expression data analysis. *Genome Biology*, 19, 15.
- Palla, G. et al. (2022). Squidpy: a scalable framework for spatial omics analysis. *Nature Methods*, 19, 171–178.
- Chen, K.H. et al. (2015). Spatially resolved, highly multiplexed RNA profiling in single cells. *Science*, 348, aaa6090. *(MERFISH)*
- Zhang, M. et al. (2021). Spatially resolved cell atlas of the mouse primary motor cortex by MERFISH. *Nature*, 598, 137–143.
- Traag, V.A. et al. (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports*, 9, 5233.
