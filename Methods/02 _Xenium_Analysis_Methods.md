# 02 — Xenium Analysis Methods

**Tutorial source:** https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_xenium.html  
**Platform:** 10x Xenium  
**Dataset:** Human Lung Cancer (FFPE)  
**Cells:** ~161,000 | **Genes:** 480 (targeted panel)

---

## 1. Data Loading and Format

### Platform background

10x Xenium is a subcellular-resolution in-situ sequencing platform. Unlike spot-based Visium, Xenium assigns individual fluorescently detected transcripts to cells using a combination of cell segmentation masks and spatial proximity. Each transcript is recorded with its (x, y, z) coordinates within the tissue section.

Xenium output consists of:
- Cell-by-gene count matrix
- Per-transcript coordinate tables (~40M+ rows)
- Cell and nucleus segmentation masks (as polygons and label images)
- Morphology focus images (multiscale TIFF)

### SpatialData format

Data was loaded using `spatialdata-io`'s `xenium()` reader, which parses the raw Xenium output directory into a `SpatialData` object — a unified container for multi-modal spatial data developed by the scverse consortium (Marconato et al., 2024).

```python
from spatialdata_io import xenium
import spatialdata as sd

sdata = xenium("data/Xenium/")
sdata.write("data/Xenium.zarr")
sdata = sd.read_zarr("data/Xenium.zarr")
```

The `SpatialData` object contains:

| Slot | Contents |
|---|---|
| `sdata.images` | Morphology focus image (multiscale) |
| `sdata.labels` | Cell and nucleus segmentation masks |
| `sdata.points` | Per-transcript coordinate table |
| `sdata.shapes` | Cell/nucleus boundary polygons |
| `sdata.tables["table"]` | AnnData: 161,000 cells × 480 genes |

### Zarr storage

The SpatialData object was written to Zarr format before analysis. Zarr is a chunked, compressed array format designed for large datasets that do not fit in memory. It allows lazy loading — only the data chunks needed for a given operation are loaded, rather than the entire dataset. This is essential for datasets at the scale of Xenium.

---

## 2. Quality Control

### Metrics used

Three QC metrics were calculated and visualised:

| Metric | Source | Interpretation |
|---|---|---|
| `transcript_counts` | Pre-computed in Xenium output | Total transcripts detected per cell |
| `cell_area` | Pre-computed in Xenium output | Area of the cell segmentation polygon (µm²) |
| `control_probe_counts` | Xenium negative controls | Non-specific/background signal |

### Filtering

Cells with fewer than 10 transcript counts were removed:

```python
sc.pp.filter_cells(adata, min_counts=10)
```

**Note:** The QC thresholds for Xenium differ from standard scRNA-seq. Because the gene panel is targeted (480 genes), the number of detected genes per cell is bounded by the panel size and is not a reliable indicator of cell quality in the same way as whole-transcriptome data. Transcript count and cell area are more informative.

Cells with very high transcript counts were not filtered, as high-expressing cells (e.g., actively secreting cells) are biologically valid in this context.

---

## 3. Preprocessing

Standard scRNA-seq preprocessing was applied to the count matrix:

```python
sc.pp.normalize_total(adata, target_sum=1e4)   # total-count normalise
sc.pp.log1p(adata)                              # log(x+1) transform
sc.pp.highly_variable_genes(adata,              # select HVGs
    flavor='seurat_v3', n_top_genes=3000)
sc.tl.pca(adata)                                # PCA
sc.pp.neighbors(adata)                          # k-NN graph
sc.tl.umap(adata)                               # UMAP embedding
sc.tl.leiden(adata, key_added='leiden')         # Leiden clustering
```

**HVG selection:** Because the panel is only 480 genes, `n_top_genes=3000` effectively retains most genes. The `seurat_v3` flavour was used as it accounts for the mean-variance relationship in count data.

**Leiden vs Louvain:** Leiden clustering was used as it produces better-connected communities than Louvain and is the current recommended algorithm in the scverse ecosystem.

---

## 4. Spatial Visualisation

Clusters were visualised both in UMAP space and in physical tissue coordinates using `sq.pl.spatial_scatter()`:

```python
sq.pl.spatial_scatter(
    adata,
    color='leiden',
    shape=None,   # plot cells as points, not polygons
    size=2,
)
```

`shape=None` plots each cell as a dot at its centroid coordinate rather than drawing the full cell boundary polygon — the latter is more accurate but computationally expensive at 161,000 cells.

---

## 5. Spatial Graph Construction

A spatial neighbour graph was constructed using Delaunay triangulation, which connects each cell to its natural geometric neighbours without requiring a fixed radius:

```python
sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True)
```

The resulting graph is stored in `adata.obsp['spatial_connectivities']` and `adata.obsp['spatial_distances']`. All subsequent spatial statistics use this graph.

---

## 6. Neighborhood Enrichment

Neighborhood enrichment tests whether pairs of cell types are spatially co-localized more than expected by chance. For each pair of clusters (i, j), a z-score is computed by comparing the observed number of i–j neighboring pairs against a null distribution generated by randomly permuting cluster labels while preserving the graph structure.

```python
sq.gr.nhood_enrichment(adata, cluster_key='leiden')
```

A positive z-score indicates spatial co-enrichment (the two cell types are neighbors more often than chance). A negative z-score indicates spatial avoidance.

---

## 7. Co-occurrence Score

The co-occurrence score measures how the probability of observing two cell types within a given spatial distance changes as a function of that distance, relative to a random baseline:

```python
sq.gr.co_occurrence(adata, cluster_key='leiden')
```

A score > 1 at distance d indicates that the two cell types are more likely to co-occur within distance d than expected by chance. This metric captures distance-dependent spatial relationships that the neighborhood enrichment (which uses a fixed graph) cannot.

---

## 8. Spatial Autocorrelation — Moran's I

Moran's I quantifies the degree of spatial autocorrelation for a continuous variable (e.g., gene expression). A value of I ≈ 1 indicates strong positive spatial autocorrelation (nearby cells have similar expression), I ≈ 0 indicates random spatial distribution, and I ≈ -1 indicates spatial dispersion.

```python
sq.gr.spatial_autocorr(adata, mode='moran')
```

The output (`adata.uns['moranI']`) ranks all genes by their Moran's I score, identifying **spatially variable genes (SVGs)** — genes whose expression is non-randomly distributed in tissue space.

---

## 9. Centrality Scores

Graph-theoretic centrality metrics were computed for each cluster to characterise its structural role in the spatial tissue graph:

| Metric | Definition |
|---|---|
| Degree centrality | Proportion of cells in the graph that are neighbors of cluster i |
| Closeness centrality | How close cluster i is on average to all other clusters in the graph |
| Betweenness centrality | How often cluster i lies on the shortest path between other clusters |

```python
sq.gr.centrality_scores(adata, cluster_key='leiden')
```

---

## 10. Interaction Matrix

The interaction matrix counts the number of observed edges between each pair of clusters in the spatial neighbor graph. Unlike neighborhood enrichment (which normalises for expected frequency), the interaction matrix reports raw interaction counts:

```python
sq.gr.interaction_matrix(adata, cluster_key='leiden')
```

---

## Key Parameters Summary

| Step | Parameter | Value | Rationale |
|---|---|---|---|
| QC filter | `min_counts` | 10 | Remove debris/empty cells |
| Normalisation | `target_sum` | 10,000 | Standard total-count normalisation |
| HVG selection | `n_top_genes` | 3,000 | Retains most of the 480-gene panel |
| HVG flavour | `seurat_v3` | — | Accounts for mean-variance relationship |
| Spatial graph | `delaunay` | True | Geometry-based, no fixed radius needed |
| Clustering | Leiden | default resolution | Community detection on k-NN graph |

---

## References

- Marconato, L. et al. (2024). SpatialData: an open and universal data framework for spatial omics. *Nature Methods*, 21, 1102–1110.
- Moran, P.A.P. (1950). Notes on continuous stochastic phenomena. *Biometrika*, 37, 17–23.
- Traag, V.A. et al. (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports*, 9, 5233.
- 10x Genomics Xenium documentation: https://www.10xgenomics.com/platforms/xenium
