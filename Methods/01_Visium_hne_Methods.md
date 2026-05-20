# 01 — Visium H&E Analysis Methods

**Tutorial source:** https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html  
**Platform:** 10x Visium  
**Dataset:** Mouse Brain coronal section (H&E stained, pre-annotated)  
**Spots:** 2,688 | **Genes:** ~32,000

---

# 1. Data Loading

The dataset was loaded via Squidpy's built-in dataset module, which downloads a pre-processed AnnData and ImageContainer:

```python
import squidpy as sq

img   = sq.datasets.visium_hne_image()    # ~314 MB on first download
adata = sq.datasets.visium_hne_adata()
```

## Visium Platform Background

10x Visium captures spatially resolved gene expression by placing a tissue section on a slide printed with a grid of ~5,000 barcoded spots (~55 µm diameter, ~100 µm center-to-center spacing). mRNA released from the tissue binds to the barcoded oligos on each spot, is reverse-transcribed, and sequenced. The resulting data is a spot-by-gene count matrix where each spot's physical coordinate on the slide is known.

The Visium slide also captures a high-resolution brightfield tissue image (H&E stained), which is aligned to the spot grid. This image is the basis for the image feature analysis in this tutorial.

## AnnData Structure

| Slot | Contents |
|---|---|
| `adata.X` | Normalised count matrix (2,688 spots × genes) |
| `adata.obs['cluster']` | Pre-annotated brain region labels (14 clusters) |
| `adata.obsm['spatial']` | `(x, y)` pixel coordinates of each spot on the tissue image |
| `adata.uns['spatial']` | Image metadata, scale factors, spot diameter |

## ImageContainer

The `ImageContainer` is Squidpy's data structure for storing and manipulating the tissue image alongside the AnnData. It supports multiscale images and lazy loading.

---

# 2. Spatial Visualisation of Cluster Annotations

The pre-annotated brain region clusters were visualised directly on the tissue image:

```python
sq.pl.spatial_scatter(adata, color='cluster')
```

`sq.pl.spatial_scatter()` overlays spot-level annotations on the tissue image using the `adata.obsm['spatial']` coordinates and the scale factors stored in `adata.uns['spatial']`.

---

# 3. Image Feature Extraction

## Rationale

A key capability of Squidpy is extracting quantitative features from the tissue image for each spot. This allows tissue morphology information — which is captured in the H&E image but not in the gene expression matrix — to be incorporated into downstream analysis. The hypothesis tested here is that image features alone can partially recapitulate the biological cluster structure identified by gene expression.

## Feature Types

Summary features compute statistical summaries (mean, standard deviation, quantiles) of pixel intensity values within each spot's image region:

```python
sq.im.calculate_image_features(
    adata,
    img.compute(),
    features='summary',
    key_added='features_summary_scale1.0',
    n_jobs=4,
    scale=1.0,
)
```

The `scale` parameter controls the spatial resolution at which features are extracted — a scale of `1.0` uses the full resolution image, while larger scales use downsampled versions that capture coarser tissue context.

## Multi-scale Extraction

Features were extracted at two scales (`1.0` and `2.0`) and concatenated:

```python
adata.obsm['features'] = pd.concat(
    [adata.obsm[f] for f in adata.obsm.keys() if 'features_summary' in f],
    axis='columns',
)
```

Multi-scale features capture both fine-grained cellular texture (scale `1.0`) and broader tissue architecture (scale `2.0`).

---

# 4. Clustering on Image Features

A temporary AnnData object was created from the image feature matrix and processed with the standard dimensionality reduction and clustering pipeline:

```python
def cluster_features(features):
    tmp = ad.AnnData(features)
    sc.pp.scale(tmp)
    sc.pp.pca(tmp, n_comps=10)
    sc.pp.neighbors(tmp)
    sc.tl.leiden(tmp, key_added='leiden_features')
    return tmp.obs['leiden_features']

adata.obs['image_cluster'] = cluster_features(adata.obsm['features'])
```

The resulting image-based clusters were compared to the gene expression-based cluster annotations to assess how much biological signal is captured by tissue morphology alone.

---

# 5. Spatial Graph Construction

A grid-based spatial neighbour graph was constructed, connecting each spot to its six nearest neighbours on the Visium hexagonal grid:

```python
sq.gr.spatial_neighbors(adata, coord_type='grid', n_neighs=6)
```

`coord_type='grid'` is appropriate for Visium data where spots are arranged on a regular hexagonal grid. This differs from the Delaunay-based graph used for Xenium, where cells are irregularly distributed.

---

# 6. Neighborhood Enrichment

See `01_xenium_methods.md — Section 6` for the general method description.

For Visium H&E, the cluster key is `'cluster'` (pre-annotated brain regions) rather than Leiden clusters:

```python
sq.gr.nhood_enrichment(adata, cluster_key='cluster')
```

The enrichment matrix reveals which brain regions are anatomically adjacent — for example, cortical layers are expected to neighbor each other, while distant regions (e.g., cortex and cerebellum) should show low or negative enrichment.

---

# 7. Co-occurrence Score

```python
sq.gr.co_occurrence(adata, cluster_key='cluster')
```

The co-occurrence score was computed for the Hippocampus cluster across a range of spatial distances (`0–800 µm`). The distance-dependent profile shows which regions are spatially proximal to the hippocampus at different scales.

---

# 8. Interaction Matrix

```python
sq.gr.interaction_matrix(adata, cluster_key='cluster')
```

---

# 9. Centrality Scores

```python
sq.gr.centrality_scores(adata, cluster_key='cluster')
```

In the context of the mouse brain, high degree centrality indicates regions that border many other regions (e.g., white matter tracts that run throughout the section). High betweenness centrality indicates regions that act as spatial bridges between otherwise separated areas.

---

# 10. Ripley’s Statistics

Ripley’s statistics are a classical spatial statistics tool for testing whether a point pattern (here: spots of a given cluster type) is more clustered, more dispersed, or randomly distributed relative to a homogeneous Poisson process.

The L statistic is a normalised transformation of Ripley’s K function that makes the expected value under complete spatial randomness (CSR) a straight line at `L(r) = r`:

```python
sq.gr.ripley(adata, cluster_key='cluster', mode='L')
```

Values of `L(r)` above the diagonal indicate spatial clustering at distance `r`; values below indicate dispersion.

---

# 11. Spatial Autocorrelation — Moran’s I

```python
sq.gr.spatial_autocorr(adata, mode='moran', n_perms=100, n_jobs=4)
```

`n_perms=100` permutations were used to compute the significance of each gene’s Moran’s I score. The top spatially variable genes were visualised on the tissue image.

---

# Key Parameters Summary

| Step | Parameter | Value | Rationale |
|---|---|---|---|
| Image features | `scale` | `1.0, 2.0` | Multi-scale captures fine and coarse morphology |
| Image features | `features` | `'summary'` | Statistical pixel intensity summaries |
| Image feature PCA | `n_comps` | `10` | Sufficient for ~20 image features |
| Spatial graph | `coord_type` | `'grid'` | Visium hexagonal grid layout |
| Spatial graph | `n_neighs` | `6` | Standard for hexagonal Visium grid |
| Moran’s I | `n_perms` | `100` | Permutation-based significance testing |

---

# References

1. Palla, G. et al. (2022). *Squidpy: a scalable framework for spatial omics analysis.* Nature Methods, 19, 171–178.

2. Ripley, B.D. (1976). *The second-order analysis of stationary point processes.* Journal of Applied Probability, 13, 255–266.

3. Moran, P.A.P. (1950). *Notes on continuous stochastic phenomena.* Biometrika, 37, 17–23.

4. 10x Genomics Visium documentation: https://www.10xgenomics.com/spatial-transcriptomics
