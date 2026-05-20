# 03 — Visium Fluorescence Analysis Methods

**Tutorial source:** https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_fluo.html  
**Platform:** 10x Visium  
**Dataset:** Mouse Brain coronal crop (fluorescence stained)  
**Key focus:** Image segmentation and segmentation-based feature extraction

---

## 1. Data Loading

```python
img   = sq.datasets.visium_fluo_image_crop()    # ~370 MB on first download
adata = sq.datasets.visium_fluo_adata_crop()
```

This dataset is a cropped region of a mouse brain coronal section stained with fluorescence markers rather than H&E. Fluorescence imaging differs from H&E in a fundamental way: instead of colour-coded morphological staining, fluorescence channels encode specific molecular targets with quantitative intensity values. This makes fluorescence images more amenable to automated cell segmentation, as foreground (cells) and background can be separated by intensity thresholding.

---

## 2. Image Processing Pipeline

### 2a. Gaussian smoothing

Before segmentation, the image was smoothed with a Gaussian filter to reduce noise and improve segmentation quality:

```python
sq.im.process(
    img,
    layer='image',
    method='smooth',
)
```

Gaussian smoothing averages each pixel with its neighbours weighted by a Gaussian kernel, suppressing high-frequency noise while preserving the overall structure of cells. The smoothed image is stored as a new layer (`'image_smooth'`) in the ImageContainer.

### 2b. Watershed segmentation

Cell segmentation was performed using the watershed algorithm:

```python
sq.im.segment(
    img=img,
    layer='image_smooth',
    method='watershed',
    channel=0,
)
```

**How watershed works:**
1. The image is treated as a topographic surface where pixel intensity represents elevation
2. Local intensity minima are identified as seeds (cell centers)
3. The algorithm "floods" from each seed outward, filling catchment basins
4. Watershed lines — boundaries between adjacent basins — become cell borders

Watershed segmentation is well-suited to fluorescence microscopy images because cell nuclei typically appear as bright, well-separated peaks in a nuclear stain channel, providing reliable seed points.

The segmentation result is stored as `img['segmented_watershed']` — a label image where each cell is assigned a unique integer label, and background pixels have label 0.

---

## 3. Segmentation Feature Extraction

Features derived from the segmentation mask were extracted for each Visium spot:

```python
sq.im.calculate_image_features(
    adata,
    img,
    layer='image',
    key_added='features_segmentation',
    n_jobs=4,
    features='segmentation',
    features_kwargs={
        'segmentation': {
            'label_layer': 'segmented_watershed',
            'props': ['label', 'area', 'mean_intensity'],
            'channels': [0, 1],
        }
    },
)
```

### Segmentation feature types

| Feature | Description | Biological interpretation |
|---|---|---|
| `label` | Number of unique cell labels per spot | Cell count per spot (cell density) |
| `area` | Mean area of segmented cells per spot | Average cell size |
| `mean_intensity` | Mean fluorescence intensity per cell per channel | Channel-specific protein abundance |

These features are computed per-spot by collecting all segmented cells whose centroids fall within each spot's image region and computing summary statistics over them.

**Cell count per spot** is particularly informative: dense cellular regions (e.g., cortical layers, hippocampus) will have higher cell counts per spot than cell-sparse regions (e.g., white matter). This signal is complementary to gene expression.

---

## 4. Multi-scale Summary Feature Extraction

In addition to segmentation features, multi-scale summary image features were extracted following the same approach as Tutorial 02:

```python
for scale in [1.0, 2.0]:
    sq.im.calculate_image_features(
        adata,
        img.compute(),
        features='summary',
        key_added=f'features_summary_scale{scale}',
        n_jobs=4,
        scale=scale,
    )

adata.obsm['features'] = pd.concat(
    [adata.obsm[f] for f in adata.obsm.keys() if 'features_summary' in f],
    axis='columns',
)
```

---

## 5. Image Feature Clustering

Spots were clustered based on the combined image feature matrix (segmentation + summary features):

```python
def cluster_features(features):
    tmp = ad.AnnData(features)
    sc.pp.scale(tmp)
    sc.pp.pca(tmp, n_comps=min(10, features.shape[1] - 1))
    sc.pp.neighbors(tmp)
    sc.tl.leiden(tmp, key_added='leiden_img')
    return tmp.obs['leiden_img']

adata.obs['image_cluster'] = cluster_features(adata.obsm['features'])
```

### Comparison to gene expression clusters

The image-based clusters were compared to the pre-annotated gene expression clusters. The degree of agreement between the two cluster assignments reflects how much biological information is encoded in the fluorescence image morphology — specifically, whether the image features capture the same cell type variation that gene expression does.

Perfect agreement would mean that tissue morphology fully predicts cell identity. Partial agreement (the typical result) means that image features capture some but not all of the variation distinguishing cell types — likely the aspects that produce visible morphological differences (e.g., cell size, density, nuclear shape) but not those that require molecular profiling (e.g., transcription factor expression).

---

## 6. Key Differences from H&E Analysis (Tutorial 02)

| Aspect | Visium H&E (Tutorial 02) | Visium Fluorescence (Tutorial 03) |
|---|---|---|
| Image type | Brightfield H&E | Fluorescence (2+ channels) |
| Image features | Summary statistics only | Summary + segmentation features |
| Cell segmentation | Not performed | Watershed segmentation |
| Per-cell features | No | Yes (area, intensity, count per spot) |
| Quantitative image data | No (staining is qualitative) | Yes (fluorescence intensity is quantitative) |
| Spatial stats | Full suite | Not performed (focus on image analysis) |

---

## Key Parameters Summary

| Step | Parameter | Value | Rationale |
|---|---|---|---|
| Smoothing | `method` | `'smooth'` (Gaussian) | Reduce noise before segmentation |
| Segmentation | `method` | `'watershed'` | Standard for fluorescence cell segmentation |
| Segmentation | `channel` | 0 | Nuclear stain channel (first channel) |
| Seg features | `props` | `['label','area','mean_intensity']` | Cell count, size, and intensity per spot |
| Summary features | `scale` | 1.0, 2.0 | Multi-scale morphological context |
| Feature PCA | `n_comps` | min(10, n_features-1) | Dimensionality reduction for clustering |

---

## References

- Palla, G. et al. (2022). Squidpy: a scalable framework for spatial omics analysis. *Nature Methods*, 19, 171–178.
- Meyer, F. & Beucher, S. (1990). Morphological segmentation. *Journal of Visual Communication and Image Representation*, 1(1), 21–46. *(Watershed algorithm)*
- van der Walt, S. et al. (2014). scikit-image: image processing in Python. *PeerJ*, 2, e453. *(Segmentation implementation)*
