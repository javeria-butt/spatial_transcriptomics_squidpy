[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1r4MjFgaYZAxcxDeTNY5BO3gUACop6ui8?usp=sharing)
# Spatial Transcriptomics Analysis with Squidpy & Scanpy

**Author:** Javeria Butt   
**Tools:** Squidpy | Scanpy | SpatialData | AnnData  
**Reference genome / platforms:** | 10x Xenium | 10x Visium | MERFISH

---

## Table of Contents

- [Overview](#overview)
- [Tutorials Covered](#tutorials-covered)
- [Datasets](#datasets)
- [Installation](#installation)
- [Usage](#usage)
- [Results Summary](#results-summary)
- [Detailed Methods](#detailed-methods)
- [Discussion](#discussion)
- [References](#references)

---

## Overview

This repository documents the completion of four spatial transcriptomics tutorials using **Squidpy** and **Scanpy**. Spatial transcriptomics allows gene expression to be measured while preserving the physical location of each cell or spot within the tissue 鈥?adding a crucial spatial dimension to single-cell analysis.

The four tutorials cover three different spatial transcriptomics platforms (Xenium, Visium H&E, Visium Fluorescence) and one foundational Scanpy-based spatial workflow, together demonstrating the complete spectrum of spatial data analysis from raw data loading to spatial statistics and neighborhood enrichment.

| # | Tutorial | Platform | Tissue | Key analysis |
|---|---|---|---|---|
| 1 | Xenium | 10x Xenium | Human Lung Cancer | QC, spatial stats, neighborhood enrichment |
| 2 | Visium H&E | 10x Visium | Mouse Brain (coronal) | Image features, spatial graph, co-occurrence |
| 3 | Visium Fluorescence | 10x Visium | Mouse Brain (coronal crop) | Image segmentation, segmentation features |
| 4 | Scanpy Spatial | 10x Visium + MERFISH | Mouse Brain / Mouse Cortex | Basic spatial QC, clustering, marker genes |

---

### 1. Analyze Visium H&E Data
**Source:** [squidpy.readthedocs.io 鈥?tutorial_visium_hne](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html)  


Analyzes a **10x Visium** H&E stained coronal section of the mouse brain. Visium captures gene expression at spots (~55 碌m diameter) arrayed on a slide, with a high-resolution tissue image captured alongside.

Key steps:
- Load pre-processed Visium AnnData and ImageContainer from Squidpy datasets
- Visualize cluster annotations in spatial context
- Extract multi-scale summary image features using `sq.im.calculate_image_features()`
- Cluster spots based on image features alone (morphology-based clustering)
- Combine gene expression and image feature embeddings
- Spatial graph analysis: neighborhood enrichment, co-occurrence, interaction matrix, centrality scores, Ripley's statistics, spatial autocorrelation (Moran's I)

---


### 2. Analyze Xenium Data
**Source:** [squidpy.readthedocs.io 鈥?tutorial_xenium](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_xenium.html)  
**Notebook:** [`Notebooks/01_xenium_analysis.ipynb`](Notebooks/01_xenium_analysis.ipynb)

Analyzes a **10x Xenium** dataset of Human Lung Cancer (161,000 cells 脳 480 genes). Xenium is a subcellular-resolution spatial transcriptomics platform that assigns individual transcripts to cells using in-situ sequencing.

Key steps:
- Load Xenium data using `spatialdata-io` and convert to Zarr format
- Calculate QC metrics (transcript counts, control probes, cell area)
- Standard scRNA-seq preprocessing: filter, normalize, log1p, HVG, PCA, UMAP, Leiden clustering
- Visualize clusters in UMAP and spatial coordinates
- Compute spatial statistics: spatial autocorrelation (Moran's I), co-occurrence scores, neighborhood enrichment, interaction matrix, centrality scores
- Interactive visualization with `napari-spatialdata`

---


### 2. Analyze Visium Fluorescence Data
**Source:** [squidpy.readthedocs.io 鈥?tutorial_visium_fluo](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_fluo.html)  


Analyzes a **10x Visium** fluorescence-stained coronal section (crop) of the mouse brain. Demonstrates Squidpy's image processing capabilities for fluorescence data.

Key steps:
- Load pre-processed Visium fluorescence AnnData and ImageContainer
- Image segmentation using `sq.im.segment()` (thresholding + watershed)
- Extract segmentation features: cell count, mean intensity, area per spot
- Extract and combine multi-scale summary features
- Cluster spots based on image features using Leiden clustering
- Correlation between image-based and gene expression-based cluster assignments

---

### 4. Scanpy Spatial Basic Analysis
**Source:** [scanpy-tutorials.readthedocs.io 鈥?spatial/basic-analysis](https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html)  
**Notebook:** [`Notebooks/04_scanpy_spatial_basic.ipynb`](Notebooks/04_scanpy_spatial_basic.ipynb)

Foundational spatial transcriptomics workflow in Scanpy covering Visium and MERFISH data.

Key steps:
- Read Visium data with `sc.read_visium()`
- QC and preprocessing (filter, normalize, log1p, HVG, PCA, UMAP, Leiden)
- Visualization in spatial coordinates using `sc.pl.spatial()`
- Cluster marker gene identification and spatial expression mapping
- MERFISH example: load, preprocess, and visualize mouse cortex data

---

## Datasets

| Tutorial | Dataset | Size | Source | Auto-download? |
|---|---|---|---|---|
| 02 Xenium | Human Lung Cancer (Xenium) | ~10 GB | [10x Genomics](https://www.10xgenomics.com/datasets/preview-data-ffpe-human-lung-cancer-with-xenium-multimodal-cell-segmentation-1-standard) | Manual |
| 01 Visium H&E | Mouse Brain coronal (H&E) | ~314 MB | Squidpy datasets | Auto |
| 03 Visium Fluo | Mouse Brain coronal crop (fluo) | ~370 MB | Squidpy datasets | Auto |
| 04 Scanpy Spatial | Mouse Brain Visium + MERFISH cortex | ~100 MB | Squidpy datasets | Auto |

### Tutorial 02 Manual Xenium download

```bash
# Download from 10x Genomics (requires registration)
# https://www.10xgenomics.com/datasets/preview-data-ffpe-human-lung-cancer-with-xenium-multimodal-cell-segmentation-1-standard

# After downloading, extract to:
mkdir -p data/Xenium
# Place extracted files in Data/Xenium/
```

### Tutorials 01, 03, 04 Auto download

These datasets download automatically when you run the Notebooks:

```python
# Tutorial 01 loads ~314 MB on first run
adata = sq.datasets.visium_hne_adata()
img   = sq.datasets.visium_hne_image()

# Tutorial 03 loads ~370 MB on first run
adata = sq.datasets.visium_fluo_adata_crop()
img   = sq.datasets.visium_fluo_image_crop()

# Tutorial 04 loads Visium + MERFISH datasets
adata = sq.datasets.visium_hne_adata()
```

---

## Installation

### Option A Conda (recommended)

```bash
conda env create -f Environment/environment.yml
conda activate spatial-sq
```

### Option B pip

```bash
pip install -r Environment/requirements.txt
```

### Verify installation

```python
import squidpy as sq
import scanpy as sc
sc.logging.print_header()
print(f"squidpy=={sq.__version__}")
```

---

## Usage

```bash
# 1. Clone the repo
git clone https://github.com/NazarMDin/spatial-transcriptomics-squidpy.git
cd spatial-transcriptomics-squidpy

# 2. Create environment
conda env create -f Environment/environment.yml
conda activate spatial-sq

# 3. Launch notebooks
jupyter notebook Notebooks/
```

Run Notebooks in order (01 to 04). Each notebook is self-contained.

---

## Results Summary

### Tutorial 01 Visium H&E (Mouse Brain)

| Metric | Value |
|---|---|
| Spots | 2,688 |
| Genes | ~32,000 |
| Pre-annotated clusters | 14 brain regions |
| Key finding | Image features alone recapitulate gene expression-based clusters |

### Tutorial 02 Xenium (Human Lung Cancer)

| Metric | Value |
|---|---|
| Cells | 161,000 |
| Genes | 480 |
| Clusters (Leiden) | ~10 |
| Key finding | Distinct spatial organization of tumor vs immune cell populations |

### Tutorial 03 Visium Fluorescence (Mouse Brain crop)

| Metric | Value |
|---|---|
| Spots | Cropped region |
| Key finding | Segmentation features capture cell density differences across brain regions |

### Tutorial 04 Scanpy Spatial (Visium + MERFISH)

| Metric | Value |
|---|---|
| Visium spots | 2,688 (mouse brain) |
| MERFISH cells | Mouse cortex |
| Key finding | Spatial gene expression patterns match known anatomical boundaries |

---

## Detailed Methods

Each tutorial has a dedicated methods file documenting the rationale behind every analytical step, parameter choices, and relevant background on the biology and algorithms used.

| File | Tutorial | Key topics covered |
|---|---|---|
| [`Methods/01_Visium_hne_Methods.md`](Methods/01_Visium_hne_Methods.md) | Visium H&E Mouse Brain | Visium platform, ImageContainer, multi-scale image features, image-based clustering, hexagonal spatial graph, Ripley's L, Moran's I |
| [`Methods/02_Xenium_Analysis_Methods.md`](Methods/02_Xenium_Analysis_Methods.md) | Xenium Human Lung Cancer | SpatialData/Zarr format, Xenium QC, Leiden clustering, Delaunay graph, Moran's I, co-occurrence, neighborhood enrichment, centrality scores |
| [`Methods/03_Visium_Flourescence_Analysis_Methods.md`](Methods/03_Visium_Flourescence_Analysis_Methods.md) | Visium Fluorescence Mouse Brain | Fluorescence vs H&E imaging, Gaussian smoothing, watershed segmentation, segmentation features (cell count, area, intensity) |
| [`Methods/04_Scanpy_Spatial_Basic_Analysis_Methods.md`](Methods/04_Scanpy_Spatial_Basic_Analysis_Mathods.md) | Scanpy Spatial (Visium + MERFISH) | Foundational spatial workflow, spatial visualisation, marker gene mapping, MERFISH platform background, cross-tutorial platform comparison |

---

## Discussion

See [`Discussion/DISCUSSION.md`](Discussion/DISCUSSION.md) for a full interpretation of all four tutorials, comparison across platforms, and key biological findings.

---

## References

- Palla, G. et al. (2022). Squidpy: a scalable framework for spatial omics analysis. *Nature Methods*, 19, 171鈥?78.
- Wolf, F.A. et al. (2018). Scanpy: large-scale single-cell gene expression data analysis. *Genome Biology*, 19, 15.
- Virshup, I. et al. (2021). anndata: Annotated data. *bioRxiv*. https://doi.org/10.1101/2021.12.16.473007
- Marconato, L. et al. (2024). SpatialData: an open and universal data framework for spatial omics. *Nature Methods*, 21, 1102鈥?110.
- Moran, P.A.P. (1950). Notes on continuous stochastic phenomena. *Biometrika*, 37, 17鈥?3. *(Moran's I)*
- 10x Genomics Xenium: https://www.10xgenomics.com/platforms/xenium
- 10x Genomics Visium: https://www.10xgenomics.com/spatial-transcriptomics
