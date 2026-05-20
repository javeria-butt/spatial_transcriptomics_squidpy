# Spatial Transcriptomics Tutorials: Platform & Workflow Comparison

## Overview

This repository covers four spatial transcriptomics tutorials spanning three platforms (**Xenium**, **Visium H&E**, **Visium Fluorescence**) and one foundational **Scanpy** spatial workflow. Together, they demonstrate how spatial context fundamentally changes the interpretation of single-cell gene expression data:
* **Context Shifts Interpretation:** Clusters that appear similar in UMAP space can occupy completely different anatomical regions.
* **Biological Proximity:** Cells that appear spatially adjacent can show strong co-localization signals that reflect known, functional biology.

---

## Detailed Tutorial Breakdowns

### Tutorial 01 — Xenium (Human Lung Cancer)
The Xenium platform represents the highest-resolution spatial transcriptomics technology covered in this repository, capturing individual transcript locations within cells at subcellular resolution. 

* **Scale and Storage:** The dataset (161,000 cells × 480 genes) is orders of magnitude larger than typical Visium datasets. Utilizing the `SpatialData`/`Zarr` storage format is essential for handling data at this scale efficiently.
* **Targeted Quality Control (QC):** The QC step for Xenium differs from standard scRNA-seq because the gene panel is targeted (480 curated genes) rather than transcriptome-wide. This means the number of genes per cell is strictly bounded by the panel size. Consequently, quality is better assessed through transcript counts and cell area metrics rather than traditional gene count thresholds.
* **Spatial Architecture & Microenvironment:** Spatial statistics (neighborhood enrichment, co-occurrence) successfully revealed the tissue architecture of the lung cancer sample. Tumor cells clustered tightly together and exhibited distinct co-localization patterns with immune cell populations in the tumor microenvironment. 
* **Autocorrelation:** Moran's I identified genes whose expression is spatially non-random, highlighting potential spatial communication signals between adjacent cell populations.

### Tutorial 02 — Visium H&E (Mouse Brain)
This stands as the most analytically rich of the four tutorials, blending histological morphology with standard transcriptomics.

* **Morphology Recapitulates Expression:** A key finding was that image features extracted from the H&E tissue image alone partially recapitulated gene expression-based cluster assignments. Because different brain regions possess distinct cellular morphologies visible via H&E staining, the image captures a significant portion of the same biological signal as the gene expression matrix.
* **Spatial Graph Analysis:** Confirmed known neuroanatomical relationships:
    * **Neighborhood Enrichment:** Demonstrated that Hippocampus spots preferentially neighbor other Hippocampus spots rather than distant cortical regions, consistent with the compact, well-defined anatomy of the hippocampus.
    * **Ripley's L Statistic:** Confirmed spatial clustering of specific cell populations beyond what would be expected from a completely random distribution.
    * **Moran's I:** Identified the most spatially variable genes (those with the highest spatial autocorrelation), which represent strong candidates for spatially regulated expression programs.
* **Analytical Synergy:** The combination of gene expression and image features is a powerful analytical direction. Morphological features are computationally cheap to extract from existing tissue images, while gene expression is more expensive. Combining both significantly boosts the biological signal available for downstream analysis.

### Tutorial 03 — Visium Fluorescence (Mouse Brain)
This tutorial demonstrated Squidpy's image processing pipeline for fluorescence microscopy data. Fluorescence differs fundamentally from H&E as it carries quantitative intensity information rather than qualitative morphological staining.

* **Segmentation Pipeline:** A robust workflow (**Gaussian smoothing → thresholding → watershed**) successfully separated individual cells within each Visium spot, enabling cell-level features to be extracted and aggregated per spot.
* **Density Mapping:** The cell count per spot varied substantially across brain regions. Dense cellular regions (cortex, hippocampus) showed higher cell counts per spot than cell-sparse regions (white matter tracts), a signal perfectly captured by the segmentation features.
* **Cluster Alignment:** Comparing image-feature-based clustering and gene expression-based clustering showed partial but imperfect agreement. Fluorescence intensity features capture different information than H&E morphology, and the degree of overlap with gene expression clusters depends directly on how much the targeted fluorescence signal correlates with underlying cell type identity.

### Tutorial 04 — Scanpy Spatial (Visium + MERFISH)
This tutorial established the foundational, underlying workflow for spatial transcriptomics in the Scanpy ecosystem.

* **Visualizing Spatial Overlays:** The core conceptual addition over standard scRNA-seq is the use of `sc.pl.spatial()` and `sq.pl.spatial_scatter()` to overlay expression values and cluster annotations directly onto the tissue image. This immediately brings to light spatial patterns that are entirely invisible in standard UMAP space.
* **MERFISH Mouse Cortex Case Study:** Illustrated how subcellular-resolution spatial data differs from spot-based Visium:
    * Individual cells, rather than multi-cellular spots, serve as the baseline unit of observation.
    * Cell type boundaries appear significantly sharper in spatial maps.
    * Spatial statistics operate at higher resolution and can detect much finer-scale organization.
* **Cortical Stratification:** The cortical layer organization was clearly visible in the MERFISH data—excitatory neuron subtypes occupied distinct, well-demarcated cortical layers, inhibitory neurons were interspersed throughout, and non-neuronal cells (astrocytes, oligodendrocytes) followed their own distinct spatial distributions.

---

## Cross-Platform & Tutorial Comparison

| Aspect | Tutorial 01: Xenium | Tutorial 02: Visium H&E | Tutorial 03: Visium Fluo | Tutorial 04: Scanpy Spatial |
| :--- | :--- | :--- | :--- | :--- |
| **Resolution** | Subcellular | ~55 µm spots | ~55 µm spots | Spot / Subcellular |
| **Gene Panel** | Targeted (480 genes) | Whole Transcriptome | Whole Transcriptome | Whole Transcriptome |
| **Image Data** | Morphology image | H&E image | Fluorescence image | H&E image |
| **Scale** | 161,000 cells | 2,688 spots | Cropped region | 2,688 spots |
| **Key Strength** | Cell resolution, large scale | Image + expression integration | Segmentation-based features | Foundational baseline workflow |

---

## Key Biological Findings

1.  **Spatial Organization is Non-Random:** Across all four datasets, spatial statistics confirmed that cell types and clusters are distributed deterministically—organized precisely according to anatomical structures or functional microenvironment relationships.
2.  **Tissue Images Carry Biological Signal:** In both Visium tutorials, image features computed from the tissue slice partially recapitulated gene expression-based cluster assignments, demonstrating that morphological information is highly complementary to transcriptomic data.
3.  **Co-Localization Reflects In Vivo Biology:** Neighborhood enrichment consistently identified cell-type pairs that are known to interact biologically, such as tumor-immune co-localization in lung cancer and region-specific neuron layers in the mouse brain.
4.  **Platform Resolution Dictates Analytical Resolution:** Xenium and MERFISH data (subcellular resolution) revealed cell-level spatial organization patterns that are naturally blurred or averaged out in spot-based Visium data. Both levels of resolution are informative and complementary for tissue indexing.

---

## Technical & Experimental Limitations

* **Xenium Dataset Size:** At 161,000 cells, the Xenium analysis is computationally demanding. Performing the full suite of spatial statistics may require significant memory allocations (**>64 GB RAM**) and lengthy compute runtime.
* **Targeted Gene Panels:** Xenium and MERFISH use pre-designed targeted gene panels, which restricts the discovery of novel biology strictly to the genes included on the panel.
* **Visium Spot Size Mixing:** At a ~55 µm diameter, each Visium spot captures a multi-cellular mixture (typically 1–10 cells depending on tissue density). Gene expression per spot therefore reflects a bulk-like mixture of cell types rather than a singular cell resolution.
* **Pre-Annotated Clusters:** Tutorials 02, 03, and 04 rely on pre-annotated cluster labels. In an unsupervised real-world analysis, cluster annotation requires rigorous marker gene review and matching against reference single-cell atlases.
