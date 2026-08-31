# 🧬 Mouse Skin Aging Single-Cell RNA-seq Analysis (Replication of GSE336400 Dataset)

![R](https://img.shields.io/badge/R-4.x-blue.svg)
![Seurat](https://img.shields.io/badge/Seurat-v5.0-orange.svg)
![Harmony](https://img.shields.io/badge/Integration-Harmony-green.svg)
![License](https://img.shields.io/badge/License-MIT-lightgrey.svg)

This project builds an end-to-end bioinformatic workflow based on the mouse skin aging single-cell RNA sequencing dataset (**GSE336400**).

By comparing the skin cell atlases of mice across different age groups (**8 weeks - young, 26 weeks - adult, and 52 weeks - aged**), this study successfully identified and deeply characterized a **pro-inflammatory/pro-aging macrophage subcluster (Mac4)** that is significantly enriched in aging skin.

---

## 📌 Biological Insights

* **Mapping the Global Landscape of Skin Aging**: Systematically analyzed cellular composition dynamics in mouse skin across different age groups (8wk / 26wk / 52wk), accurately identifying major cell types including fibroblasts, keratinocytes, endothelial cells, and macrophages.
* **Discovery of an Aging-Specific Macrophage Subcluster (Mac4)**: Pinpointed a critical pro-inflammatory subcluster, **Mac4**, within the broader macrophage population (Mac1–Mac6). This subcluster exhibits an **explosive expansion in 52-week-old (aged) mouse skin**, serving as a core hallmark of aging microenvironment remodeling.
* **Characterization of Mac4 Molecular Signature**: Identified specific high expression of classic pro-inflammatory markers and biomarker genes (*S100a8*, *S100a9*, *Cd14*, and *Il1b*) in the Mac4 subcluster.
* **Elucidation of Molecular Mechanisms Driving Skin Immunosenescence**: GO functional enrichment revealed that the Mac4 subcluster strongly activates **LPS-like inflammatory response, leukocyte recruitment, and IL-1 signaling pathways**, driving chronic inflammation and immunosenescence in the skin microenvironment by continuously releasing pro-inflammatory factors to recruit immune cells.

---



## 📂 Repository Structure

```text
.
├── BioPratice10_scRNA.Rproj                # RStudio project file
├── README.md                               # Project documentation in English
├── README_CN.md                            # Project documentation in Chinese
├── data/                                   # Raw 10x Genomics scRNA-seq data directory
├── processed_data/                         # Processed data objects & outputs (.rds, .csv)
├── plots/                                  # High-resolution output figures directory
│   ├── 01_QC_vlnplot_before_filtering.png
│   ├── 01_QC_vlnplot_after_filtering.png
│   ├── 02_UMAP_global_clusters.png
│   ├── 02_UMAP_annotated.png
│   ├── 02_UMAP_split_by_age.png
│   ├── 02_Marker_genes_dotplot.png
│   ├── 03_Macrophage_UMAP_subclusters.png
│   ├── 03_Mac4_Marker_Bubble.png
│   ├── 03_Mac4_Proportion_by_Age.png
│   ├── 04_Mac4_DEG_Volcano.png
│   └── 04_Mac4_GO_Enrichment.png
└── scripts/                                # Core analysis R scripts directory
    ├── 01_qc_and_filtering.R               # Step 01: Quality control & cell filtering
    ├── 02_integration_and_clustering.R     # Step 02: Harmony integration & global annotation
    ├── 03_macrophage_subclustering.R       # Step 03: Macrophage subclustering & Mac4 locking
    └── 04_deg_and_go_enrichment.R          # Step 04: Mac4 DEG & GO enrichment analysis
```



# Result



### 1. Global Cell Landscape & Sample Integration

Integrated 8wk, 26wk, and 52wk mouse skin cells using Harmony for dimensionality reduction and batch effect removal, followed by major cell type annotation.

<table>
  <tr>
    <td align="center"><b>Global Cell Type Annotation UMAP</b></td>
    <td align="center"><b>Split-by-Age Comparison UMAP</b></td>
  </tr>
  <tr>
    <td align="center"><img src="plots/02_UMAP_annotated.png" width="100%"></td>
    <td align="center"><img src="plots/02_UMAP_split_by_age.png" width="100%"></td>
  </tr>
</table>


### 2. Macrophage Subclustering & Mac4 Lock-in

Re-clustered isolated macrophages, revealing that the *S100a8+/S100a9+* high-expressing **Mac4** subcluster drastically expands in proportion in the 52-week aged group.

<table>
  <tr>
    <td align="center"><b>Mac4 Marker Gene Expression Dotplot</b></td>
    <td align="center"><b>Macrophage Subcluster Proportion Dynamics</b></td>
  </tr>
  <tr>
    <td align="center"><img src="plots/03_Mac4_Marker_Bubble.png" width="100%"></td>
    <td align="center"><img src="plots/03_Mac4_Proportion_by_Age.png" width="100%"></td>
  </tr>
</table>


### 3. Molecular Mechanisms & GO Enrichment Analysis

The Mac4 DEG volcano plot and GO dotplot confirm that this subcluster primarily functions in pro-inflammation, cell recruitment, and IL-1 signaling in aged skin.

<table>
  <tr>
    <td align="center"><b>Mac4 DEG Volcano Plot</b></td>
    <td align="center"><b>GO Biological Process Enrichment Dotplot</b></td>
  </tr>
  <tr>
    <td align="center"><img src="plots/04_Mac4_DEG_Volcano.png" width="100%"></td>
    <td align="center"><img src="plots/04_Mac4_GO_Enrichment.png" width="100%"></td>
  </tr>
</table>



## 🚀 Quick Start & Usage

###  Data Preparation

The raw single-cell RNA-seq datasets used in this project can be obtained from the NCBI GEO database:

1. **Download Raw Data**: Visit the GEO dataset [GSE336400](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE336400) and download the supplementary files.

2. **Directory Structure**: Place and organize the downloaded 10x Genomics output files (`matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`) into the `data/` directory:
   ```text
   data/
   ├── 8wk_1/
   │   ├── matrix.mtx.gz
   │   ├── features.tsv.gz
   │   └── barcodes.tsv.gz
   ├── 26wk_1/
   │   └── ...
   └── 52wk_1/
       └── ...
### Environment Setup

Ensure R (>= 4.2) and the following core dependency packages are installed:

```
install.packages(c("tidyverse", "ggrepel"))
BiocManager::install(c("Seurat", "harmony", "clusterProfiler", "org.Mm.eg.db"))
```




### Execution Steps

Run the R scripts sequentially in the `scripts/` directory:



```
# Step 1: Data loading, quality control, and filtering
Rscript scripts/01_qc_and_filtering.R

# Step 2: Global data integration, dimension reduction, clustering, and annotation
Rscript scripts/02_integration_and_clustering.R

# Step 3: Macrophage subset extraction and Mac4 lock-in
Rscript scripts/03_macrophage_subclustering.R

# Step 4: Mac4 differentially expressed gene (DEG) & GO enrichment analysis
Rscript scripts/04_deg_and_go_enrichment.R
```



## 🛠️ Troubleshooting & Technical Pitfalls Log

### 1. Seurat V5 Unjoined Layers Error (`JoinLayers`)

- **Issue**: When running `FindMarkers()` in Script 04 to calculate Mac4 DEGs, the following error was triggered: `Error in FindMarkers.StdAssay: data layers are not joined. Please run JoinLayers`
- **Root Cause**: Seurat V5 introduced a new `Layer` mechanism. After sample integration and subset extraction (`subset`), expression matrices are split across multiple independent layers, preventing direct differential expression calculation.
- **Solution**: Explicitly execute `seu_mac <- JoinLayers(seu_mac)` to join all data layers prior to calling `FindMarkers()`.

### 2. Decoupled Customization of DimPlot Center Cluster Labels & Legend Text

- **Issue**: Desired outcome was "clean numeric cluster IDs in the center of the UMAP plot, with `Number_CellType` displayed in the right legend". Using `RenameIdents()` directly forces the center plot labels to be replaced with long English names, causing severe text overlap and crowd.
- **Root Cause**: `DimPlot(label = TRUE)` defaults to reading the Active Idents of the object; altering Idents updates both the plot labels and the legend simultaneously.
- **Solution**: Maintain Active Idents as purely numeric `seurat_clusters`, using `scale_color_discrete(labels = ...)` to strictly modify the right-side legend. Additionally, use `unname()` when writing annotations into `meta.data` to eliminate vector name conflicts (preventing the `No cell overlap between new meta data and Seurat object` error).





