# DGE-Analysis-on-human-breast-cancer-cells-GSE56615
DGE Analysis on human breast cancer cells GSE56615
# 🧬 DESeq2 DGE Pipeline: RRAS2 Knockdown in Breast Cancer (GSE56615)

This R script automates a robust pipeline for **Differential Gene Expression (DGE)** analysis of the **GSE56615** microarray dataset. The study investigates the transcriptional effects of **RRAS2 gene knockdown** in highly aggressive **MDA-MB-231-Luc human breast cancer cells**.

The pipeline adapts the **`DESeq2`** package for robust statistical testing on microarray data and includes critical workarounds for accurate visualization and quality control (QC).

## 🚀 Key Features

* **Advanced DESeq2 Application:** Runs a custom, four-step `DESeq2` workflow (size factor estimation, gene-wise dispersion estimation, nbinom Wald test) tailored for microarray data where typical `DESeq2` methods may encounter convergence issues.
* **Automated Data Retrieval:** Downloads expression data and metadata directly from **GEO (GSE56615)**.
* **Handling Continuous Data:** Correctly **rounds continuous microarray data to integers** before analysis, as required by `DESeq2`.
* **Critical Workarounds:** Implements manual `log2` transformation and **manual PCA calculation** to successfully generate QC plots despite transformation challenges often found when using `DESeq2` on non-count data.
* **Integrated Visualization:** Generates essential plots for interpretation: **Volcano Plot**, **Heatmap**, and **PCA Plot**.

---

## 🔬 Analysis Overview

| Component | Method / Test | Purpose |
| :--- | :--- | :--- |
| **Dataset** | GSE56615 | RRAS2 Knockdown in MDA-MB-231-Luc human breast cancer cells. |
| **DGE Tool** | `DESeq2` (Adapted Workflow) | Statistical method for robust comparison of gene expression levels. |
| **Comparison** | RRAS2 Knockdown vs. Control | Identifies genes regulated by RRAS2 in the context of breast cancer. |
| **Data Transformation** | $\log_2(\text{Normalized Counts} + 1)$ | Manual stabilization for visualization plots (Heatmap, PCA). |
| **Significance** | $\text{padj} < 0.05$, $|\text{log2FC}| > 1.0$ | Used for filtering and highlighting results. |

---

## 🛠️ Prerequisites and Setup

### 📦 Packages

The script automatically installs and loads the necessary Bioconductor and CRAN packages:
* `DESeq2`
* `EnhancedVolcano`
* `pheatmap`
* `GEOquery`
* `ggplot2`

### ⚙️ Execution

1.  **Download** the `DGE Analysis on human breast cancer cells GSE56615.R` file.
2.  **Optional:** The output path is set to `D:/DOWNLOADS/` by default (Step 2). You can change this path if needed.
3.  **Execute** the script in your R environment:
    ```R
    source("DGE Analysis on human breast cancer cells GSE56615.R")
    ```

---

## 📁 Output Files (3 Plots + 1 CSV)

All output files are saved to the specified `save_path` (default: `D:/DOWNLOADS/`).

### Statistical Results

| Filename | Type | Description |
| :--- | :--- | :--- |
| `Significant_DEGs_GSE56615.csv` | CSV | Table containing all genes with an adjusted p-value (padj) $< 0.05$. |

### Visualization and QC Plots

| Filename | Analysis Stage | Description |
| :--- | :--- | :--- |
| `volcano_plot_GSE56615.png` | Results | **Volcano Plot** showing $\log_2 \text{Fold Change}$ vs. $\log_{10}(\text{padj})$, highlighting DEGs. |
| `heatmap_top30_GSE56615.png` | Results | **Heatmap** of the **Top 30 Most Significant DE Genes** (using manual $\log_2$ transformed data). |
| `pca_plot_GSE56615.png` | QC / Results | **Principal Component Analysis (PCA)** plot demonstrating sample clustering and the effect of RRAS2 knockdown. |
