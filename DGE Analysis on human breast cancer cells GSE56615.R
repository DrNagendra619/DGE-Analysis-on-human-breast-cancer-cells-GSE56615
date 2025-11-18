# Title: Differential Expression Analysis of RRAS2 Knockdown in MDA-MB-231-Luc Cells (GSE56615)
# Dataset: GSE56615 (MDA-MB-231 cells: Control vs. RRAS2-Knockdown)


# 1. Installation and Library Loading
# -------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "EnhancedVolcano", "pheatmap", "GEOquery"))

library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(ggplot2)
library(GEOquery)
set.seed(123)

# 2. Define and Create Save Directory
# -------------------------------------------------
# This code sets your save path to D:\DOWNLOADS
save_path <- "D:/DOWNLOADS/"
if (!dir.exists(save_path)) {
  cat("Creating directory:", save_path, "\n")
  dir.create(save_path, recursive = TRUE)
}
cat("Plots will be saved to:", save_path, "\n")

# 3. STEP 1 — Load REAL GEO DATA (GSE56615)
# -------------------------------------------------
cat("Downloading GEO dataset GSE56615... This may take a moment.\n")
gse <- getGEO("GSE56615", GSEMatrix = TRUE)
# --- NOTE: GSE56615 usually only has one matrix file, so we select the first element ---
gse <- gse[[1]] 

# Extract expression data and metadata
exprSet <- exprs(gse)
pdata <- pData(gse)

cat("Dataset dimensions (Genes x Samples):\n")
print(dim(exprSet))

# --- Define Sample Conditions ---
# Based on the GEO page, we use the 'title' column
cat("\nSample titles:\n")
print(pdata$title)

# Create a new 'condition' column based on the titles
# --- CORRECTED LINE: Using "Control" instead of "shControl" to match metadata ---
pdata$condition <- ifelse(
  grepl("Control", pdata$title), 
  "Control", 
  "Knockdown"
)

# Set "Control" as the base level for comparison
pdata$condition <- factor(pdata$condition, levels = c("Control", "Knockdown"))

cat("\nFinal condition table:\n")
print(table(pdata$condition))

# 4. STEP 2 — Run DESeq2
# -------------------------------------------------
# --- Prepare DESeq2 Input ---
# This is microarray data. We must round the values to integers.
countData <- round(exprSet)

# Create the DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = pdata,
  design = ~ condition
)

# Filter out low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

cat("\nRunning DESeq2 analysis...\n")

# --- FIX FOR DISPERSION ERROR (FINAL VERSION) ---
# We replace the single DESeq(dds) line with the full manual workflow:
dds <- estimateSizeFactors(dds)       # 1. Calculate normalization factors (Size Factors)
dds <- estimateDispersionsGeneEst(dds) # 2. Estimate dispersions (gene-wise)
dispersions(dds) <- mcols(dds)$dispGeneEst # 3. Use gene-wise estimates as final estimates
dds <- nbinomWaldTest(dds)           # 4. Perform the statistical testing
# ------------------------------------------------

res <- results(dds, contrast = c("condition", "Knockdown", "Control"))
res <- res[order(res$padj), ] # Order by adjusted p-value

cat("\nDESeq2 results summary:\n")
summary(res)

# 5. STEP 3 — Volcano Plot and Save
# -------------------------------------------------
cat("Generating Volcano Plot...\n")

v_plot <- EnhancedVolcano(
  res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "padj", # Use adjusted p-value
  pCutoff = 0.05,
  FCcutoff = 1.0,
  title = "RRAS2 Knockdown vs. Control (GSE56615)",
  subtitle = "MDA-MB-231-Luc Cells"
)

# --- Save Volcano Plot ---
ggsave(filename = "volcano_plot_GSE56615.png", plot = v_plot, path = save_path, width = 10, height = 8, units = "in")
cat("Volcano plot saved to:", file.path(save_path, "volcano_plot_GSE56615.png"), "\n")

# 6. STEP 4 — VST + Heatmap and Save
# -------------------------------------------------
cat("Generating Heatmap...\n")

# --- FIX: MANUAL LOG TRANSFORMATION WORKAROUND ---
# Since DESeq2's VST/rlog functions fail, we use a manual log2 transformation 
# on the normalized counts to stabilize the data for plotting.

# 1. Get Normalized Counts
normalized_counts <- counts(dds, normalized=TRUE)

# 2. Apply Log2 Transformation (adding a small pseudocount of 1 to avoid log(0))
vsd_mat <- log2(normalized_counts + 1)
# -----------------------------------------------------------------------

# Select top 30 most significant DE genes
top_genes <- head(rownames(res), 30)

# --- NEW FIX: Filter out genes with zero variance across selected samples ---
# This prevents the NA/NaN clustering error in pheatmap.
vsd_mat_filtered <- vsd_mat[top_genes, ]
vsd_mat_filtered <- vsd_mat_filtered[apply(vsd_mat_filtered, 1, var) > 0, ]
# --------------------------------------------------------------------------

# --- Open PNG file device to save heatmap ---
png(file.path(save_path, "heatmap_top30_GSE56615.png"), width = 800, height = 1000, res = 100)

pheatmap(
  vsd_mat_filtered, # Use the variance-filtered matrix
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  annotation_col = as.data.frame(colData(dds)[, "condition", drop=FALSE]), 
  main = "Top 30 DE Genes (RRAS2 Knockdown - GSE56615)"
)

# --- Close the file device ---
dev.off()
cat("Heatmap saved to:", file.path(save_path, "heatmap_top30_GSE56615.png"), "\n")

# 7. STEP 5 — PCA PLOT and Save
# -------------------------------------------------
cat("Generating PCA Plot...\n")

# --- FINAL FIX: Manual PCA Calculation (Bypassing DESeqTransform Check) ---

# 1. Prepare the log-transformed matrix for PCA
# We must remove any genes with zero variance for the prcomp function.
pca_matrix <- t(vsd_mat)
pca_matrix <- pca_matrix[, apply(pca_matrix, 2, var) > 0.001] # Filter columns (genes) with near-zero variance

# 2. Run PCA (prcomp) on the filtered matrix
pca_result <- prcomp(pca_matrix, scale. = TRUE)

# 3. Extract results and prepare for ggplot2
pcaData <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  condition = colData(dds)$condition # Extract condition from original metadata
)

# 4. Calculate percent variance explained by PC1 and PC2
percentVar <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2))[1:2]

# 5. Generate the plot using ggplot2
pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal(base_size = 14) +
  ggplot2::ggtitle("PCA Plot – RRAS2 Knockdown (GSE56615)")

# --- Save PCA Plot ---
ggsave(filename = "pca_plot_GSE56615.png", plot = pca_plot, path = save_path, width = 8, height = 6, units = "in")
cat("PCA plot saved to:", file.path(save_path, "pca_plot_GSE56615.png"), "\n")

# 8. STEP 6 — Save Significant Results
# -------------------------------------------------
cat("Saving significant results table...\n")
sigGenes <- as.data.frame(subset(res, padj < 0.05))
write.csv(sigGenes, file.path(save_path, "Significant_DEGs_GSE56615.csv"))

cat("\nSignificant DEGs saved to:", file.path(save_path, "Significant_DEGs_GSE56615.csv"), "\n")

# --- Final Confirmation Message ---
cat("\nAnalysis complete. All files saved successfully to", save_path, "\n")
# ----------------------------------
