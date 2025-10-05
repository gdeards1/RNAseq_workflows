sample_table <- data.frame(
  sample = c("CL1", "CL2", "CL3",
             "IL101", "IL102", "IL103",
             "IL10CL1", "IL10CL2", "IL10CL3"),
  condition = c("CL", "CL", "CL",
                "IL10", "IL10", "IL10",
                "IL10CL", "IL10CL", "IL10CL")
)
rownames(sample_table) <- sample_table$sample

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_table,
  design = ~ condition
)

dds <- DESeq(dds)

# IL10 vs CL
res_IL10_vs_CL <- results(dds, contrast = c("condition", "IL10", "CL"), alpha = 0.05)

# IL10CL vs CL
res_IL10CL_vs_CL <- results(dds, contrast = c("condition", "IL10CL", "CL"), alpha = 0.05)

# IL10CL vs IL10
res_IL10CL_vs_IL10 <- results(dds, contrast = c("condition", "IL10CL", "IL10"), alpha = 0.05)

# 1. Perform Variance Stabilizing Transformation
library(DESeq2)
library(matrixStats)
library(biomaRt)

vst_result <- vst(dds, blind = FALSE)

# 2. Prepare annotation for heatmap
vst_mat <- assay(vst_result)
annotation_col <- sample_table
rownames(annotation_col) <- colnames(vst_mat)

# 3. Identify top 50 most variable genes
top_genes <- head(order(rowVars(vst_mat), decreasing = TRUE), 50)
vst_mat_top <- vst_mat[top_genes, ]

# Annotate genes
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
gene_ids <- rownames(vst_mat_top)
annotations <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name'),
  filters = 'ensembl_gene_id',
  values = gene_ids,
  mart = ensembl
)
rownames(annotations) <- annotations$ensembl_gene_id
rownames(vst_mat_top) <- annotations[rownames(vst_mat_top), 'external_gene_name']

# 4. Plot heatmap
library(pheatmap)
library(RColorBrewer)
group_colors <- list(
  condition = c(CL = "#4B9CD3", IL10 = "#E69F00", IL10CL = "#009E73")
)

pheatmap(
  vst_mat_top,
  annotation_col = annotation_col,
  annotation_colors = group_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 10,
  color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
  main = "Top 50 Most Variable Genes",
  filename = "heatmap_top50_genes_New.png",  # this saves it directly
  width = 10,
  height = 12)
)

# Annotate result tables for volcano plots
annotate_res <- function(res_obj) {
  gene_ids <- rownames(res_obj)
  annotations <- getBM(
    attributes = c('ensembl_gene_id', 'external_gene_name', 'description'),
    filters = 'ensembl_gene_id',
    values = gene_ids,
    mart = ensembl
  )
  res_df <- as.data.frame(res_obj)
  res_df$ensembl_gene_id <- rownames(res_df)dg
  res_annotated <- merge(res_df, annotations, by = "ensembl_gene_id")
  return(res_annotated)
}

res_IL10_vs_CL_annotated <- annotate_res(res_IL10_vs_CL)
res_IL10CL_vs_CL_annotated <- annotate_res(res_IL10CL_vs_CL)
res_IL10CL_vs_IL10_annotated <- annotate_res(res_IL10CL_vs_IL10)

# Save annotated results to CSV files
write.csv(res_IL10_vs_CL_annotated, "res_IL10_vs_CL_annotated.csv", row.names = FALSE)
write.csv(res_IL10CL_vs_CL_annotated, "res_IL10CL_vs_CL_annotated.csv", row.names = FALSE)
write.csv(res_IL10CL_vs_IL10_annotated, "res_IL10CL_vs_IL10_annotated.csv", row.names = FALSE)

# 5. Volcano plots
library(EnhancedVolcano)

png("IL10_vs_CL.png", width = 1200, height = 1000, res = 150)
EnhancedVolcano(res_IL10_vs_CL_annotated,
                lab = res_IL10_vs_CL_annotated$external_gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'IL10 vs CL',
                pCutoff = 0.05,
                FCcutoff = 1.0
)
dev.off()

png("IL10CL_vs_CL.png", width = 1200, height = 1000, res = 150)
EnhancedVolcano(res_IL10CL_vs_CL_annotated,
                lab = res_IL10CL_vs_CL_annotated$external_gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'IL10CL vs CL',
                pCutoff = 0.05,
                FCcutoff = 1.0
)
dev.off()

png("IL10CL_vs_IL10.png", width = 1200, height = 1000, res = 150)
EnhancedVolcano(res_IL10CL_vs_IL10_annotated,
                lab = res_IL10CL_vs_IL10_annotated$external_gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'IL10CL vs IL10',
                pCutoff = 0.05,
                FCcutoff = 1.0
)
dev.off()
