#Setting local work directory
setwd("C:/Users/Gabriel/Documents/30-1205430437/00_fastq/")

library(ShortRead)
library(Rsubread)
library(DESeq2)
library(EnhancedVolcano)
library(Rqc)
library(pheatmap)
library(biomaRt)


ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
gene_ids <- rownames(res_05)

# Fetch gene symbols and descriptions
annotations <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name', 'description'),
  filters = 'ensembl_gene_id',
  values = gene_ids,
  mart = ensembl
)


# Analysis

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



EnhancedVolcano(res_05,
                lab = rownames(res_IL10_vs_CL),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                title = 'Differential Expression: Treatment vs Control',
                subtitle = 'Genes with padj < 0.05 and |log2FC| > 1',
                caption = paste0("Total genes: ", nrow(res_IL10_vs_CL))
)


#total_genes_before <- nrow(dds)
#keep <- rowSums(counts(dds)) >= 10

# Count how many are kept and how many are filtered out
#genes_kept <- sum(keep)
#genes_removed <- total_genes_before - genes_kept

# Print summary
#cat("Total genes before filtering:", total_genes_before, "\n")
#cat("Genes kept (≥10 total counts):", genes_kept, "\n")
#cat("Genes removed (<10 total counts):", genes_removed, "\n")


# View results with 0.05 p-value
res_05 <- results(dds, alpha = 0.05)

# Convert rownames to column for joining
res_df <- as.data.frame(res_05)
res_df$ensembl_gene_id <- rownames(res_df)

# Merge
res_annotated <- merge(res_df, annotations, by = "ensembl_gene_id")

# View annotated table
head(res_annotated)



# IL10 vs CL
res_IL10_vs_CL <- results(dds, contrast = c("condition", "IL10", "CL"), alpha = 0.05)

# IL10CL vs CL
res_IL10CL_vs_CL <- results(dds, contrast = c("condition", "IL10CL", "CL"), alpha = 0.05)

# IL10CL vs IL10
res_IL10CL_vs_IL10 <- results(dds, contrast = c("condition", "IL10CL", "IL10"), alpha = 0.05)


