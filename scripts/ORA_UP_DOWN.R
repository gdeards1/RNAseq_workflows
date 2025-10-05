# Use one of your annotated results (e.g., IL10 vs CL)
ranked_genes <- res_IL10_vs_CL_annotated

# Remove NAs and keep gene symbols and log2FC
ranked_genes <- ranked_genes[!is.na(ranked_genes$log2FoldChange) & !is.na(ranked_genes$external_gene_name), ]
ranked_genes <- ranked_genes[!duplicated(ranked_genes$external_gene_name), ]

# Create named vector
ranked_list <- setNames(ranked_genes$log2FoldChange, ranked_genes$external_gene_name)

# Sort from high to low (optional: sort by statistic instead)
ranked_list <- sort(ranked_list, decreasing = TRUE)

# Export as a 2-column file (gene + rank)
ranked_df <- data.frame(Gene=names(ranked_list), Score=ranked_list)
write.table(ranked_df, file = "IL10_vs_CL_ranked_list.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



plotPCA(vst_result, intgroup = "condition")
dev.off()






# Function to export up- and down-regulated genes
export_up_down_ORA <- function(res_df, prefix, lfc_threshold = 0, padj_cutoff = 0.05) {
  # Filter valid gene symbols
  res_df <- res_df[!is.na(res_df$external_gene_name), ]
  
  # Upregulated
  up <- subset(res_df, padj < padj_cutoff & log2FoldChange >  lfc_threshold)
  writeLines(unique(up$external_gene_name), paste0(prefix, "_UP_ORA.txt"))
  
  # Downregulated
  down <- subset(res_df, padj < padj_cutoff & log2FoldChange < -lfc_threshold)
  writeLines(unique(down$external_gene_name), paste0(prefix, "_DOWN_ORA.txt"))
}

# Apply to each comparison
export_up_down_ORA(res_IL10_vs_CL_annotated,     "IL10_vs_CL-all")
export_up_down_ORA(res_IL10CL_vs_CL_annotated,   "IL10CL_vs_CL-all")
export_up_down_ORA(res_IL10CL_vs_IL10_annotated, "IL10CL_vs_IL10-al")