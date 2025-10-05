# ---------- WebGestalt exports, Volcano plots, Heatmaps & PCA of significant genes ----------

# parameters for "significant"
lfc_cutoff_plot <- 1.0            # used for volcano, heatmap & PCA thresholds
top_label_n     <- 12             # number of gene labels on volcano plots

# Helpers -------------------------------------------------------------

# 1) Export files for WebGestalt (ORA list and ranked list)
export_for_webgestalt <- function(res_df, tag, out_dir, symbol_map, alpha = alpha_thresh) {
  rdf <- add_symbols_to_results(res_df, symbol_map)
  rdf <- rdf[!is.na(rdf$symbol) & nzchar(rdf$symbol), ]
  
  # ORA: one symbol per line (padj < alpha)
  ora <- unique(rdf$symbol[!is.na(rdf$padj) & rdf$padj < alpha])
  ora_file <- file.path(out_dir, paste0("WebGestalt_ORA_", tag, "_symbols.txt"))
  writeLines(ora, con = ora_file)
  
  # RANKED: two columns: GeneSymbol \t Score (use Wald statistic as rank score)
  rank_df <- rdf[!is.na(rdf$stat), c("symbol", "stat")]
  # handle duplicates by keeping the entry with the largest absolute statistic
  rank_df <- rank_df[order(abs(rank_df$stat), decreasing = TRUE), ]
  rank_df <- rank_df[!duplicated(rank_df$symbol), ]
  rank_file <- file.path(out_dir, paste0("WebGestalt_RANK_", tag, "_symbols.tsv"))
  write.table(rank_df, file = rank_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  message("Saved WebGestalt files for ", tag, ":")
  message(" - ORA:   ", ora_file)
  message(" - RANK:  ", rank_file)
}

# 2) Volcano plot with clean labeling in gene SYMBOLS
make_volcano_plot <- function(res_df, tag, out_dir, alpha = alpha_thresh, lfc_cut = lfc_cutoff_plot) {
  rdf <- add_symbols_to_results(res_df, symbol_map)
  rdf$neglog10padj <- -log10(rdf$padj)
  rdf$pass <- with(rdf, !is.na(padj) & padj < alpha & !is.na(log2FoldChange) & abs(log2FoldChange) >= lfc_cut)
  
  # choose labels: most significant among those passing both cutoffs
  label_df <- rdf[rdf$pass, ]
  label_df <- label_df[order(label_df$padj), ]
  label_df$label <- ifelse(!is.na(label_df$symbol) & nzchar(label_df$symbol), label_df$symbol, label_df$ensembl_gene_id)
  label_df <- head(label_df, top_label_n)
  
  p <- ggplot(rdf, aes(x = log2FoldChange, y = neglog10padj)) +
    geom_point(aes(color = pass), alpha = 0.6, size = 1.2, na.rm = TRUE) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
    geom_text_repel(data = label_df, aes(label = label), size = 3, max.overlaps = Inf) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "steelblue")) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste0("Volcano plot: ", tag),
      x = "log2 fold change",
      y = expression(-log[10]("adjusted p-value")),
      color = paste0("padj<", alpha, " & |LFC|>=", lfc_cut)
    )
  
  out_file <- file.path(out_dir, paste0("Volcano_", tag, ".png"))
  ggsave(out_file, p, width = 8.5, height = 6.5, dpi = 300)
  message("Saved volcano: ", out_file)
}

# 3) Heatmap for each contrast (reuses existing function but runs for all tags)
make_de_heatmap_all <- function(vsd, res_list, padj_cutoff = alpha_thresh, lfc_cutoff = lfc_cutoff_plot) {
  for (nm in names(res_list)) {
    make_de_heatmap(vsd, as.data.frame(res_list[[nm]]), nm,
                    padj_cutoff = padj_cutoff, lfc_cutoff = lfc_cutoff,
                    max_genes = 50, plot_dir = plot_dir, symbol_map = symbol_map)
  }
}

# 4) PCA on significant genes only (by contrast)
make_sig_pca <- function(vsd, res_df, tag, out_dir, padj_cut = alpha_thresh, lfc_cut = lfc_cutoff_plot) {
  keep <- !is.na(res_df$padj) & res_df$padj < padj_cut &
    !is.na(res_df$log2FoldChange) & abs(res_df$log2FoldChange) >= lfc_cut
  sig_ids <- intersect(rownames(res_df)[keep], rownames(vsd))
  if (length(sig_ids) < 3) { message("Too few significant genes for PCA: ", tag); return(invisible(NULL)) }
  vsd_sig <- vsd[sig_ids, ]
  
  p <- plotPCA(vsd_sig, intgroup = "condition") +
    ggtitle(paste0("PCA (significant genes): ", tag)) +
    theme_minimal(base_size = 12)
  
  out_file <- file.path(out_dir, paste0("PCA_sig_", tag, ".png"))
  ggsave(out_file, p, width = 8.5, height = 6.5, dpi = 300)
  message("Saved PCA (sig genes): ", out_file)
}

# Run exports and plots for every contrast ----------------------------

for (nm in names(res_list)) {
  cur_res <- as.data.frame(res_list[[nm]])
  
  # WebGestalt files
  export_for_webgestalt(cur_res, nm, deseq_dir, symbol_map, alpha = alpha_thresh)
  
  # Volcano plot
  make_volcano_plot(cur_res, nm, plot_dir, alpha = alpha_thresh, lfc_cut = lfc_cutoff_plot)
  
  # PCA on significant genes
  make_sig_pca(vsd, cur_res, nm, plot_dir, padj_cut = alpha_thresh, lfc_cut = lfc_cutoff_plot)
}

# Heatmaps for all contrasts (top significant genes, labeled by symbol)
make_de_heatmap_all(vsd, res_list, padj_cutoff = alpha_thresh, lfc_cutoff = lfc_cutoff_plot)




