# ============================================
# Script: 05_featurecounts_DESeq2.R
# Purpose: Quantify genes (featureCounts), run DESeq2, add gene symbols, plot
# Inputs:  *_aligned.bam    (from 02_alignment.R)
# Outputs: counts w/ symbols, DE tables w/ symbols, PCA/heatmap (symbols), session info
# ============================================

suppressPackageStartupMessages({
  library(Rsubread)
  library(DESeq2)
  library(pheatmap)
  library(ggplot2)
  library(biomaRt)
  library(matrixStats)
  library(ggrepel)
})

# ---------- Paths / Params ----------
base_dir <- "C:/Users/Gabriel/Documents/RProjects/RNAseq_workflows"
#b_dir <- "D:/BST2/30-1016475569"
bam_dir    <- file.path("D:/BST2/30-1016475569/00_fastq/trimmed/aligned_bam")
gtf_file   <- "C:/Users/Gabriel/Documents/30-1205430437/00_fastq/genome_mouse/Mus_musculus.GRCm39.110.gtf"

res_dir   <- file.path(base_dir, "results/02_results")
fc_dir    <- file.path(res_dir, "featureCounts")
deseq_dir <- file.path(res_dir, "DESeq2")
plot_dir  <- file.path(deseq_dir, "plots")
wg_dir    <- file.path(deseq_dir, "webgestalt_ORA")  # <— NEW: clean WebGestalt exports here
dir.create(fc_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(deseq_dir,recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(wg_dir,   recursive = TRUE, showWarnings = FALSE)

threads         <- 8
strand_specific <- 0   # 0=unstranded, 2=reverse for dUTP (change if needed)
alpha_thresh    <- 0.05
primary_contrast_tag <- "KO_vs_WT"  # single heatmap will be made for this tag

# ---------- Helpers ----------
symbol_cache_csv <- file.path(fc_dir, "gene_symbol_map.csv")
entrez_cache_csv <- file.path(fc_dir, "gene_entrez_map.csv")  # <— NEW: cache for Entrez IDs

get_symbol_map <- function(gene_ids) {
  # Try cache first
  if (file.exists(symbol_cache_csv)) {
    cache <- read.csv(symbol_cache_csv, stringsAsFactors = FALSE)
    cache_map <- setNames(cache$symbol, cache$ensembl_gene_id)
  } else {
    cache_map <- setNames(character(0), character(0))
  }
  
  need <- setdiff(gene_ids, names(cache_map))
  if (length(need) > 0) {
    message("Fetching ", length(need), " new symbols from Ensembl...")
    ens <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
    ann <- getBM(attributes = c("ensembl_gene_id","external_gene_name"),
                 filters = "ensembl_gene_id", values = need, mart = ens)
    new_map <- setNames(ann$external_gene_name, ann$ensembl_gene_id)
    
    # merge with cache and write back
    all_ids <- unique(c(names(cache_map), names(new_map)))
    all_map <- setNames(character(length(all_ids)), all_ids)
    all_map[names(cache_map)] <- cache_map
    all_map[names(new_map)]   <- new_map
    
    out <- data.frame(ensembl_gene_id = names(all_map),
                      symbol = unname(all_map),
                      stringsAsFactors = FALSE)
    write.csv(out, symbol_cache_csv, row.names = FALSE)
    return(all_map)
  }
  cache_map
}

# NEW: Entrez map helper (cached) — gives a single consistent ID type if you prefer Entrez
get_entrez_map <- function(gene_ids) {
  if (file.exists(entrez_cache_csv)) {
    cache <- read.csv(entrez_cache_csv, stringsAsFactors = FALSE)
    cache_map <- setNames(cache$entrez_id, cache$ensembl_gene_id)
  } else {
    cache_map <- setNames(character(0), character(0))
  }
  need <- setdiff(gene_ids, names(cache_map))
  if (length(need) > 0) {
    message("Fetching ", length(need), " new Entrez IDs from Ensembl...")
    ens <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
    ann <- getBM(attributes = c("ensembl_gene_id","entrezgene_id"),
                 filters = "ensembl_gene_id", values = need, mart = ens)
    ann$entrezgene_id <- as.character(ann$entrezgene_id)
    new_map <- setNames(ann$entrezgene_id, ann$ensembl_gene_id)
    all_ids <- unique(c(names(cache_map), names(new_map)))
    all_map <- setNames(character(length(all_ids)), all_ids)
    all_map[names(cache_map)] <- cache_map
    all_map[names(new_map)]   <- new_map
    out <- data.frame(ensembl_gene_id = names(all_map),
                      entrez_id = unname(all_map),
                      stringsAsFactors = FALSE)
    write.csv(out, entrez_cache_csv, row.names = FALSE)
    return(all_map)
  }
  cache_map
}

make_de_heatmap <- function(vsd, res_df, tag, padj_cutoff, lfc_cutoff, max_genes, plot_dir, symbol_map) {
  if (is.null(res_df)) { message("No results for ", tag); return(invisible(NULL)) }
  
  keep <- !is.na(res_df$padj) & res_df$padj < padj_cutoff &
    !is.na(res_df$log2FoldChange) & abs(res_df$log2FoldChange) > lfc_cutoff
  sig  <- res_df[keep, , drop = FALSE]
  sig  <- sig[order(sig$padj, na.last = NA), , drop = FALSE]
  if (nrow(sig) == 0) {
    sig <- res_df[order(res_df$padj, na.last = NA), , drop = FALSE]
  }
  
  top <- head(rownames(sig), 50)
  top <- top[top %in% rownames(vsd)]
  if (length(top) < 2) { message("Too few genes to plot for ", tag); return(invisible(NULL)) }
  
  mat <- assay(vsd)[top, ]
  mat <- t(scale(t(mat)))
  ann <- as.data.frame(colData(vsd)[, "condition", drop = FALSE])
  
  disp <- ifelse(top %in% names(symbol_map) & nzchar(symbol_map[top]), symbol_map[top], top)
  rownames(mat) <- disp
  
  pal <- colorRampPalette(c("navy","white","firebrick3"))(100)
  outfile <- file.path(plot_dir, paste0("heatmap_top", length(top), "_DE_", tag, "_symbols.png"))
  png(outfile, width = 1400, height = 1600, res = 150)
  pheatmap(mat, annotation_col = ann, color = pal,
           cluster_rows = TRUE, cluster_cols = TRUE,
           show_rownames = TRUE, fontsize_row = 7, fontsize_col = 9,
           main = paste0("Top DE genes (", tag, "): padj<", padj_cutoff, ", |LFC|>", lfc_cutoff))
  dev.off()
  message("Saved heatmap: ", outfile)
}

add_symbols_to_results <- function(res_df, symbol_map) {
  out <- as.data.frame(res_df)
  out$ensembl_gene_id <- rownames(out)
  out$symbol <- ifelse(out$ensembl_gene_id %in% names(symbol_map) & nzchar(symbol_map[out$ensembl_gene_id]),
                       symbol_map[out$ensembl_gene_id], out$ensembl_gene_id)
  out
}

# ---------- Collect BAMs ----------
bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)
stopifnot(length(bam_files) > 0)
bam_names <- sub("_aligned\\.bam$", "", basename(bam_files))

# ---------- featureCounts ----------
message("Running featureCounts...")
fc <- featureCounts(files = bam_files,
                    annot.ext = gtf_file,
                    isGTFAnnotationFile = TRUE,
                    GTF.featureType = "exon",
                    GTF.attrType = "gene_id",
                    useMetaFeatures = TRUE,
                    isPairedEnd = TRUE,
                    nthreads = threads,
                    strandSpecific = strand_specific)

counts_mat <- fc$counts
colnames(counts_mat) <- bam_names
write.csv(counts_mat, file = file.path(fc_dir, "gene_counts_matrix.csv"))
write.csv(fc$stat,     file = file.path(fc_dir, "featurecounts_summary.csv"))
saveRDS(fc,            file = file.path(fc_dir, "featurecounts_object.rds"))

# Cache symbols for ALL genes we just counted
symbol_map <- get_symbol_map(rownames(counts_mat))

# Also save an annotated counts table for convenience
counts_annot <- data.frame(
  ensembl_gene_id = rownames(counts_mat),
  symbol = ifelse(rownames(counts_mat) %in% names(symbol_map) & nzchar(symbol_map[rownames(counts_mat)]),
                  symbol_map[rownames(counts_mat)],
                  rownames(counts_mat)),
  counts_mat,
  check.names = FALSE
)
write.csv(counts_annot, file = file.path(fc_dir, "gene_counts_with_symbols.csv"), row.names = FALSE)

# ---------- Sample table ----------
cond <- sub("^([A-Za-z0-9]+-NE)-\\d+$", "\\1", bam_names)
cond <- ifelse(grepl("-NE", cond), cond, sub("^([A-Za-z0-9]+)-\\d+$", "\\1", bam_names))
sample_table <- data.frame(sample = bam_names, condition = cond, row.names = bam_names, check.names = FALSE)

# ---------- DESeq2 ----------
dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = sample_table, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

vsd <- vst(dds, blind = FALSE)
saveRDS(dds, file = file.path(deseq_dir, "dds.rds"))
saveRDS(vsd, file = file.path(deseq_dir, "vsd.rds"))

# ---------- PCA ----------
png(file.path(plot_dir, "PCA_condition.png"), width = 1200, height = 1000, res = 150)
print(plotPCA(vsd, intgroup = "condition") + ggtitle("PCA by condition"))
dev.off()

# ---------- Contrasts ----------
res_list <- list()
run_contrast <- function(groupA, groupB, tag) {
  if (all(c(groupA, groupB) %in% levels(colData(dds)$condition))) {
    r <- results(dds, contrast = c("condition", groupA, groupB), alpha = alpha_thresh)
    r <- r[order(r$padj, na.last = NA), ]
    res_list[[tag]] <<- r
    # Save with symbols
    r_sym <- add_symbols_to_results(r, symbol_map)
    write.csv(r_sym, file = file.path(deseq_dir, paste0("DE_", tag, "_with_symbols.csv")), row.names = FALSE)
  }
}

dds$condition <- droplevels(dds$condition)
run_contrast("KO",    "WT",    "KO_vs_WT")
run_contrast("OE",    "WT",    "OE_vs_WT")
run_contrast("KO-NE", "WT-NE", "KO-NE_vs_WT-NE")
run_contrast("OE-NE", "WT-NE", "OE-NE_vs_WT-NE")

# ---------- Single DE heatmap (symbols) ----------
res_primary <- res_list[[primary_contrast_tag]]
make_de_heatmap(vsd, as.data.frame(res_primary), primary_contrast_tag,
                padj_cutoff = alpha_thresh, lfc_cutoff = 1.0, max_genes = 50,
                plot_dir = plot_dir, symbol_map = symbol_map)

# --------- Quick MA plots (base + improved labels) ----------
for (nm in names(res_list)) {
  png(file.path(plot_dir, paste0("MA_", nm, ".png")), width = 1200, height = 900, res = 150)
  plotMA(
    res_list[[nm]],
    main = paste("MA plot:", nm),
    xlab = "Mean of normalized counts (log scale)",
    ylab = "Log2 fold change"
  )
  abline(h = c(-1, 0, 1), col = c("red", "black", "red"), lty = c(2, 1, 2))
  dev.off()
}

for (nm in names(res_list)) {
  res_df <- as.data.frame(res_list[[nm]])
  res_df$ensembl_gene_id <- rownames(res_df)
  res_df$symbol <- ifelse(res_df$ensembl_gene_id %in% names(symbol_map) & nzchar(symbol_map[res_df$ensembl_gene_id]),
                          symbol_map[res_df$ensembl_gene_id], res_df$ensembl_gene_id)
  res_df$significant <- ifelse(!is.na(res_df$padj) & res_df$padj < alpha_thresh, "Significant", "Not significant")
  top_genes <- head(res_df[order(res_df$padj), ], 10)
  p <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = significant)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_x_log10() +
    geom_hline(yintercept = 0, color = "black") +
    geom_hline(yintercept = c(-1, 1), color = "red", linetype = "dashed") +
    geom_text_repel(data = top_genes, aes(label = symbol), size = 3, color = "black") +
    theme_minimal() +
    labs(
      title = paste("MA plot:", nm),
      x = "Mean of normalized counts (log10 scale)",
      y = "log2 fold change"
    ) +
    scale_color_manual(values = c("gray", "blue"))
  ggsave(file.path(plot_dir, paste0("MA_", nm, "_ggplot.png")), p, width = 8, height = 6, dpi = 300)
}

# ---------- NEW: Clean WebGestalt ORA exports (single ID type per file) ----------
export_webgestalt_ora_lists <- function(res_list, symbol_map, out_dir = wg_dir, padj_cutoff = alpha_thresh) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  for (nm in names(res_list)) {
    r <- as.data.frame(res_list[[nm]])
    r <- r[!is.na(r$padj) & r$padj < padj_cutoff, , drop = FALSE]
    if (nrow(r) == 0) next
    ens_ids <- unique(rownames(r))
    # ENSEMBL: always consistent and accepted — use this if auto-detect fails
    ensembl_file <- file.path(out_dir, paste0("WebGestalt_ORA_", nm, "_ensembl.txt"))
    writeLines(ens_ids, con = ensembl_file, useBytes = TRUE)
    # SYMBOLS: drop any rows where symbol is missing/empty or still looks like an ENSMUSG*
    syms <- unname(symbol_map[ens_ids])
    syms <- unique(syms[!is.na(syms) & nzchar(syms) & !grepl("^ENSMUSG", syms)])
    symbol_file <- file.path(out_dir, paste0("WebGestalt_ORA_", nm, "_symbols.txt"))
    writeLines(syms, con = symbol_file, useBytes = TRUE)
    # OPTIONAL: Entrez IDs (often safest for enrichment tools)
    entrez_map <- get_entrez_map(ens_ids)
    entrez <- unique(unname(entrez_map[ens_ids]))
    entrez <- entrez[!is.na(entrez) & nzchar(entrez)]
    if (length(entrez)) {
      entrez_file <- file.path(out_dir, paste0("WebGestalt_ORA_", nm, "_entrez.txt"))
      writeLines(entrez, con = entrez_file, useBytes = TRUE)
    }
  }
  message("Saved clean WebGestalt ORA lists to: ", out_dir,
          " ( *_ensembl.txt , *_symbols.txt , and *_entrez.txt where available )")
}

export_webgestalt_ora_lists(res_list, symbol_map)

# ---------- Session info ----------
sink(file.path(deseq_dir, "sessionInfo.txt")); print(sessionInfo()); sink()

message("Done. Key outputs:")
message(" - Counts (raw):   ", file.path(fc_dir, "gene_counts_matrix.csv"))
message(" - Counts+symbols: ", file.path(fc_dir, "gene_counts_with_symbols.csv"))
message(" - DE tables w/ symbols in: ", deseq_dir)
message(" - WebGestalt ORA lists: ", wg_dir)
message(" - Plots: ", plot_dir)
