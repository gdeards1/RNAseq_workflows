#Setting local work directory
setwd("C:/Users/Gabriel/Documents/30-1205430437/00_fastq/")

install.packages("BiocManager")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("Rsubread")

BiocManager::install(c("ShortRead", "Rqc", "Rsubread", "DESeq2", "EnhancedVolcano", "pheatmap"), force = TRUE)

library(ShortRead)
library(Rsubread)
library(DESeq2)
library(EnhancedVolcano)
library(Rqc)
library(pheatmap)


# Check to see how quickly my computer can handle these files.
fastq_file <- "CL1_R1_001.fastq.gz"
chunk_size <- 100000  # reads per chunk

stream <- FastqStreamer(fastq_file, n = chunk_size)

total_reads <- 0
chunk_count <- 0

tic("Streaming entire file")

repeat {
  fq_chunk <- yield(stream)
  if (length(fq_chunk) == 0) break  # end of file
  total_reads <- total_reads + length(fq_chunk)
  chunk_count <- chunk_count + 1
}

toc()
# Streaming entire file: 110.22 sec elapsed

close(stream)

cat("Total reads read:", total_reads, "\n")
cat("Chunks processed:", chunk_count, "\n")

# Total reads read: 32187686 
# Chunks processed: 322 

# Clear that all from memory
rm(list=ls())


# With that done, let's move on and set up our analysis pipeline

# performed fastQC and multiQC using bash.

# Now let's build our genomic index. 

library(Rsubread)
buildindex(
  basename = "C:/Users/Gabriel/Documents/30-1205430437/00_fastq/genome_mouse/GRCm39_index",
  reference = "C:/Users/Gabriel/Documents/30-1205430437/00_fastq/genome_mouse/Mus_musculus.GRCm39.dna.primary_assembly.fa",
  gappedIndex = TRUE
  )

library(Rsubread)

# Path to the genome index (no file extension!)
index_path <- "C:/Users/Gabriel/Documents/30-1205430437/00_fastq/genome_mouse/GRCm39_index"

# Output directory for BAM files
output_dir <- "aligned_bam"
dir.create(output_dir, showWarnings = FALSE)

# List all FASTQ files
fastq_files <- list.files(pattern = "_R[12]_001\\.fastq\\.gz$", full.names = TRUE)
fastq_files <- sort(fastq_files)

# Extract base sample names (remove _R1/_R2 suffix)
sample_names <- unique(gsub("_R[12]_001\\.fastq\\.gz$", "", basename(fastq_files)))

# Loop through each sample
for (sample in sample_names) {
  r1_file <- paste0(sample, "_R1_001.fastq.gz")
  r2_file <- paste0(sample, "_R2_001.fastq.gz")
  
  # Full paths
  r1_path <- normalizePath(r1_file, winslash = "/", mustWork = TRUE)
  r2_path <- normalizePath(r2_file, winslash = "/", mustWork = TRUE)
  output_bam <- file.path(output_dir, paste0(basename(sample), "_aligned.bam"))
  
  cat("\n▶️ Aligning:", basename(sample), "\n")
  
  align(
    index = index_path,
    readfile1 = r1_path,
    readfile2 = r2_path,
    input_format = "gzFASTQ",
    output_file = output_bam,
    nthreads = 4,
    phredOffset = 33
  )
  
  cat("✅ Finished:", basename(sample), "\n")
  gc()  # Clean up memory
}


# Set paths again
bam_dir <- "C:/Users/Gabriel/Documents/30-1205430437/00_fastq/aligned_bam"
gtf_file <- "C:/Users/Gabriel/Documents/30-1205430437/00_fastq/genome_mouse/Mus_musculus.GRCm39.110.gtf"

# List all BAM files
bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)

# Run featureCounts
fc_results <- featureCounts(
  files = bam_files,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE,
  isPairedEnd = TRUE,
  nthreads = 4
)

# Save count matrix to CSV
write.csv(fc_results$counts, file = "gene_counts_matrix.csv")

# Save assignment stats (how many reads counted/unassigned)
write.csv(fc_results$stat, file = "featurecounts_summary.csv")


# Actual Analysis

sample_table <- data.frame(
  sample = c("CL1", "CL2", "CL3",
             "IL101", "IL102", "IL103",
             "IL10CL1", "IL10CL2", "IL10CL3"),
  condition = c("control", "control", "control",
                "treatment", "treatment", "treatment",
                "treatment", "treatment", "treatment")
)

# Make condition a factor (very important!)
sample_table$condition <- factor(sample_table$condition)

# Check the table
print(sample_table)


# DESEq

# Load the counts from CSV (first column is gene names)
counts <- read.csv("gene_counts_matrix.csv", row.names = 1)

# Confirm counts column names match sample names
colnames(counts)
# Should match: CL1, CL2, CL3, IL101, ...

# Make sure the sample_table from Step 1 is loaded
# sample_table should already be in your session

# Sanity check — do names match?
all(sample_table$sample %in% colnames(counts))  # Should return TRUE


# Remove _aligned.bam from column names
colnames(counts) <- gsub("_aligned\\.bam$", "", colnames(counts))

# Check again
all(sample_table$sample %in% colnames(counts))  # Should now return TRUE

# Reorder columns to match sample_table
counts <- counts[, sample_table$sample]



# Reorder columns of counts to match sample_table
counts <- counts[, sample_table$sample]

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_table,
  design = ~ condition  # this tells DESeq2 to test for differences by condition
)


# Total number of genes before filtering
total_genes_before <- nrow(dds)

# Pre-filtering: remove rows with low counts (optional but helps performance)

# Determine which rows will be kept
keep <- rowSums(counts(dds)) >= 10

# Count how many are kept and how many are filtered out
genes_kept <- sum(keep)
genes_removed <- total_genes_before - genes_kept

# Print summary
cat("Total genes before filtering:", total_genes_before, "\n")
cat("Genes kept (≥10 total counts):", genes_kept, "\n")
cat("Genes removed (<10 total counts):", genes_removed, "\n")


# finally running DESeq2 analysis

# Run the differential expression analysis
dds <- DESeq(dds)

# View results (defaults to comparing treatment vs control)
res <- results(dds)

# View results with 0.05 p-value
res_05 <- results(dds, alpha = 0.05)


# Quick summary
summary(res)
summary(res_05)

# Order by adjusted p-value (most significant genes first)
res_ordered <- res[order(res$padj), ]

# Save to CSV
write.csv(as.data.frame(res_ordered), "deseq2_results.csv")

res_05_ordered <- res_05[order(res_05$padj), ]
write.csv(as.data.frame(res_05_ordered), "deseq2_results_padj_0.05.csv")

saveRDS(dds, file = "dds_processed.rds")




# visualization

library(EnhancedVolcano)

EnhancedVolcano(res_05,
                lab = rownames(res_05),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                title = 'Differential Expression: Treatment vs Control',
                subtitle = 'Genes with padj < 0.05 and |log2FC| > 1',
                caption = paste0("Total genes: ", nrow(res_05))
)



# making the genes readable

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

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


# Convert rownames to column for joining
res_df <- as.data.frame(res_05)
res_df$ensembl_gene_id <- rownames(res_df)

# Merge
res_annotated <- merge(res_df, annotations, by = "ensembl_gene_id")

# View annotated table
head(res_annotated)

# new dataframe
res_df <- as.data.frame(res_05)
res_df$ensembl_gene_id <- rownames(res_df)

res_annotated <- merge(res_df, annotations, by = "ensembl_gene_id")
write.csv(res_annotated, "deseq2_results_annotated.csv", row.names = FALSE)



# heatmap
vsd <- vst(dds, blind = FALSE)  # variance-stabilizing transformation

plotPCA(vsd, intgroup = "condition")


library(pheatmap)

# Extract top 50 genes by adjusted p-value
top_genes <- head(order(res_05$padj), 50)

# Get normalized expression values
mat <- assay(vsd)[top_genes, ]

# Z-score transform rows (optional but standard)
mat <- t(scale(t(mat)))

# Plot heatmap
pheatmap(mat, annotation_col = as.data.frame(colData(vsd)))

# Replace rownames of the vsd matrix with gene symbols
# First, make sure symbols are in the same order as your vst matrix
symbol_map <- annotations
rownames(symbol_map) <- symbol_map$ensembl_gene_id

# Get gene symbols in correct order
gene_symbols <- symbol_map[rownames(vsd), "external_gene_name"]
rownames(vsd) <- ifelse(is.na(gene_symbols), rownames(vsd), gene_symbols)


# enhanced heatmap:

top_genes <- head(order(res_05$padj), 50)
mat <- assay(vsd)[top_genes, ]
mat <- t(scale(t(mat)))  # z-score

pheatmap(mat, annotation_col = as.data.frame(colData(vsd)))


# volcano
EnhancedVolcano(res_annotated,
                lab = res_annotated$external_gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Volcano Plot: Treatment vs Control',
                pCutoff = 0.05,
                FCcutoff = 1.0
)


png("volcano_plot.png", width = 1200, height = 1000, res = 150)
EnhancedVolcano(res_annotated,
                lab = res_annotated$external_gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Volcano Plot: Treatment vs Control',
                pCutoff = 0.05,
                FCcutoff = 1.0
)
dev.off()


pheatmap(mat,
         annotation_col = as.data.frame(colData(vsd)),
         filename = "heatmap_top50_genes3.png",  # this saves it directly
         width = 10,
         height = 12)



EnhancedVolcano(res,
                lab = res$gene_symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'p-value'),
                title = 'Treatment vs Control',
                subtitle = 'Differential expression analysis',
                caption = paste0('Total = ', nrow(res), ' genes'),
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.5,
                labSize = 3.5,
                max.overlaps = 10,
                colAlpha = 0.6,
                col = c("grey50", "forestgreen", "royalblue", "red2"),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                gridlines.major = FALSE,
                gridlines.minor = FALSE
)



#############

#oops! CL1-3 are not control. Need to redo a lot of things.

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

total_genes_before <- nrow(dds)

# Pre-filtering: remove rows with low counts (optional but helps performance)

# Determine which rows will be kept
keep <- rowSums(counts(dds)) >= 10

# Count how many are kept and how many are filtered out
genes_kept <- sum(keep)
genes_removed <- total_genes_before - genes_kept

# Print summary
cat("Total genes before filtering:", total_genes_before, "\n")
cat("Genes kept (≥10 total counts):", genes_kept, "\n")
cat("Genes removed (<10 total counts):", genes_removed, "\n")


# View results (defaults to comparing treatment vs control)
res <- results(dds)

# View results with 0.05 p-value
res_05 <- results(dds, alpha = 0.05)

# Convert rownames to column for joining
res_df <- as.data.frame(res_05)
res_df$ensembl_gene_id <- rownames(res_df)

# Merge
res_annotated <- merge(res_df, annotations, by = "ensembl_gene_id")

# View annotated table
head(res_annotated)

# new dataframe
res_df <- as.data.frame(res_05)
res_df$ensembl_gene_id <- rownames(res_df)

# IL10 vs CL
res_IL10_vs_CL <- results(dds, contrast = c("condition", "IL10", "CL"), alpha = 0.05)

# IL10CL vs CL
res_IL10CL_vs_CL <- results(dds, contrast = c("condition", "IL10CL", "CL"), alpha = 0.05)

# IL10CL vs IL10
res_IL10CL_vs_IL10 <- results(dds, contrast = c("condition", "IL10CL", "IL10"), alpha = 0.05)



