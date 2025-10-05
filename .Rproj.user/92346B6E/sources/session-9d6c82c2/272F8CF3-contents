# ============================================
# Script: 02_alignment.R
# Purpose: Align trimmed paired-end RNA-seq reads
# Author: Gabriel
# Date: YYYY-MM-DD
# Notes:
#   - Input: *_R1_001.fastq_trimmed.fastq.gz + *_R1_001.fastq_R2_trimmed.fastq.gz
#   - Output: *_aligned.bam in aligned_bam/
#   - Aligner: Rsubread
# ============================================

library(Rsubread)

input_dir  <- "D:/BST2/30-1016475569/00_fastq/trimmed"
index_path <- "C:/Users/Gabriel/Documents/30-1205430437/00_fastq/genome_mouse/GRCm39_index"
out_dir    <- file.path(input_dir, "aligned_bam")
dir.create(out_dir, showWarnings = FALSE)

# Find R1 files
r1s <- list.files(input_dir, pattern = "_R1_001\\.fastq_trimmed\\.fastq\\.gz$", full.names = TRUE)

# Derive sample names
samples <- sub("_R1_001\\.fastq_trimmed\\.fastq\\.gz$", "", basename(r1s))

run_align <- function(r1_file) {
  sample <- sub("_R1_001\\.fastq_trimmed\\.fastq\\.gz$", "", basename(r1_file))
  
  # Construct R2 filename
  r2_file <- sub("_R1_001\\.fastq_trimmed\\.fastq\\.gz$", "_R1_001.fastq_R2_trimmed.fastq.gz", r1_file)
  
  # Output BAM
  bam <- file.path(out_dir, paste0(sample, "_aligned.bam"))
  
  cat("▶️ Aligning sample:", sample, "\n")
  
  Rsubread::align(index=index_path,
                  readfile1=r1_file,
                  readfile2=r2_file,
                  input_format="gzFASTQ",
                  output_file=bam,
                  nthreads=8,
                  phredOffset=33)
  return(bam)
}

# Run all alignments sequentially (8 threads each)
bam_files <- vapply(r1s, run_align, character(1))
