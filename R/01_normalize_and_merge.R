#!/usr/bin/env Rscript
# 01_normalize_and_merge.R   —   STAGE 1 of the reproducible pipeline
#
# Load the two raw time courses, NORMALIZE each to a common count-like scale,
# test each for significance across time, and MERGE them onto a shared set of
# transcripts. The output feeds the alignment stage (02_align_timecourse.R).
#
#   Scott RNA-seq (E15,P1,P7,P14)
#       htseq counts --DESeq2 size-factor normalize--> normalized counts
#                    --DESeq2 contrasts (P7-P14, P1-P14, E15-P7)--> min adj p
#   Hatten TRAP microarray (P0,P7,P12,P18,P21,P56)
#       RMA (log2) --2^--> linear intensity
#                  --limma adj p (precomputed) --> min adj p
#   merge on Ensembl gene<->transcript annotation --> shared transcripts
#
# OUTPUT: data/normalized_merged.RData  containing
#   RNAseq_list  list(scott = normalized counts, hatten = linear intensity),
#                both matrices on the same shared transcript rows
#   meta_list    per-replicate time labels for each dataset
#   pval_all     annotation + per-dataset significance (min adjusted p-value)
#
# Usage:  Rscript R/01_normalize_and_merge.R

suppressWarnings(suppressMessages({
  library(DESeq2); library(matrixStats)
  library(stringr); library(dplyr); library(tibble)
  args_full <- commandArgs(trailingOnly = FALSE)
  here <- dirname(normalizePath(sub("^--file=", "",
            grep("^--file=", args_full, value = TRUE)[1])))
  source(file.path(here, "lib", "common.R"))
  mt_source_lib()
}))

read_cutoff <- 5   # keep genes with >=5 normalized counts in >=3 samples (for DESeq)

# ---- 1. Scott RNA-seq counts ----------------------------------------------
message("[1/4] Loading Scott RNA-seq counts")
scott_RNAseq <- read.table(MT_RAW_SCOTT_COUNTS, sep = "\t", header = TRUE, row.names = 1)
scott_RNAseq$X <- NULL
colnames(scott_RNAseq) <- matrix(unlist(strsplit(colnames(scott_RNAseq), "[.]")),
                                 ncol = 3, byrow = TRUE)[, 2]
scott_time_datasets <- read.table(MT_RAW_SCOTT_TIMES, sep = "\t", header = FALSE)
scott_time_datasets <- gsub(".*[.]", "", scott_time_datasets$V1)
scott_RNAseq <- scott_RNAseq[1:(nrow(scott_RNAseq) - 5), ]   # drop htseq summary tail
scott_GNP <- scott_RNAseq[, colnames(scott_RNAseq) %in% scott_time_datasets]
scott_GNP <- scott_GNP[, match(scott_time_datasets, colnames(scott_GNP))]

# ---- 2. DESeq2 size-factor normalization + time contrasts -----------------
message("[2/4] DESeq2 normalization and contrasts")
scott_GNP <- as.matrix(scott_GNP)
colData <- data.frame(cols = colnames(scott_GNP), condition = MT_SCOTT_META)
dds <- DESeqDataSetFromMatrix(scott_GNP, colData, design = ~ condition)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
scott_GNP_counts <- counts(dds, normalized = TRUE)   # <- normalized expression

idx <- rowSums(counts(dds, normalized = TRUE) >= read_cutoff) >= 3
dds <- dds[idx, ]
dds <- DESeq(dds)
mkres <- function(a, b) {
  r <- results(dds, contrast = c("condition", a, b))
  colnames(r) <- str_c(a, "_", b, "_", colnames(r)); r
}
scott_gnp_pval <- cbind(mkres("P7","P14"), mkres("P1","P14"), mkres("E15","P7"))
scott_gnp_pval$scott_gnp_min_pval <-
  rowMins(as.matrix(scott_gnp_pval[, grep("padj", colnames(scott_gnp_pval))]), na.rm = TRUE)

# ---- 3. Hatten TRAP microarray --------------------------------------------
message("[3/4] Loading Hatten TRAP microarray")
hatten_trap <- read.table(MT_RAW_HATTEN_EXPR, sep = "\t", header = TRUE, row.names = 1)
hatten_trap_pval <- read.table(MT_RAW_HATTEN_PVAL, sep = "\t", header = TRUE, row.names = 1)
hatten_trap <- 2^(hatten_trap)                        # RMA values are log2 -> linear
hatten_trap_pval$hatten_trap_min_pval <-
  rowMins(as.matrix(hatten_trap_pval[, grep("adj.P.Val", colnames(hatten_trap_pval))]), na.rm = TRUE)

# ---- 4. Merge on annotation (gene <-> transcript) -------------------------
message("[4/4] Merging Scott + Hatten on shared annotation")
anno <- read.table(MT_RAW_ANNO, sep = "\t", header = TRUE)
colnames(anno)[1] <- "ensembl_transcript_id"

common_genes <- intersect(anno$ensembl_gene_id, row.names(hatten_trap))
hatten_trap <- hatten_trap[row.names(hatten_trap) %in% common_genes, ]
anno <- anno %>%
  dplyr::filter(ensembl_gene_id %in% common_genes, !is.na(ensembl_gene_id)) %>%
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)

hatten_ids <- data.frame(ensembl_gene_id = row.names(hatten_trap))
hatten_ids$order <- row.names(hatten_ids)
hatten_gene2tx <- hatten_ids %>%
  dplyr::left_join(anno, by = "ensembl_gene_id") %>%
  dplyr::filter(!is.na(ensembl_transcript_id)) %>%
  dplyr::select(ensembl_gene_id, ensembl_transcript_id, order)
row.names(hatten_trap) <- hatten_gene2tx$ensembl_transcript_id

common_tx <- intersect(row.names(scott_GNP_counts), row.names(hatten_trap))
hatten_trap      <- hatten_trap[match(common_tx, row.names(hatten_trap)), ]
scott_GNP_counts <- scott_GNP_counts[match(common_tx, row.names(scott_GNP_counts)), ]

RNAseq_list <- lapply(list(scott_GNP_counts, hatten_trap), as.matrix)
names(RNAseq_list) <- c("scott_GNP", "hatten_GNP")
meta_list <- mt_meta_list()

scott_gnp_pval <- as_tibble(scott_gnp_pval, rownames = "ensembl_transcript_id")
pval_all <- anno %>%
  dplyr::left_join(hatten_trap_pval, by = "ensembl_gene_id") %>%
  dplyr::left_join(scott_gnp_pval, by = "ensembl_transcript_id")

# ---- save -----------------------------------------------------------------
out <- file.path(MT_DATA, "normalized_merged.RData")
save(RNAseq_list, meta_list, pval_all, file = out)
message(sprintf("Wrote %s  (%d shared transcripts; Scott %d cols, Hatten %d cols)",
                out, nrow(RNAseq_list[[1]]), ncol(RNAseq_list[[1]]), ncol(RNAseq_list[[2]])))
