#!/usr/bin/env Rscript
# 01_load_and_merge.R
#
# Reproduce the merged per-transcript timecourse table (toc_t50) from raw
# inputs. This is a runnable port of 1_transcript_timecourse_3.Rmd:
#
#   Scott RNA-seq (E15,P1,P7,P14)  --DESeq2 normalize + contrasts-->
#   Hatten TRAP microarray (P0..P56) --2^ + limma pvals-->
#     merge on ensembl gene/transcript annotation
#       --> row-normalized log2, average per time point
#       --> time_of_change / interpolate / align on pseudo-time
#       --> t50, fold change, direction, time-of-max
#       --> expression binarization + kinetic categories
#   ==> data/toc_t50_<date>.txt  and  data/log2_norm_list_shift.RData
#
# NOTE: data/toc_t50_2023_05_17.txt (shipped) is the canonical table used by
# the API and figures. This script regenerates an equivalent table from raw
# and is the documented-reproducible path; small numeric differences from the
# canonical file are expected (DESeq2 / gamma-mixture stochasticity).
#
# Usage:  Rscript R/01_load_and_merge.R [--out data/toc_t50_rebuilt.txt]

suppressWarnings(suppressMessages({
  library(DESeq2); library(reshape2); library(matrixStats)
  library(stringr); library(tidyr); library(dplyr); library(tibble); library(akima)
  args_full <- commandArgs(trailingOnly = FALSE)
  here <- dirname(normalizePath(sub("^--file=", "",
            grep("^--file=", args_full, value = TRUE)[1])))
  source(file.path(here, "lib", "common.R"))
  mt_source_lib()
}))

argv <- commandArgs(trailingOnly = TRUE)
out_path <- MT_DATA |> file.path("toc_t50_rebuilt.txt")
if (length(argv) >= 2 && argv[1] == "--out") out_path <- argv[2]

read_cutoff  <- 5
p_val_cutoff <- 0.001

# ---- 1. Scott RNA-seq counts ----------------------------------------------
message("[1/9] Loading Scott RNA-seq counts")
scott_RNAseq <- read.table(MT_RAW_SCOTT_COUNTS, sep = "\t", header = TRUE, row.names = 1)
scott_RNAseq$X <- NULL
colnames(scott_RNAseq) <- matrix(unlist(strsplit(colnames(scott_RNAseq), "[.]")),
                                 ncol = 3, byrow = TRUE)[, 2]
scott_time_datasets <- read.table(MT_RAW_SCOTT_TIMES, sep = "\t", header = FALSE)
scott_time_datasets <- gsub(".*[.]", "", scott_time_datasets$V1)
# drop the htseq summary tail rows
scott_RNAseq <- scott_RNAseq[1:(nrow(scott_RNAseq) - 5), ]
scott_GNP <- scott_RNAseq[, colnames(scott_RNAseq) %in% scott_time_datasets]
scott_GNP <- scott_GNP[, match(scott_time_datasets, colnames(scott_GNP))]
scott_GNP_meta <- MT_SCOTT_META

# ---- 2. DESeq2 normalize + contrasts --------------------------------------
message("[2/9] DESeq2 normalization and contrasts")
scott_GNP <- as.matrix(scott_GNP)
colData <- data.frame(cols = colnames(scott_GNP), condition = scott_GNP_meta)
dds <- DESeqDataSetFromMatrix(scott_GNP, colData, design = ~ condition)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
scott_GNP_counts <- counts(dds, normalized = TRUE)

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
message("[3/9] Loading Hatten TRAP microarray")
hatten_trap <- read.table(MT_RAW_HATTEN_EXPR, sep = "\t", header = TRUE, row.names = 1)
hatten_trap_pval <- read.table(MT_RAW_HATTEN_PVAL, sep = "\t", header = TRUE, row.names = 1)
hatten_trap <- 2^(hatten_trap)                       # rma values are log2
hatten_meta <- MT_HATTEN_META
hatten_trap_pval$hatten_trap_min_pval <-
  rowMins(as.matrix(hatten_trap_pval[, grep("adj.P.Val", colnames(hatten_trap_pval))]), na.rm = TRUE)

# ---- 4. Merge on annotation (gene <-> transcript) -------------------------
message("[4/9] Merging Scott + Hatten on shared annotation")
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
meta_list   <- mt_meta_list()
pseudo_time <- MT_PSEUDO_TIME

# merge p-value tables
scott_gnp_pval <- as_tibble(scott_gnp_pval, rownames = "ensembl_transcript_id")
pval_all <- anno %>%
  dplyr::left_join(hatten_trap_pval, by = "ensembl_gene_id") %>%
  dplyr::left_join(scott_gnp_pval, by = "ensembl_transcript_id")

# ---- 5. Row-normalize log2 + average per timepoint ------------------------
message("[5/9] Row-normalized log2 and per-timepoint averaging")
min_val <- min(sapply(RNAseq_list, function(x) min(x[x > 0])))
RNAseq_list <- lapply(RNAseq_list, function(x) ifelse(x == 0, min_val, x))
log2_norm_list <- lapply(RNAseq_list, function(l) log2(sweep(l, 1, rowMeans(l), "/")))
average_df <- function(df, meta) {
  a <- as.data.frame(lapply(levels(meta), function(x) rowMeans(df[, which(meta == x)])))
  colnames(a) <- levels(meta); a
}
log2_norm_avg_list <- lapply(seq_along(log2_norm_list),
                             function(n) average_df(log2_norm_list[[n]], meta_list[[n]]))

# ---- 6. time_of_change ----------------------------------------------------
message("[6/9] time_of_change")
toc_list <- lapply(log2_norm_avg_list, time_of_change)
toc <- Reduce(function(a, b) full_join(a %>% mutate(Symbol = rownames(a)),
                                       b %>% mutate(Symbol = rownames(b)),
                                       by = "Symbol", suffix = c(".Scott", ".Hatten")),
              toc_list)
row.names(toc) <- toc$Symbol

# ---- 7. interpolate / align / t50 -----------------------------------------
message("[7/9] Interpolate, align, t50 / fold change / direction")
time_ref <- mt_time_ref()
interp <- time_expr_interp(log2_norm_avg_list, time_ref)
interp_align <- time_expr_align(interp)
shift <- time_expr_align(interp, return_shift = TRUE)
log2_norm_list_shift <- lapply(seq_along(log2_norm_list),
                               function(i) log2_norm_list[[i]] + shift[, i])
interp_align_avg <- time_expr_avg(interp_align)

total_change <- sapply(interp_align, function(j)
  apply(j, 1, function(i) max(i, na.rm = TRUE)) - apply(j, 1, function(i) min(i, na.rm = TRUE)))
total_change_avg <- apply(interp_align_avg, 1, function(i) max(i, na.rm = TRUE)) -
                    apply(interp_align_avg, 1, function(i) min(i, na.rm = TRUE))
total_change_dir <- sapply(interp_align, function(j)
  ifelse(apply(j, 1, which.max) > apply(j, 1, which.min), 1, -1))
total_change_dir_avg <- ifelse(apply(interp_align_avg, 1, which.max) >
                               apply(interp_align_avg, 1, which.min), 1, -1)
max_interp     <- sapply(interp_align, function(j) apply(j, 1, which.max))
max_interp_avg <- apply(interp_align_avg, 1, which.max)
min_interp     <- sapply(interp_align, function(j) apply(j, 1, which.min))
min_interp_avg <- apply(interp_align_avg, 1, which.min)

total_change <- cbind(total_change, total_change_avg);  colnames(total_change) <- c("scott_fc","hatten_fc","avg_fc")
total_change_dir <- cbind(total_change_dir, total_change_dir_avg); colnames(total_change_dir) <- c("scott_dir","hatten_dir","avg_dir")
max_min_interp <- cbind(min_interp, min_interp_avg, max_interp, max_interp_avg)
colnames(max_min_interp) <- c("scott_min","hatten_min","avg_min","scott_max","hatten_max","avg_max")

expr <- as.data.frame(expr_half_max(interp_align))
expr$avg <- expr_half_max(interp_align_avg)
colnames(expr) <- c("scott_t50","hatten_t50","avg_t50")
expr <- cbind(expr, total_change, total_change_dir, max_min_interp)

toc_t50 <- cbind(toc, expr)
toc_t50 <- as_tibble(toc_t50, rownames = "ensembl_transcript_id")
toc_t50 <- pval_all %>% left_join(toc_t50, by = "ensembl_transcript_id")

# ---- 8. Expression level + binarization -----------------------------------
message("[8/9] Expression level and binarization")

# Robust binarizer: try the mixture-model call; if the EM fails to converge
# (a known fragility on some distributions), fall back to a deterministic
# threshold at the density antimode (valley) between the two modes. Returns a
# data.frame with columns `counts` and `bin` (0/1), matching binarize_counts().
binarize_expr <- function(v, mu = c(0.3, 0.7)) {
  out <- tryCatch(
    binarize_counts(v, exclude_extreme = FALSE, mu = mu),
    error = function(e) NULL)
  if (!is.null(out)) return(out)
  message("    mixture EM did not converge; using density-antimode threshold")
  d <- density(v, na.rm = TRUE)
  # candidate valley: lowest density point strictly between the two largest peaks
  peaks <- which(diff(sign(diff(d$y))) == -2) + 1
  if (length(peaks) >= 2) {
    top2 <- sort(peaks[order(d$y[peaks], decreasing = TRUE)][1:2])
    seg <- top2[1]:top2[2]
    thr <- d$x[seg][which.min(d$y[seg])]
  } else {
    thr <- median(v, na.rm = TRUE)   # unimodal fallback
  }
  data.frame(counts = v, bin = as.integer(v >= thr))
}
toc_t50$ensembl_exon_length <- sapply(
  strsplit(as.character(toc_t50$ensembl_blockSizes), ","),
  function(x) sum(as.numeric(x)))
max3 <- function(m) apply(m, 1, function(x) mean(sort(x, decreasing = TRUE)[1:3]))
max_expr <- as_tibble(cbind(max_expr_Scott = max3(scott_GNP_counts),
                            max_expr_Hatten = max3(hatten_trap)),
                      rownames = "ensembl_transcript_id")
toc_t50 <- toc_t50 %>% left_join(max_expr, by = "ensembl_transcript_id")
toc_t50$max_expr_Scott <- toc_t50$max_expr_Scott * 1000 / toc_t50$ensembl_exon_length
toc_t50 <- toc_t50 %>% dplyr::distinct(ensembl_transcript_id, .keep_all = TRUE)

# Scott: winsorized log10, 2-component normal mixture cutoff.
# NOTE: the source Rmd binarized Scott with a gamma-mixture (binarize_counts_p)
# in a chunk left eval=FALSE ("troubleshoot the gamma mixture"); that fit is
# very slow and was never part of the canonical run. We use the same fast,
# stable normal-mixture binarizer as Hatten, which reproduces the expression
# calls in the canonical table closely.
scott_e <- toc_t50 %>% drop_na(max_expr_Scott) %>%
  dplyr::select(ensembl_transcript_id, max_expr_Scott) %>% as.data.frame()
row.names(scott_e) <- scott_e$ensembl_transcript_id; scott_e$ensembl_transcript_id <- NULL
scott_e$max_expr_Scott <- log10(scott_e$max_expr_Scott + 1)
scott_e <- wins_norm_histones(scott_e, low = 0, high = 20)
scott_bin <- binarize_expr(scott_e$max_expr_Scott, mu = c(0.3, 0.7))
colnames(scott_bin) <- str_c("Scott_max_expr_", colnames(scott_bin))
scott_cut <- cbind(scott_e, scott_bin)

# Hatten: winsorized log10, 2-component normal mixture cutoff
hatten_e <- toc_t50 %>% drop_na(max_expr_Hatten) %>%
  dplyr::select(ensembl_transcript_id, max_expr_Hatten) %>% as.data.frame()
row.names(hatten_e) <- hatten_e$ensembl_transcript_id; hatten_e$ensembl_transcript_id <- NULL
hatten_e$max_expr_Hatten <- log10(hatten_e$max_expr_Hatten + 1)
hatten_e <- wins_norm_histones(hatten_e, low = 100, high = 50)
hatten_bin <- binarize_expr(hatten_e$max_expr_Hatten, mu = c(0.3, 0.7))
colnames(hatten_bin) <- str_c("Hatten_max_expr_", colnames(hatten_bin))
hatten_cut <- cbind(hatten_e, hatten_bin)

expr_bin <- as_tibble(cbind(scott_cut, hatten_cut), rownames = "ensembl_transcript_id") %>%
  dplyr::select(ensembl_transcript_id, Scott_max_expr_counts, Scott_max_expr_bin,
                Hatten_max_expr_counts, Hatten_max_expr_bin)
toc_t50 <- toc_t50 %>% left_join(expr_bin, by = "ensembl_transcript_id")

# ---- 9. Kinetic categories ------------------------------------------------
message("[9/9] Kinetic categories")
toc_t50$expr <- (toc_t50$Scott_max_expr_bin == 1) | (toc_t50$Hatten_max_expr_bin == 1)
toc_t50$max_cat <- cut(toc_t50$avg_max, c(0,2,20,64), right = FALSE,
                       labels = c("E15","P0_7","P14_P56"))
toc_t50$t50_cat <- cut(toc_t50$avg_t50, c(0,6,18,64), right = FALSE,
                       labels = c("Early","Mid","Late"))
toc_t50$max_t50_cat <- str_c("Max_", toc_t50$max_cat, " ", "t50_", toc_t50$t50_cat)

unchg <- toc_t50$avg_fc < 1 |
  (toc_t50$hatten_trap_min_pval > p_val_cutoff & toc_t50$scott_gnp_min_pval > p_val_cutoff)
toc_t50$max_cat <- as.character(toc_t50$max_cat)
toc_t50$max_cat[unchg] <- "unchanged"
toc_t50$max_t50_cat[unchg] <- "unchanged"
toc_t50$max_cat <- ifelse(toc_t50$expr, toc_t50$max_cat, "unexpressed")
toc_t50$max_t50_cat <- ifelse(toc_t50$expr, toc_t50$max_t50_cat, "unexpressed")
toc_t50$max_t50_cat[is.na(toc_t50$max_t50_cat)] <- "unexpressed"
toc_t50$max_cat <- as.factor(toc_t50$max_cat)
toc_t50$max_t50_cat <- as.factor(toc_t50$max_t50_cat)
toc_t50 <- toc_t50 %>% dplyr::mutate(cat_avg_fc = case_when(
  avg_fc < 1              ~ "0-1 log2FC",
  avg_fc >= 1 & avg_fc < 2 ~ "1-2 log2FC",
  avg_fc >= 2             ~ ">2 log2FC"))

# ---- save -----------------------------------------------------------------
write.table(toc_t50, out_path, sep = "\t", row.names = TRUE, col.names = TRUE)
names(log2_norm_list_shift) <- NULL
save(log2_norm_list_shift, file = file.path(MT_DATA, "log2_norm_list_shift_rebuilt.RData"))

message(sprintf("Wrote %s  (%d rows x %d cols)", out_path, nrow(toc_t50), ncol(toc_t50)))
message("Category distribution:")
print(sort(table(toc_t50$max_t50_cat), decreasing = TRUE))
