#!/usr/bin/env Rscript
# 02_align_timecourse.R   —   STAGE 2 of the reproducible pipeline
#
# Turn the normalized, merged expression (from 01_normalize_and_merge.R) into
# the ALIGNMENT DATASETS and the per-transcript kinetic table.
#
#   normalized counts / intensity
#     --row-normalize to log2 fold change vs each gene's own mean-->  log2_norm_list
#     --average replicates within each time point-->                  log2_norm_avg_list
#     --interpolate onto the E15..P56 day axis (approx)-->            *_interp
#     --align: shift the 2nd series so overlapping days agree,
#              then centre each gene on zero-->                       *_interp_align
#     --carry the per-gene shift back onto the per-replicate data-->  log2_norm_list_shift
#     --t50 / fold change / direction / time-of-max / expression-->   toc_t50
#
# OUTPUTS
#   data/alignment_datasets.RData   the alignment datasets (documented bundle,
#                                   incl. log2_norm_list_shift used for plotting)
#   data/log2_norm_list_shift.RData the per-replicate aligned matrices (kept as
#                                   a standalone file for the plotting layer)
#   data/toc_t50_rebuilt.txt        the rebuilt per-transcript kinetic table
#
# Usage:  Rscript R/02_align_timecourse.R [--out data/toc_t50_rebuilt.txt]

suppressWarnings(suppressMessages({
  library(reshape2); library(matrixStats)
  library(stringr); library(tidyr); library(dplyr); library(tibble); library(akima)
  args_full <- commandArgs(trailingOnly = FALSE)
  here <- dirname(normalizePath(sub("^--file=", "",
            grep("^--file=", args_full, value = TRUE)[1])))
  source(file.path(here, "lib", "common.R"))
  mt_source_lib()
}))

argv <- commandArgs(trailingOnly = TRUE)
out_path <- file.path(MT_DATA, "toc_t50_rebuilt.txt")
if (length(argv) >= 2 && argv[1] == "--out") out_path <- argv[2]

p_val_cutoff <- 0.001

in_file <- file.path(MT_DATA, "normalized_merged.RData")
if (!file.exists(in_file))
  stop("Missing ", in_file, ". Run:  Rscript R/01_normalize_and_merge.R  first.")
load(in_file)                        # RNAseq_list, meta_list, pval_all
scott_GNP_counts <- RNAseq_list[[1]]
hatten_trap      <- RNAseq_list[[2]]

# ---- 5. Row-normalize to log2 + average per time point --------------------
message("[5/9] Row-normalized log2 fold change and per-timepoint averaging")
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

# ---- 7. Interpolate onto the day axis, align, shift -----------------------
message("[7/9] Interpolate, align, shift  (the alignment datasets)")
time_ref <- mt_time_ref()
log2_norm_avg_list_interp       <- time_expr_interp(log2_norm_avg_list, time_ref)
log2_norm_avg_list_interp_align <- time_expr_align(log2_norm_avg_list_interp)
shift <- time_expr_align(log2_norm_avg_list_interp, return_shift = TRUE)
log2_norm_list_shift <- lapply(seq_along(log2_norm_list),
                               function(i) log2_norm_list[[i]] + shift[, i])
names(log2_norm_list_shift) <- NULL
log2_norm_avg_list_interp_align_avg <- time_expr_avg(log2_norm_avg_list_interp_align)

# kinetics derived from the aligned data
ia  <- log2_norm_avg_list_interp_align
iaa <- log2_norm_avg_list_interp_align_avg
total_change <- sapply(ia, function(j)
  apply(j, 1, function(i) max(i, na.rm = TRUE)) - apply(j, 1, function(i) min(i, na.rm = TRUE)))
total_change_avg <- apply(iaa, 1, function(i) max(i, na.rm = TRUE)) -
                    apply(iaa, 1, function(i) min(i, na.rm = TRUE))
total_change_dir <- sapply(ia, function(j)
  ifelse(apply(j, 1, which.max) > apply(j, 1, which.min), 1, -1))
total_change_dir_avg <- ifelse(apply(iaa, 1, which.max) > apply(iaa, 1, which.min), 1, -1)
max_interp <- sapply(ia, function(j) apply(j, 1, which.max)); max_interp_avg <- apply(iaa, 1, which.max)
min_interp <- sapply(ia, function(j) apply(j, 1, which.min)); min_interp_avg <- apply(iaa, 1, which.min)

total_change <- cbind(total_change, total_change_avg); colnames(total_change) <- c("scott_fc","hatten_fc","avg_fc")
total_change_dir <- cbind(total_change_dir, total_change_dir_avg); colnames(total_change_dir) <- c("scott_dir","hatten_dir","avg_dir")
max_min_interp <- cbind(min_interp, min_interp_avg, max_interp, max_interp_avg)
colnames(max_min_interp) <- c("scott_min","hatten_min","avg_min","scott_max","hatten_max","avg_max")

expr <- as.data.frame(expr_half_max(ia)); expr$avg <- expr_half_max(iaa)
colnames(expr) <- c("scott_t50","hatten_t50","avg_t50")
expr <- cbind(expr, total_change, total_change_dir, max_min_interp)

toc_t50 <- cbind(toc, expr)
toc_t50 <- as_tibble(toc_t50, rownames = "ensembl_transcript_id")
toc_t50 <- pval_all %>% left_join(toc_t50, by = "ensembl_transcript_id")

# ---- save the alignment datasets ------------------------------------------
alignment_datasets <- list(
  log2_norm_list        = log2_norm_list,                     # per-replicate log2FC
  log2_norm_avg_list    = log2_norm_avg_list,                 # per-timepoint mean
  interp                = log2_norm_avg_list_interp,          # interpolated onto day axis
  interp_align          = log2_norm_avg_list_interp_align,    # aligned + zero-centred
  interp_align_avg      = log2_norm_avg_list_interp_align_avg,# cross-dataset mean
  shift                 = shift,                              # per-gene shift applied
  log2_norm_list_shift  = log2_norm_list_shift,               # per-replicate aligned (plots)
  meta_list             = meta_list,
  pseudo_time           = MT_PSEUDO_TIME,
  time_ref              = time_ref
)
save(alignment_datasets, file = file.path(MT_DATA, "alignment_datasets.RData"))
save(log2_norm_list_shift, file = file.path(MT_DATA, "log2_norm_list_shift.RData"))
message("    wrote data/alignment_datasets.RData and data/log2_norm_list_shift.RData")

# ---- 8. Expression level + binarization -----------------------------------
message("[8/9] Expression level and binarization")

# Robust binarizer: try the mixture-model call; if the EM fails to converge
# (a known fragility on some distributions), fall back to a deterministic
# threshold at the density antimode (valley) between the two modes. Returns a
# data.frame with columns `counts` and `bin` (0/1), matching binarize_counts().
binarize_expr <- function(v, mu = c(0.3, 0.7)) {
  out <- tryCatch(binarize_counts(v, exclude_extreme = FALSE, mu = mu),
                  error = function(e) NULL)
  if (!is.null(out)) return(out)
  message("    mixture EM did not converge; using density-antimode threshold")
  d <- density(v, na.rm = TRUE)
  peaks <- which(diff(sign(diff(d$y))) == -2) + 1
  if (length(peaks) >= 2) {
    top2 <- sort(peaks[order(d$y[peaks], decreasing = TRUE)][1:2])
    seg <- top2[1]:top2[2]; thr <- d$x[seg][which.min(d$y[seg])]
  } else thr <- median(v, na.rm = TRUE)
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

# Scott: winsorized log10, normal-mixture cutoff (see 01 note on the Rmd's
# eval=FALSE gamma-mixture step -- we use the stable normal-mixture / fallback).
scott_e <- toc_t50 %>% drop_na(max_expr_Scott) %>%
  dplyr::select(ensembl_transcript_id, max_expr_Scott) %>% as.data.frame()
row.names(scott_e) <- scott_e$ensembl_transcript_id; scott_e$ensembl_transcript_id <- NULL
scott_e$max_expr_Scott <- log10(scott_e$max_expr_Scott + 1)
scott_e <- wins_norm_histones(scott_e, low = 0, high = 20)
scott_bin <- binarize_expr(scott_e$max_expr_Scott, mu = c(0.3, 0.7))
colnames(scott_bin) <- str_c("Scott_max_expr_", colnames(scott_bin))
scott_cut <- cbind(scott_e, scott_bin)

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
toc_t50$max_cat <- cut(toc_t50$avg_max, c(0,2,20,64), right = FALSE, labels = c("E15","P0_7","P14_P56"))
toc_t50$t50_cat <- cut(toc_t50$avg_t50, c(0,6,18,64), right = FALSE, labels = c("Early","Mid","Late"))
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
  avg_fc < 1               ~ "0-1 log2FC",
  avg_fc >= 1 & avg_fc < 2 ~ "1-2 log2FC",
  avg_fc >= 2              ~ ">2 log2FC"))

# ---- save -----------------------------------------------------------------
write.table(toc_t50, out_path, sep = "\t", row.names = TRUE, col.names = TRUE)
message(sprintf("Wrote %s  (%d rows x %d cols)", out_path, nrow(toc_t50), ncol(toc_t50)))
print(sort(table(toc_t50$max_t50_cat), decreasing = TRUE))
