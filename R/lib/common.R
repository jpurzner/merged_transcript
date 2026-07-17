# common.R -- shared constants, paths, and loaders for the merged_transcript project
#
# This project merges two granule-neuron-precursor (GNP) developmental
# transcriptomic time courses of the mouse cerebellum:
#
#   * Scott RNA-seq   : E15, P1, P7, P14   (11 samples)
#   * Hatten TRAP array: P0, P7, P12, P18, P21, P56 (18 samples)
#
# The two series are interpolated onto a common postnatal-day pseudo-time
# axis, aligned, and reduced to per-transcript kinetic metrics (t50, fold
# change, direction, time-of-max category). See docs/SKILL.md for the full
# description of every derived column.

# ---- repo layout -----------------------------------------------------------

# Locate this file (R/lib/common.R). When sourced, sys.frame's `ofile` holds
# its own path; that is preferred over commandArgs, which points at the entry
# script (R/<something>.R), not this file.
.mt_this_file <- function() {
  for (i in seq_len(sys.nframe())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(normalizePath(of))
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(normalizePath(sub("^--file=", "", file_arg[1])))
  NA_character_
}

# common.R lives at <root>/R/lib/common.R, so the root is two directories up.
MT_ROOT <- local({
  env <- Sys.getenv("MERGED_TRANSCRIPT_ROOT")
  if (nzchar(env)) return(normalizePath(env))
  f <- .mt_this_file()
  if (!is.na(f)) return(normalizePath(file.path(dirname(f), "..", "..")))
  # last resort: current working directory
  normalizePath(getwd())
})

MT_DATA    <- file.path(MT_ROOT, "data")
MT_LIB     <- file.path(MT_ROOT, "R", "lib")
MT_FIGURES <- file.path(MT_ROOT, "figures")

# Shipped, precomputed artifacts (source of truth for the API + plots)
MT_TOC_T50      <- file.path(MT_DATA, "toc_t50_2023_05_17.txt")
MT_QUERY_INDEX  <- file.path(MT_DATA, "timecourse_query_index.tsv")
MT_EXPR_RDATA   <- file.path(MT_DATA, "log2_norm_list_shift.RData")

# Alignment datasets produced by 02_align_timecourse.R
MT_ALIGN_RDATA    <- file.path(MT_DATA, "alignment_datasets.RData")
MT_NORMMERGE_RDATA<- file.path(MT_DATA, "normalized_merged.RData")

# Raw inputs for the reproducible merge (01_load_and_merge.R). These live in
# data/raw/ as symlinks to their original locations and are git-ignored; see
# data/raw/README.md. Any path can be overridden with the matching env var.
MT_RAW <- file.path(MT_DATA, "raw")
mt_raw <- function(name, env) {
  e <- Sys.getenv(env)
  if (nzchar(e)) return(e)
  file.path(MT_RAW, name)
}
MT_RAW_SCOTT_COUNTS <- mt_raw("scott_rnaseq_counts.txt",        "MT_SCOTT_COUNTS")
MT_RAW_SCOTT_TIMES  <- mt_raw("scott_time_series_name.txt",     "MT_SCOTT_TIMES")
MT_RAW_HATTEN_EXPR  <- mt_raw("hatten_trap_microarray.txt",     "MT_HATTEN_EXPR")
MT_RAW_HATTEN_PVAL  <- mt_raw("hatten_trap_microarray_pval.txt","MT_HATTEN_PVAL")
MT_RAW_ANNO         <- mt_raw("gnp_tss_model_anno.txt",         "MT_ANNO")

# ---- experimental design constants (from 1_transcript_timecourse_3.Rmd) ----

# Per-replicate time labels for each dataset. Order matches the columns of
# the matrices in log2_norm_list_shift.
MT_SCOTT_META  <- factor(
  c("E15", "E15", "E15", "P1", "P1", "P7", "P7", "P7", "P14", "P14", "P14"),
  levels = c("E15", "P1", "P7", "P14")
)
MT_HATTEN_META <- factor(
  c("P0","P0","P0","P7","P7","P7","P12","P12","P12",
    "P18","P18","P18","P21","P21","P21","P56","P56","P56"),
  levels = c("P0", "P7", "P12", "P18", "P21", "P56")
)

# Ordered union of all developmental time points (the pseudo-time axis).
MT_PSEUDO_TIME <- c("E15", "P0", "P1", "P7", "P12", "P14", "P18", "P21", "P56")

# meta_list as consumed by plot_group_timeline()
mt_meta_list <- function() {
  ml <- list(MT_SCOTT_META, MT_HATTEN_META)
  names(ml) <- c("scott_GNP", "hatten_GNP")
  ml
}

# time_ref: maps each day on the E15..P56 axis to which dataset sampled it.
mt_time_ref <- function() {
  tr <- data.frame(time = c(paste0("E", 15:20), paste0("P", 0:56)),
                   stringsAsFactors = FALSE)
  tr$scott  <- tr$time %in% c("E15", "P1", "P7", "P14")
  tr$hatten <- tr$time %in% c("P0", "P7", "P12", "P18", "P21", "P56")
  rownames(tr) <- tr$time
  tr$time <- NULL
  tr
}

# ---- loaders ---------------------------------------------------------------

# Source every helper function in R/lib (except this file).
mt_source_lib <- function(quiet = TRUE) {
  files <- list.files(MT_LIB, pattern = "\\.R$", full.names = TRUE)
  files <- files[basename(files) != "common.R"]
  for (f in files) {
    if (!quiet) message("source: ", basename(f))
    suppressWarnings(source(f))
  }
  invisible(files)
}

# Load the merged per-transcript table (toc_t50). Cached in an option so
# repeated calls in one session are cheap.
mt_load_toc <- function(path = MT_TOC_T50) {
  read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
             quote = "\"", comment.char = "")
}

# Load the aligned, row-normalized log2 expression matrices used for plotting.
# Returns list(scott = <matrix 11 cols>, hatten = <matrix 18 cols>).
# Prefers the standalone log2_norm_list_shift.RData; falls back to pulling it
# out of the full alignment_datasets.RData bundle if that is all that's present.
mt_load_expr <- function(path = MT_EXPR_RDATA) {
  e <- new.env()
  if (file.exists(path)) {
    load(path, envir = e); L <- e$log2_norm_list_shift
  } else if (file.exists(MT_ALIGN_RDATA)) {
    load(MT_ALIGN_RDATA, envir = e); L <- e$alignment_datasets$log2_norm_list_shift
  } else {
    stop("No aligned expression found. Run:  Rscript R/run_pipeline.R --step merge")
  }
  names(L) <- c("scott_GNP", "hatten_GNP")
  L
}

# Load the full alignment-datasets bundle produced by 02_align_timecourse.R.
mt_load_alignment <- function(path = MT_ALIGN_RDATA) {
  if (!file.exists(path))
    stop("Missing ", path, ". Run:  Rscript R/run_pipeline.R --step merge")
  e <- new.env(); load(path, envir = e); e$alignment_datasets
}
