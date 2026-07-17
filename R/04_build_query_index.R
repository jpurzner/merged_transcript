#!/usr/bin/env Rscript
# 04_build_query_index.R
#
# Reduce the 115-column merged table (toc_t50) to a slim, query-optimized
# index: one row per transcript, gene symbol as the lookup key, only the
# interpretable kinetic / expression / significance columns that an agent
# would ask about. Writes data/timecourse_query_index.tsv.
#
# Usage:  Rscript R/03_build_query_index.R

suppressWarnings(suppressMessages({
  args <- commandArgs(trailingOnly = FALSE)
  here <- dirname(normalizePath(sub("^--file=", "",
            grep("^--file=", args, value = TRUE)[1])))
  source(file.path(here, "lib", "common.R"))
}))

message("Loading merged table: ", MT_TOC_T50)
d <- mt_load_toc()

# Columns to carry into the slim index, with clearer names where useful.
keep <- c(
  gene              = "mgi_symbol",
  transcript_id     = "ensembl_transcript_id",
  gene_id           = "ensembl_gene_id",
  chr               = "ensembl_chr",
  start             = "ensembl_chromStart",
  end               = "ensembl_chromEnd",
  strand            = "ensembl_strand",
  # kinetics ---------------------------------------------------------------
  t50               = "avg_t50",       # time to 50% of max change (pseudo-time index)
  t50_scott         = "scott_t50",
  t50_hatten        = "hatten_t50",
  log2fc            = "avg_fc",         # magnitude of change (log2), unsigned
  log2fc_scott      = "scott_fc",
  log2fc_hatten     = "hatten_fc",
  direction         = "avg_dir",        # +1 increasing, -1 decreasing over development
  time_of_max       = "avg_max",        # pseudo-time index of peak expression
  time_of_min       = "avg_min",
  # categorical calls ------------------------------------------------------
  max_category      = "max_cat",        # E15 / P0_7 / P14_P56 / unchanged / unexpressed
  t50_category      = "t50_cat",        # Early / Mid / Late
  kinetic_class     = "max_t50_cat",    # combined Max_* t50_* label
  fc_class          = "cat_avg_fc",     # 0-1 / 1-2 / >2 log2FC
  expressed         = "expr",
  # expression level -------------------------------------------------------
  max_expr_scott    = "max_expr_Scott",
  max_expr_hatten   = "max_expr_Hatten",
  # significance -----------------------------------------------------------
  min_padj_scott    = "scott_gnp_min_pval",
  min_padj_hatten   = "hatten_trap_min_pval"
)

missing <- setdiff(unname(keep), colnames(d))
if (length(missing)) stop("toc_t50 is missing expected columns: ",
                          paste(missing, collapse = ", "))

idx <- d[, unname(keep)]
colnames(idx) <- names(keep)

# Human-readable direction, and a signed fold change for convenience.
idx$direction_label <- ifelse(idx$direction > 0, "up",
                       ifelse(idx$direction < 0, "down", "flat"))
idx$signed_log2fc   <- idx$log2fc * idx$direction

# Human-readable developmental-day labels for the day-axis indices (t50,
# time_of_max, time_of_min live on the 63-point E15..P56 daily axis).
idx$t50_day         <- mt_day_label(idx$t50)
idx$time_of_max_day <- mt_day_label(idx$time_of_max)
idx$time_of_min_day <- mt_day_label(idx$time_of_min)

# Round the floating-point metrics for a compact, readable file.
num_cols <- c("t50","t50_scott","t50_hatten","log2fc","log2fc_scott",
              "log2fc_hatten","time_of_max","time_of_min","max_expr_scott",
              "max_expr_hatten","min_padj_scott","min_padj_hatten","signed_log2fc")
for (cc in num_cols) idx[[cc]] <- round(idx[[cc]], 4)

# Drop rows with no usable gene symbol (can't be queried by name).
idx <- idx[!is.na(idx$gene) & idx$gene != "", ]

# Order by gene symbol then descending magnitude so the "representative"
# (largest-change) transcript for a gene comes first.
idx <- idx[order(idx$gene, -abs(idx$log2fc)), ]

write.table(idx, MT_QUERY_INDEX, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

message(sprintf("Wrote %s  (%d transcripts, %d unique genes, %d columns)",
                MT_QUERY_INDEX, nrow(idx), length(unique(idx$gene)), ncol(idx)))
