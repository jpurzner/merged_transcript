# mt_values.R -- extract the numerical timecourse behind the plots.
#
# The per-time-point aligned log2 fold-change values (mean +/- spread across
# replicates) for a gene, from the aligned expression matrices
# (log2_norm_list_shift) -- i.e. the numbers a plot draws as points/ribbon.
# Requires common.R to have been sourced.

# Resolve a gene symbol to the transcript used for its timecourse: the
# representative (largest |log2fc|) transcript that is present in the aligned
# matrices. Returns NA if the gene has no transcript in the aligned set.
mt_resolve_aligned_tx <- function(gene, idx, aligned_ids) {
  hit <- idx[idx$gene == gene | tolower(idx$gene) == tolower(gene) |
             idx$transcript_id == gene, , drop = FALSE]
  hit <- hit[hit$transcript_id %in% aligned_ids, , drop = FALSE]
  if (nrow(hit) == 0) return(NA_character_)
  hit <- hit[order(-abs(hit$log2fc)), ]
  hit$transcript_id[1]
}

# Long-format numerical timecourse for one or more genes.
# Columns: gene, transcript_id, dataset, time, pseudo_time_index, n,
#          mean, sd, sem  (+ list-column `replicates` if with_replicates=TRUE)
mt_timecourse_values <- function(genes, expr = NULL, idx = NULL,
                                 with_replicates = FALSE) {
  if (is.null(expr)) expr <- mt_load_expr()
  if (is.null(idx))  idx  <- mt_load_index()
  metas       <- mt_meta_list()
  aligned_ids <- rownames(expr[[1]])
  pt_index    <- setNames(seq_along(MT_PSEUDO_TIME), MT_PSEUDO_TIME)

  rows <- list()
  attr_out <- list(found = character(0), missing = character(0),
                   not_aligned = character(0))
  for (g in genes) {
    tid <- mt_resolve_aligned_tx(g, idx, aligned_ids)
    if (is.na(tid)) {
      # distinguish "unknown gene" from "known but no aligned transcript"
      if (any(idx$gene == g | tolower(idx$gene) == tolower(g)))
        attr_out$not_aligned <- c(attr_out$not_aligned, g)
      else
        attr_out$missing <- c(attr_out$missing, g)
      next
    }
    attr_out$found <- c(attr_out$found, g)
    for (d in seq_along(expr)) {
      v <- expr[[d]][tid, ]
      m <- metas[[d]]
      ds <- names(metas)[d]
      for (lv in levels(m)) {
        vals <- as.numeric(v[m == lv])
        n <- length(vals)
        row <- data.frame(
          gene = g, transcript_id = tid, dataset = ds, time = lv,
          pseudo_time_index = unname(pt_index[lv]), n = n,
          mean = mean(vals), sd = stats::sd(vals),
          sem = stats::sd(vals) / sqrt(n),
          stringsAsFactors = FALSE)
        if (with_replicates) row$replicates <- I(list(vals))
        rows[[length(rows) + 1]] <- row
      }
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else
    data.frame(gene = character(0), transcript_id = character(0),
               dataset = character(0), time = character(0),
               pseudo_time_index = integer(0), n = integer(0),
               mean = numeric(0), sd = numeric(0), sem = numeric(0))
  if (nrow(out)) out <- out[order(out$gene, out$pseudo_time_index, out$dataset), ]
  attr(out, "resolution") <- attr_out
  out
}
