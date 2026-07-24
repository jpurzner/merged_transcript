# mt_plot.R -- timecourse plotting wrappers built on plot_group_timeline()
#
# Requires that common.R and the R/lib helpers have been sourced
# (mt_source_lib()).

# Build a gene2cluster vector (named: transcript_id -> label) for a set of
# gene symbols, using the query index to resolve symbols to transcript ids.
# Each gene becomes its own facet, labelled with the symbol (+ kinetic class).
mt_gene2cluster <- function(genes, idx = NULL, all_transcripts = FALSE,
                            label_with_class = TRUE) {
  if (is.null(idx)) idx <- mt_load_index()
  hits <- idx[idx$gene %in% genes, , drop = FALSE]
  if (nrow(hits) == 0) return(setNames(character(0), character(0)))
  if (!all_transcripts) {
    # keep the representative (largest |log2fc|) transcript per gene
    hits <- hits[order(hits$gene, -abs(hits$log2fc)), ]
    hits <- hits[!duplicated(hits$gene), ]
  }
  label <- hits$gene
  if (label_with_class) {
    label <- ifelse(is.na(hits$kinetic_class), hits$gene,
                    paste0(hits$gene, "\n", hits$kinetic_class))
  }
  setNames(label, hits$transcript_id)
}

# Load the slim query index (cached in an option for the session).
mt_load_index <- function(path = MT_QUERY_INDEX) {
  cached <- getOption("mt.index")
  if (!is.null(cached)) return(cached)
  idx <- read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                    quote = "", comment.char = "")
  options(mt.index = idx)
  idx
}

# Plot developmental timelines for one or more genes. Returns a ggplot.
#   genes           : character vector of gene symbols
#   expr            : list from mt_load_expr(); loaded if NULL
#   ylimits         : y-axis (log2 fold change) range
#   add_t50         : draw the median t50 marker on each facet
mt_plot_genes <- function(genes, expr = NULL, idx = NULL,
                          ylimits = c(-3, 3), all_transcripts = FALSE,
                          add_t50 = TRUE, col_num = NULL, text_size = 14) {
  if (is.null(expr)) expr <- mt_load_expr()
  if (is.null(idx))  idx  <- mt_load_index()

  g2c <- mt_gene2cluster(genes, idx = idx, all_transcripts = all_transcripts)
  if (length(g2c) == 0)
    stop("None of the requested genes were found in the query index: ",
         paste(genes, collapse = ", "))

  found <- unique(idx$gene[idx$transcript_id %in% names(g2c)])
  missing <- setdiff(genes, found)
  if (length(missing))
    message("Not found (skipped): ", paste(missing, collapse = ", "))

  if (is.null(col_num)) col_num <- min(4, length(unique(g2c)))

  # Optional t50 marker table: data.frame(transcript_id, avg_t50)
  t50_tbl <- NULL
  if (add_t50) {
    sub <- idx[idx$transcript_id %in% names(g2c), c("transcript_id", "t50")]
    colnames(sub) <- c("gene", "avg_t50")   # plot_group_timeline expects these
    t50_tbl <- sub
  }

  plot_group_timeline(
    expr_list    = expr,
    meta_list    = mt_meta_list(),
    gene2cluster = g2c,
    p_time       = MT_PSEUDO_TIME,
    ylimits      = ylimits,
    include_n    = FALSE,
    col_num      = col_num,
    text_size    = text_size,
    x_label      = "remove",
    add_t50      = t50_tbl,
    time_ref     = mt_time_ref()
  )
}

# Save a gene timeline plot to a file (pdf or png). Returns the path.
mt_save_gene_plot <- function(genes, file = NULL, format = c("pdf", "png"),
                              width = NULL, height = NULL, ...) {
  format <- match.arg(format)
  if (is.null(file)) {
    tag <- if (length(genes) <= 3) paste(genes, collapse = "_") else
      paste0(length(genes), "genes")
    file <- file.path(MT_FIGURES, paste0("timeline_", tag, ".", format))
  }
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  p <- mt_plot_genes(genes, ...)
  n <- length(unique(mt_gene2cluster(genes)))
  cols <- min(4, n); rows <- ceiling(n / cols)
  if (is.null(width))  width  <- 2.6 * cols + 1
  if (is.null(height)) height <- 2.6 * rows + 1
  ggplot2::ggsave(file, p, width = width, height = height, limitsize = FALSE,
                  bg = "white")
  file
}
