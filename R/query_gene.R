#!/usr/bin/env Rscript
# query_gene.R -- command-line API over the merged GNP transcript timecourse.
#
# Other agents (or humans) query a gene's developmental kinetics without
# touching R internals: pass one or more gene symbols, get structured output.
#
# USAGE
#   Rscript R/query_gene.R <GENE> [<GENE> ...] [options]
#
# OPTIONS
#   --format json|table|tsv   output format (default: json)
#   --fields a,b,c            restrict to these columns
#   --all-transcripts         return every transcript, not just the
#                             representative (largest |log2fc|) one per gene
#   --plot                    also render a timeline plot to figures/
#   --plot-format pdf|png     plot file type (default: pdf)
#   --index <path>            override the query index path
#   --list-fields             print available fields and exit
#   --help                    show this help
#
# EXAMPLES
#   Rscript R/query_gene.R Ezh2
#   Rscript R/query_gene.R Ezh2 Gabra6 Sox2 --format table
#   Rscript R/query_gene.R Atoh1 --fields gene,t50,log2fc,direction_label,kinetic_class
#   Rscript R/query_gene.R Ezh2 Neurod1 --plot --plot-format png
#
# EXIT STATUS
#   0  at least one gene resolved
#   2  no gene resolved / usage error

suppressWarnings(suppressMessages({
  .argv_full <- commandArgs(trailingOnly = FALSE)
  .here <- dirname(normalizePath(sub("^--file=", "",
             grep("^--file=", .argv_full, value = TRUE)[1])))
  source(file.path(.here, "lib", "common.R"))
}))

# ---- tiny arg parser -------------------------------------------------------
argv <- commandArgs(trailingOnly = TRUE)
opts <- list(format = "json", fields = NULL, all_transcripts = FALSE,
             plot = FALSE, plot_format = "pdf", index = MT_QUERY_INDEX,
             list_fields = FALSE, help = FALSE)
genes <- character(0)
i <- 1
while (i <= length(argv)) {
  a <- argv[i]
  if (a == "--help" || a == "-h") { opts$help <- TRUE
  } else if (a == "--list-fields") { opts$list_fields <- TRUE
  } else if (a == "--all-transcripts") { opts$all_transcripts <- TRUE
  } else if (a == "--plot") { opts$plot <- TRUE
  } else if (a == "--format") { i <- i + 1; opts$format <- argv[i]
  } else if (a == "--fields") { i <- i + 1; opts$fields <- strsplit(argv[i], ",")[[1]]
  } else if (a == "--plot-format") { i <- i + 1; opts$plot_format <- argv[i]
  } else if (a == "--index") { i <- i + 1; opts$index <- argv[i]
  } else if (grepl("^--", a)) { stop("Unknown option: ", a)
  } else { genes <- c(genes, a) }
  i <- i + 1
}

help_text <- function() {
  cat(readLines(file.path(.here, "query_gene.R"))[grep("^# ", readLines(file.path(.here, "query_gene.R")))][1:40],
      sep = "\n")
}
if (opts$help) { help_text(); quit(status = 0) }

idx <- read.table(opts$index, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                  quote = "", comment.char = "")

if (opts$list_fields) {
  cat(paste(colnames(idx), collapse = "\n"), "\n", sep = "")
  quit(status = 0)
}

if (length(genes) == 0) {
  message("No gene supplied. Try:  Rscript R/query_gene.R Ezh2\n",
          "See --help for options.")
  quit(status = 2)
}

# ---- resolve genes (case-insensitive on symbol, also match transcript id) --
resolve <- function(g, idx) {
  hit <- idx[idx$gene == g, , drop = FALSE]
  if (nrow(hit) == 0)
    hit <- idx[tolower(idx$gene) == tolower(g), , drop = FALSE]
  if (nrow(hit) == 0)
    hit <- idx[idx$transcript_id == g | idx$gene_id == g, , drop = FALSE]
  hit
}

res_list <- lapply(genes, resolve, idx = idx)
names(res_list) <- genes
found_mask <- vapply(res_list, nrow, integer(1)) > 0
res <- do.call(rbind, res_list[found_mask])
if (!is.null(res)) rownames(res) <- NULL
missing <- genes[!found_mask]

if (!is.null(res) && nrow(res) > 0 && !opts$all_transcripts) {
  res <- res[order(res$gene, -abs(res$log2fc)), ]
  res <- res[!duplicated(res$gene), ]
}

if (!is.null(opts$fields) && !is.null(res) && nrow(res) > 0) {
  bad <- setdiff(opts$fields, colnames(res))
  if (length(bad)) stop("Unknown field(s): ", paste(bad, collapse = ", "),
                        "\nUse --list-fields to see options.")
  res <- res[, opts$fields, drop = FALSE]
}

# ---- emit ------------------------------------------------------------------
if (is.null(res) || nrow(res) == 0) {
  if (opts$format == "json") {
    cat(jsonlite::toJSON(list(found = FALSE, missing = missing),
                         auto_unbox = TRUE, pretty = TRUE), "\n")
  } else {
    message("No matches. Missing: ", paste(missing, collapse = ", "))
  }
  quit(status = 2)
}

if (opts$format == "json") {
  out <- list(found = TRUE, n = nrow(res))
  if (length(missing)) out$missing <- missing
  out$results <- res
  cat(jsonlite::toJSON(out, auto_unbox = TRUE, pretty = TRUE, na = "null"), "\n")
} else if (opts$format == "tsv") {
  write.table(res, stdout(), sep = "\t", quote = FALSE, row.names = FALSE)
  if (length(missing)) message("# missing: ", paste(missing, collapse = ", "))
} else if (opts$format == "table") {
  print(res, row.names = FALSE)
  if (length(missing)) message("Missing: ", paste(missing, collapse = ", "))
} else {
  stop("Unknown --format: ", opts$format)
}

# ---- optional plot ---------------------------------------------------------
if (opts$plot) {
  suppressWarnings(suppressMessages({
    library(ggplot2)
    mt_source_lib()
    source(file.path(MT_LIB, "mt_plot.R"))
  }))
  plot_genes <- unique(res$gene)
  f <- mt_save_gene_plot(plot_genes, format = opts$plot_format)
  message("Plot written: ", f)
}
