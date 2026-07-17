#!/usr/bin/env Rscript
# timecourse_values.R -- numerical timecourse API over the merged GNP data.
#
# Returns the per-time-point aligned log2 fold-change values behind the plots,
# so an agent can reason about the shape of a gene's developmental trajectory
# without seeing the image. Values are mean +/- spread across replicates at
# each developmental time point, for each dataset.
#
# USAGE
#   Rscript R/timecourse_values.R <GENE> [<GENE> ...] [options]
#
# OPTIONS
#   --format json|table|tsv   output format (default: json)
#   --replicates              include the raw per-replicate values
#   --series                  JSON only: also emit a compact per-dataset series
#                             {time, mean} — handy for a quick shape read
#   --index <path>            override the query index path
#   --help                    show this help
#
# EXAMPLES
#   Rscript R/timecourse_values.R Ccnd1
#   Rscript R/timecourse_values.R Ezh2 Gabra6 --format table
#   Rscript R/timecourse_values.R Atoh1 --replicates
#
# OUTPUT (json): per gene, the summary metrics plus a `timecourse` array of
#   {dataset, time, pseudo_time_index, n, mean, sd, sem}. Units: log2 fold
#   change vs the gene's own mean, aligned across the two datasets on the
#   E15..P56 pseudo-time axis.
#
# EXIT STATUS
#   0  at least one gene resolved to an aligned timecourse
#   2  none resolved / usage error

suppressWarnings(suppressMessages({
  .argv_full <- commandArgs(trailingOnly = FALSE)
  .here <- dirname(normalizePath(sub("^--file=", "",
             grep("^--file=", .argv_full, value = TRUE)[1])))
  source(file.path(.here, "lib", "common.R"))
  source(file.path(MT_LIB, "mt_plot.R"))    # mt_load_index()
  source(file.path(MT_LIB, "mt_values.R"))
}))

# ---- arg parsing -----------------------------------------------------------
argv <- commandArgs(trailingOnly = TRUE)
opts <- list(format = "json", replicates = FALSE, series = FALSE,
             index = MT_QUERY_INDEX, help = FALSE)
genes <- character(0)
i <- 1
while (i <= length(argv)) {
  a <- argv[i]
  if (a %in% c("--help", "-h")) opts$help <- TRUE
  else if (a == "--replicates") opts$replicates <- TRUE
  else if (a == "--series") opts$series <- TRUE
  else if (a == "--format") { i <- i + 1; opts$format <- argv[i] }
  else if (a == "--index")  { i <- i + 1; opts$index <- argv[i] }
  else if (grepl("^--", a)) stop("Unknown option: ", a)
  else genes <- c(genes, a)
  i <- i + 1
}
if (opts$help) {
  cat(grep("^#", readLines(file.path(.here, "timecourse_values.R")), value = TRUE)[1:36],
      sep = "\n"); cat("\n"); quit(status = 0)
}
if (length(genes) == 0) {
  message("No gene supplied. Try:  Rscript R/timecourse_values.R Ccnd1")
  quit(status = 2)
}

idx  <- read.table(opts$index, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                   quote = "", comment.char = "")
expr <- mt_load_expr()
tc   <- mt_timecourse_values(genes, expr = expr, idx = idx,
                             with_replicates = opts$replicates)
resn <- attr(tc, "resolution")

if (nrow(tc) == 0) {
  if (opts$format == "json")
    cat(jsonlite::toJSON(list(found = FALSE, missing = resn$missing,
                              not_aligned = resn$not_aligned),
                         auto_unbox = TRUE, pretty = TRUE), "\n")
  else
    message("No aligned timecourse. Unknown: ",
            paste(resn$missing, collapse = ", "),
            if (length(resn$not_aligned))
              paste0("  |  no aligned transcript: ",
                     paste(resn$not_aligned, collapse = ", ")) else "")
  quit(status = 2)
}

# round for readability
for (cc in c("mean", "sd", "sem")) tc[[cc]] <- round(tc[[cc]], 4)

# ---- emit ------------------------------------------------------------------
if (opts$format == "tsv") {
  flat <- tc; flat$replicates <- NULL
  write.table(flat, stdout(), sep = "\t", quote = FALSE, row.names = FALSE)
} else if (opts$format == "table") {
  flat <- tc; flat$replicates <- NULL
  print(flat, row.names = FALSE)
} else if (opts$format == "json") {
  meta <- read.table(opts$index, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                     quote = "", comment.char = "")
  per_gene <- lapply(unique(tc$gene), function(g) {
    sub <- tc[tc$gene == g, ]
    s   <- meta[meta$transcript_id == sub$transcript_id[1], ]
    obj <- list(
      gene = g,
      transcript_id = sub$transcript_id[1],
      summary = list(t50 = s$t50, t50_day = mt_day_label(s$t50),
                     log2fc = s$log2fc,
                     direction = s$direction_label,
                     time_of_max = mt_day_label(s$time_of_max),
                     kinetic_class = s$kinetic_class),
      timecourse = sub[, c("dataset","time","pseudo_time_index","n","mean","sd","sem",
                           intersect("replicates", names(sub)))]
    )
    if (opts$series) {
      obj$series <- lapply(unique(sub$dataset), function(d) {
        ss <- sub[sub$dataset == d, ]
        list(dataset = d, time = ss$time, mean = ss$mean)
      })
    }
    obj
  })
  out <- list(found = TRUE,
              unit = "log2 fold change vs gene mean, aligned across datasets",
              pseudo_time = MT_PSEUDO_TIME,
              n = length(per_gene))
  if (length(resn$missing))     out$missing     <- resn$missing
  if (length(resn$not_aligned)) out$not_aligned <- resn$not_aligned
  out$genes <- per_gene
  cat(jsonlite::toJSON(out, auto_unbox = TRUE, pretty = TRUE, na = "null"), "\n")
} else stop("Unknown --format: ", opts$format)

if (opts$format != "json" && length(c(resn$missing, resn$not_aligned)))
  message("# unresolved: ", paste(c(resn$missing, resn$not_aligned), collapse = ", "))
