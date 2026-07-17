#!/usr/bin/env Rscript
# run_pipeline.R -- orchestrate the merged_transcript pipeline.
#
# Usage:
#   Rscript R/run_pipeline.R --step <step>
#
# Steps:
#   merge     rebuild toc_t50 from raw inputs        (R/01_load_and_merge.R)
#   plots     generate the standard figures          (R/02_timecourse_plots.R)
#   index     build the slim query index             (R/03_build_query_index.R)
#   all       index -> plots  (uses the shipped merged table; does NOT rerun
#             the heavy DESeq2 merge -- run `--step merge` explicitly for that)
#
# The API (R/query_gene.R) is a standalone entry point, not a pipeline step.

args_full <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "",
          grep("^--file=", args_full, value = TRUE)[1])))

argv <- commandArgs(trailingOnly = TRUE)
step <- NULL
if (length(argv) >= 2 && argv[1] == "--step") step <- argv[2]
if (is.null(step)) {
  message("Usage: Rscript R/run_pipeline.R --step {merge|plots|index|all}")
  quit(status = 2)
}

run <- function(script) {
  message("\n=== ", script, " ===")
  code <- system2("Rscript", file.path(here, script))
  if (code != 0) stop(script, " failed (exit ", code, ")")
}

switch(step,
  merge = run("01_load_and_merge.R"),
  plots = run("02_timecourse_plots.R"),
  index = run("03_build_query_index.R"),
  all   = { run("03_build_query_index.R"); run("02_timecourse_plots.R") },
  stop("Unknown step: ", step,
       "  (choose merge, plots, index, or all)")
)
message("\nDone: step '", step, "'")
