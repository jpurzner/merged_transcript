#!/usr/bin/env Rscript
# run_pipeline.R -- orchestrate the merged_transcript pipeline.
#
# Usage:
#   Rscript R/run_pipeline.R --step <step>
#
# Steps:
#   normalize  load raw, DESeq2-normalize Scott + 2^ Hatten, merge on
#              annotation                              (R/01_normalize_and_merge.R)
#   align      interpolate + align onto the pseudo-time axis, build the
#              alignment datasets + toc_t50            (R/02_align_timecourse.R)
#   merge      normalize -> align  (the full raw rebuild)
#   plots      generate the standard figures           (R/03_timecourse_plots.R)
#   index      build the slim query index              (R/04_build_query_index.R)
#   all        index -> plots  (uses the shipped data; does NOT rerun the heavy
#              normalize/align rebuild -- use `--step merge` for that)
#
# The API (R/query_gene.R) is a standalone entry point, not a pipeline step.

args_full <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "",
          grep("^--file=", args_full, value = TRUE)[1])))

argv <- commandArgs(trailingOnly = TRUE)
step <- NULL
if (length(argv) >= 2 && argv[1] == "--step") step <- argv[2]
if (is.null(step)) {
  message("Usage: Rscript R/run_pipeline.R --step {normalize|align|merge|plots|index|all}")
  quit(status = 2)
}

run <- function(script) {
  message("\n=== ", script, " ===")
  code <- system2("Rscript", file.path(here, script))
  if (code != 0) stop(script, " failed (exit ", code, ")")
}

switch(step,
  normalize = run("01_normalize_and_merge.R"),
  align     = run("02_align_timecourse.R"),
  merge     = { run("01_normalize_and_merge.R"); run("02_align_timecourse.R") },
  plots     = run("03_timecourse_plots.R"),
  index     = run("04_build_query_index.R"),
  all       = { run("04_build_query_index.R"); run("03_timecourse_plots.R") },
  stop("Unknown step: ", step,
       "  (choose normalize, align, merge, plots, index, or all)")
)
message("\nDone: step '", step, "'")
