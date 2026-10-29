#!/usr/bin/env Rscript
# export_fixtures.R -- regenerate the Python validation fixtures from the R data.
#
# Writes the R-normalized Scott + Hatten matrices, their sample->timepoint maps,
# and the canonical toc_t50 kinetics into python/tests/fixtures/. Feeding the
# SAME normalized matrices into the Python merge isolates the port validation to
# the alignment + kinetics (see validate_against_r.py).
#
# Usage (from repo root):  Rscript python/tests/export_fixtures.R

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])))
root <- normalizePath(file.path(here, "..", ".."))
source(file.path(root, "R", "lib", "common.R"))

if (!file.exists(MT_NORMMERGE_RDATA))
  stop("Run `Rscript R/run_pipeline.R --step normalize` first.")
load(MT_NORMMERGE_RDATA)                 # RNAseq_list, meta_list, pval_all
fix <- file.path(here, "fixtures"); dir.create(fix, showWarnings = FALSE, recursive = TRUE)

write_mat <- function(m, path)
  write.table(data.frame(transcript_id = rownames(m), m, check.names = FALSE),
              path, sep = "\t", quote = FALSE, row.names = FALSE)
write_meta <- function(cols, meta, path)
  write.table(data.frame(sample = cols, timepoint = as.character(meta)),
              path, sep = "\t", quote = FALSE, row.names = FALSE)

write_mat(RNAseq_list[[1]], file.path(fix, "scott_normalized.tsv"))
write_mat(RNAseq_list[[2]], file.path(fix, "hatten_normalized.tsv"))
write_meta(colnames(RNAseq_list[[1]]), MT_SCOTT_META,  file.path(fix, "scott_meta.tsv"))
write_meta(colnames(RNAseq_list[[2]]), MT_HATTEN_META, file.path(fix, "hatten_meta.tsv"))

d <- mt_load_toc()
write.table(d[, c("ensembl_transcript_id","avg_t50","avg_fc","avg_dir",
                  "avg_max","avg_min","max_cat","t50_cat")],
            file.path(fix, "expected_kinetics.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
message("Fixtures written to ", fix)
