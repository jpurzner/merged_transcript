# Raw inputs

These files are the raw inputs to `R/01_load_and_merge.R`. They are **not
committed** — in this working copy they are symlinks to their original
locations, and they are git-ignored. To rebuild the merged table elsewhere,
place files with these names here (or point the matching env var at them).

| File (in this dir)              | Env override      | Original location | Description |
|---------------------------------|-------------------|-------------------|-------------|
| `scott_rnaseq_counts.txt`       | `MT_SCOTT_COUNTS` | `GNP_MB_integrate/tss_model/Scott_RNAseq_tss_model_uniq_htseq_counts.txt` | Scott GNP RNA-seq htseq counts (E15,P1,P7,P14) |
| `scott_time_series_name.txt`    | `MT_SCOTT_TIMES`  | `GNP_MB_integrate/2016_06_18_integration/table_2017_12_11/time_series_name.txt` | which Scott columns form the time series |
| `hatten_trap_microarray.txt`    | `MT_HATTEN_EXPR`  | `GNP_MB_integrate/Ezh2_2022/hatten_trap_microarray_2023_04_18.txt` | Hatten TRAP array RMA values (P0…P56), log2 |
| `hatten_trap_microarray_pval.txt`| `MT_HATTEN_PVAL` | `GNP_MB_integrate/Ezh2_2022/hatten_trap_microarray_pval_2023_04_18.txt` | Hatten limma per-contrast adjusted p-values |
| `gnp_tss_model_anno.txt`        | `MT_ANNO`         | `GNP_MB_integrate/tss_model/GNP_tss_model_anno_2018_12_08.txt` | Ensembl transcript↔gene annotation + exon blocks |

The shipped `data/toc_t50_2023_05_17.txt` was produced from exactly these
inputs and is the canonical merged table; you only need the raw files to
re-derive it.
