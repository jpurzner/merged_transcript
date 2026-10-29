# analysis/ — source code for the RNA-seq figures

Verbatim source R Markdown for the **RNA-seq-based figures** of Purzner et al.
2025 (eLife 10.7554/eLife.105339.1; see `../references/`). Copied from
`GNP_MB_integrate/Ezh2_2022/`. The ChIP-seq / histone figures (Fig 2
co-occupancy, Fig 8 MB histones, the promoter-histone analyses) are **not**
RNA-seq and live with the ChIP work in the original directory and other repos.

These `.Rmd` files are the authoritative figure code. They were written to run
interactively and read data from hard-coded `~/Dropbox/...` paths, so most are
**reference code**, not runnable as-is in this repo — only the bulk-timecourse
figures have been ported to standalone scripts that run from the shipped data
(see "Ported & runnable" below).

## Files → figures → inputs

| Rmd | Paper figure(s) | What it does | Key inputs |
|-----|-----------------|--------------|-----------|
| `1_transcript_timecourse_3.Rmd` | Fig 1 timecourse | Build the merged Scott+Hatten timecourse → `toc_t50` (normalize, align, t50/fold-change, kinetic classes); timecourse profile plots | raw Scott RNA-seq counts, Hatten TRAP array, annotation |
| `3_gene_classification.Rmd` | **Fig S1** (Sankey) + differentiation t50-vs-log2FC scatter | Stepwise gene-program classification; the Sankey; labelled scatter of P14–P56 genes | `sce_contaminant_full.txt` (= `toc_t50` + scRNA-seq lineage cols) |
| `3_gene_exclusion_classification.Rmd` | (feeds Fig S1) | Lineage / cell-cycle / contamination exclusion that produces the `group_1` / `exclude_cat` annotations | Taylor scRNA-seq classification + `toc_t50` |
| `2_scRNAseq_2.Rmd` | scRNA-seq panels | Taylor cerebellar scRNA-seq cell-type assignment used to flag non-GNP-lineage / cell-cycle genes | Seurat/`SingleCellExperiment` objects (not in repo) |
| `7_GNP_diff_go_terms_2024_10_29.Rmd` | GO-term panels | GO enrichment of the embryonic / proliferation / differentiation programs | the gene programs + GO databases |
| `EZH2_diff_fig3.Rmd` | Fig 3 | EZH2 / H3K27me3 relationship to differentiation-gene activation (RNA-seq + marks) | `toc_t50` + mark data |
| `6_Human_MB.Rmd` | Fig 5 | Expression of the gene programs across human MB subtypes | human MB expression tables |
| `survival_analysis.Rmd` | Fig 5 survival | SHH-MB survival by H3K27me3-bound differentiation genes | clinical + expression + mark data |
| `t50window_MF_2024_10_24.Rmd` | supp. | t50-window gene extraction / sensitivity analysis | `toc_t50` |
| `trim_upper_tail_rationale.Rmd` | methods | Rationale for the winsorization / upper-tail trim in binarization | — |
| `sample_RNAseq.Rmd` | methods | Sample RNA-seq processing pipeline | raw sequencing |
| `human2mouse.Rmd` | helper | Human↔mouse ortholog mapping | ortholog tables |

## Ported & runnable (from this repo's shipped data)

- **Fig 1 timecourse** → `R/01_normalize_and_merge.R` + `R/02_align_timecourse.R`
  (rebuild `toc_t50`) and `R/03_timecourse_plots.R` (profile plots, t50-vs-FC
  scatter).
- **Fig S1 Sankey** → `R/figures/sankey_classification.R` (reproduces the
  published 783 / 440 / 787 gene-set sizes).

The remaining Rmds depend on data that is not shipped here (scRNA-seq objects,
ChIP marks, human MB / clinical tables); they are included as the complete
source for those figures. Several are discussed across other repos in this
project.
