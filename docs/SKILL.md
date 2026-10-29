---
name: gnp-transcript-timecourse
description: >
  Query the developmental expression kinetics of a mouse gene during
  cerebellar granule-neuron-precursor (GNP) differentiation. Given a gene
  symbol, returns when it turns on/off (t50), how much it changes (log2 fold
  change), direction, time of peak expression, and a kinetic class — merged
  across two independent time courses (Scott RNA-seq E15–P14 + Hatten TRAP
  array P0–P56). Can also render a timecourse plot. Use when asked whether a
  gene is early/late, rising/falling, or expressed during GNP/cerebellar
  development, or to compare the timing of several genes.
---

# GNP transcript timecourse

## What this is

Two independent mouse cerebellar developmental transcriptomic time courses,
merged onto a common postnatal-day **pseudo-time axis** and reduced to
per-transcript kinetic metrics:

| Dataset      | Platform         | Time points              | Samples |
|--------------|------------------|--------------------------|---------|
| `scott_GNP`  | RNA-seq          | E15, P1, P7, P14         | 11      |
| `hatten_GNP` | TRAP microarray  | P0, P7, P12, P18, P21, P56 | 18    |

Data sources: `hatten_GNP` is the Hatten-lab granule-neuron NeuroD1-TRAP data,
GEO GSE74400 (Zhu et al., *Neuron* 89:100–112, 2016, doi 10.1016/j.neuron.2015.11.021);
`scott_GNP` RNA-seq was generated for Purzner et al. 2025. See `references/`.

Pseudo-time axis: **E15 → P0 → P1 → P7 → P12 → P14 → P18 → P21 → P56**.

The two series are interpolated onto this axis, aligned (shifted so their
overlapping days agree), row-normalized to log2 fold change vs each gene's own
mean, and summarized. Coverage: **24,379 genes / 24,525 transcripts**.

## How to query it (the API)

A dependency-light CLI — the intended entry point for agents:

```bash
Rscript R/query_gene.R <GENE> [<GENE> ...] [options]
```

Default output is JSON. Examples:

```bash
Rscript R/query_gene.R Ezh2
Rscript R/query_gene.R Ezh2 Gabra6 Sox2 --format table
Rscript R/query_gene.R Atoh1 --fields gene,t50,log2fc,direction_label,kinetic_class
Rscript R/query_gene.R Ezh2 Neurod1 --plot --plot-format png
```

Options: `--format json|table|tsv` · `--fields a,b,c` · `--all-transcripts`
(every transcript, not just the largest-change one per gene) · `--plot`
`--plot-format pdf|png` · `--list-fields` · `--help`.

Exit status `0` if any gene resolved, `2` if none did. Lookup is
case-insensitive on the gene symbol and also matches an Ensembl transcript or
gene id.

Example JSON (trimmed):

```json
{
  "found": true, "n": 1,
  "results": [{
    "gene": "Ezh2", "transcript_id": "ENSMUST00000114623",
    "t50": 22.01, "log2fc": 3.14, "direction": -1, "direction_label": "down",
    "time_of_max": 1, "max_category": "E15", "t50_category": "Late",
    "kinetic_class": "Max_E15 t50_Late", "fc_class": ">2 log2FC",
    "expressed": true
  }]
}
```

For bulk / programmatic access, read the slim table directly:
`data/timecourse_query_index.tsv` (one row per transcript, tab-separated,
gene symbol in column `gene`).

### Numerical timecourse values

When you need the actual shape of the trajectory (not just the summary), get
the per-time-point values behind the plot:

```bash
Rscript R/timecourse_values.R <GENE> [<GENE> ...] [options]
```

Options: `--format json|table|tsv` · `--replicates` (raw per-replicate values)
· `--series` (JSON: compact `{time, mean}` per dataset) · `--help`.

Each gene returns its `summary` plus a `timecourse` array of
`{dataset, time, pseudo_time_index, n, mean, sd, sem}` — the aligned **log2
fold change vs the gene's own mean** at each of the 9 sampled time points, for
`scott_GNP` (E15,P1,P7,P14) and `hatten_GNP` (P0,P7,P12,P18,P21,P56). Read the
`time` labels (E15…P56), not `pseudo_time_index`, as the axis. Example:

```json
{ "unit": "log2 fold change vs gene mean, aligned across datasets",
  "genes": [{ "gene": "Ccnd1", "transcript_id": "uc009kqt.1",
    "summary": { "t50_day": "P12", "time_of_max": "P1", "direction": "down" },
    "timecourse": [
      {"dataset":"scott_GNP","time":"E15","mean":0.22,"sd":0.08,"n":3},
      {"dataset":"scott_GNP","time":"P1","mean":1.36,"sd":0.38,"n":2},
      {"dataset":"hatten_GNP","time":"P56","mean":-1.36,"sd":0.07,"n":3} ] }]}
```

Exit status `0` if any gene resolved to an aligned timecourse, `2` otherwise.
A gene with no transcript in the aligned matrices is reported under
`not_aligned` (its summary metrics are still available via `query_gene.R`).

## Field dictionary

**Two time axes — do not confuse them:**

- The interpolated metrics `t50`, `time_of_max`, `time_of_min` are indices on
  the **63-point daily axis** E15,E16,…,E20,P0,P1,…,P56, i.e.
  **1 = E15, 7 = P0, 8 = P1, 14 = P7, 19 = P12, 21 = P14, 25 = P18, 28 = P21,
  63 = P56** (fractional between days). So `t50 = 22` ≈ **P15**, and
  `time_of_max = 8` = **P1**. Use the companion `*_day` fields
  (`t50_day`, `time_of_max_day`, `time_of_min_day`) for the human-readable
  label — you rarely need the raw index.
- The 9 actually-sampled time points (`MT_PSEUDO_TIME`: E15, P0, P1, P7, P12,
  P14, P18, P21, P56) are what the observed timecourse values (below) are
  reported at.

| Field | Meaning |
|-------|---------|
| `gene` | MGI gene symbol (the lookup key) |
| `transcript_id`, `gene_id` | Ensembl transcript / gene id |
| `chr`, `start`, `end`, `strand` | genomic location (mm, Ensembl coords) |
| `t50` / `t50_day` | daily-axis index / day label at which the gene reaches 50% of its total change — *when it switches*. Lower = earlier. `t50_scott`/`t50_hatten` are the per-dataset raw indices. |
| `log2fc` | magnitude of change over development, log2, **unsigned**. `log2fc_scott`/`log2fc_hatten` per dataset. |
| `signed_log2fc` | `log2fc × direction` (positive = up, negative = down) |
| `direction` / `direction_label` | `+1`/`"up"` = increases with differentiation; `-1`/`"down"` = decreases |
| `time_of_max` / `time_of_max_day`, `time_of_min` / `time_of_min_day` | daily-axis index / day label of peak / trough expression |
| `max_category` | time-of-peak bin: `E15` (peaks embryonic), `P0_7` (peaks early postnatal), `P14_P56` (peaks late), or `unchanged` / `unexpressed` |
| `t50_category` | `Early` (t50 < 6 → switches by ~E20/P0), `Mid` (6–18 → ~P0–P12), `Late` (>18 → after ~P12) |
| `kinetic_class` | combined `Max_<maxcat> t50_<t50cat>` label — the primary grouping used in figures |
| `gene_program` | developmental program by time of peak (Purzner et al. 2025): `embryonic` (peak E15), `proliferation` (peak P0–P7; e.g. Gli1, Ccnd1), `differentiation` (peak P14–P56; e.g. Cacna1a, Cbx7), else `unchanged`/`unexpressed`. For the paper's strict gene sets (783/440/787) also require `fc_class == ">2 log2FC"`. |
| `fc_class` | `0-1 log2FC` / `1-2 log2FC` / `>2 log2FC` |
| `expressed` | TRUE if called expressed in either dataset (binarized max expression) |
| `max_expr_scott`, `max_expr_hatten` | max expression level (length-normalized counts / array intensity) |
| `min_padj_scott`, `min_padj_hatten` | smallest adjusted p-value across that dataset's contrasts (significance of change) |

A gene is **"changing"** if `log2fc >= 1` and it is significant
(`min_padj < 0.001`) in at least one dataset; otherwise it is `unchanged`.
Genes not called expressed are `unexpressed`. These two labels fill
`max_category` / `kinetic_class` for non-changing genes.

## Biological reading

- **`direction = down`, `t50` late, peaks E15** → progenitor / proliferation
  program shutting off during differentiation (e.g. *Ezh2*, PRC2).
- **`direction = up`, peaks P14–P56** → maturation / mature-neuron program
  (e.g. *Gabra6*, a mature granule-neuron marker).
- Comparing `t50` across genes orders them along the differentiation
  trajectory — smaller `t50` switches earlier.

## Generating figures

```bash
Rscript R/03_timecourse_plots.R          # standard figure set -> figures/
Rscript R/query_gene.R <genes> --plot    # ad-hoc per-gene timelines
```

Standard outputs in `figures/`:
- `kinetic_class_timelines.pdf` — mean ±SD trajectory of every changing
  transcript grouped by kinetic class, both datasets overlaid, median t50 marked.
- `example_genes_timeline.pdf` — illustrative marker genes.
- `t50_vs_log2fc.pdf` — t50 vs signed fold change by time-of-max category.

Gene-classification Sankey (Figure S1 of Purzner et al. 2025), the stepwise
filter that defines the three gene programs — reproduces the published set
sizes (783 embryonic / 440 proliferation / 787 differentiation):

```bash
Rscript R/figures/sankey_classification.R   # -> figures/sankey_gene_classification.{pdf,png}
```

Needs `ggsankey` (`devtools::install_github("davidsjoberg/ggsankey")`) and
`data/gene_classification.tsv`. The verbatim source Rmds are in `analysis/`.

## Provenance & reproducibility

- Canonical merged table: `data/toc_t50_2023_05_17.txt` (115 columns; the slim
  index is derived from it). This is the source of truth for the API and figures.
- Aligned expression matrices for plotting: `data/log2_norm_list_shift.RData`;
  the full alignment datasets (normalization → interpolation → alignment
  intermediates) are bundled in `data/alignment_datasets.RData` — see
  `data/README.md`.
- Rebuild from raw in two stages: `Rscript R/run_pipeline.R --step normalize`
  (DESeq2-normalize + merge) then `--step align` (interpolate + align + t50),
  or `--step merge` for both. Needs the raw inputs symlinked in `data/raw/`
  (see `data/raw/README.md`). The rebuild
  reproduces the kinetic metrics of the canonical table **exactly** (t50,
  fold change, time-of-max: Pearson r = 1; direction: 100% agreement) and
  ~93% of the categorical `kinetic_class` calls — the difference is confined
  to the expression-binarization cutoff (the source Rmd's gamma-mixture step
  was left `eval=FALSE`; the rebuild uses a stable normal-mixture / density
  fallback), which only reclassifies borderline genes near the
  expressed/unchanged boundary.
- Ported from `GNP_MB_integrate/Ezh2_2022/1_transcript_timecourse_3.Rmd`
  (also copied verbatim into `analysis/`).
- This dataset underlies **Purzner et al. 2025**, *Ezh2 Delays Activation of
  Differentiation Genes During Normal Cerebellar Granule Neuron Development and
  in Medulloblastoma*, eLife (DOI 10.7554/eLife.105339.1; bioRxiv
  2024.11.21.624171). See `references/` for the PDF, citation, and how
  `gene_program` maps to the paper's gene sets.
