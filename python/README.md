# merged_timecourse (Python)

A Python recreation of the R merged_transcript pipeline, **generalized to any
number of datasets** and to **incrementally adding** a dataset to an existing
merge. Normalize → interpolate onto a shared pseudo-time axis → align (N
datasets) → per-gene kinetics (t50, fold change, direction, program categories).

Validated against the R pipeline: on the same normalized inputs it reproduces
the canonical `toc_t50` kinetics exactly — `avg_max` identical, `avg_fc` / `t50`
Pearson r = 1, direction / t50_cat / (changing-gene) max_cat 100% agreement.

## Install

```bash
pip install -r python/requirements.txt   # numpy, pandas, scipy, matplotlib
```

The example, CLI and validation read the Scott+Hatten fixtures, which are large
and git-ignored — generate them once from the R data:

```bash
Rscript python/tests/export_fixtures.R
```

## The three capabilities

```python
from merged_timecourse import Dataset, merge_datasets, add_dataset

POSITIONS = {"E15":1,"P0":7,"P1":8,"P7":14,"P12":19,"P14":21,"P18":25,"P21":28,"P56":63}

scott  = Dataset.from_files("scott_GNP",  "scott_normalized.tsv",  "scott_meta.tsv",  POSITIONS)
hatten = Dataset.from_files("hatten_GNP", "hatten_normalized.tsv", "hatten_meta.tsv", POSITIONS)

# 1. recreate the merge (and 3. handle N datasets) ------------------------
merged = merge_datasets([scott, hatten])          # or [ds1, ds2, ds3, ...]
merged.kinetics.loc["ENSMUST00000114623"]         # t50 / fold change / program

# 2. add a dataset to an existing merge -----------------------------------
merged = merged.add_dataset(third_dataset)        # anchors onto the existing scale
```

`merge_datasets([...])` and `add_dataset(...)` are the same routine: the first
dataset is the alignment reference, and each later dataset is interpolated onto
the shared axis and shifted (per gene) to agree with the datasets already placed
— so adding one more is just extending that sequence. Genes are intersected
across datasets; the pseudo-time axis is built from the union of all datasets'
time positions.

### A Dataset

- `expr`  — genes (rows) × samples (columns), **normalized** expression
  (DESeq2-normalized counts, `2**rma` array intensity, …). `normalize.py` has
  `size_factor_normalize` (median-of-ratios) and `unlog2` helpers; or supply
  already-normalized matrices.
- `sample_times` — the timepoint label of each sample column.
- `positions` — `{timepoint: numeric axis position}` (e.g. developmental day).
  This is what lets arbitrary datasets with different timepoints share one axis.

### Gene programs

`core.apply_gene_program(kinetics, expressed, significant=None)` overlays the
`unchanged` / `unexpressed` calls and adds the `gene_program` column
(`embryonic` / `proliferation` / `differentiation`), matching the R annotation
(Purzner et al. 2025).

## CLI

```bash
python python/query_gene.py ENSMUST00000114623            # Ezh2 kinetics (JSON)
python python/query_gene.py ENSMUST00000155218 --values   # + aligned timecourse
```

## Layout

```
python/
├── merged_timecourse/
│   ├── dataset.py    Dataset (expression + timepoints + axis positions)
│   ├── normalize.py  median-of-ratios / unlog2 helpers
│   ├── core.py       row_log2fc, interpolate, align (N datasets), kinetics, gene_program
│   ├── merge.py      MergedTimecourse, merge_datasets(), add_dataset()
│   └── plot.py       per-gene timecourse plot (matplotlib)
├── query_gene.py     query CLI (counterpart of R/query_gene.R)
├── examples/merge_example.py
└── tests/
    ├── export_fixtures.R      regenerate fixtures from the R data
    ├── fixtures/              R-normalized Scott+Hatten + expected kinetics
    └── validate_against_r.py  reproduce toc_t50 from the same inputs
```

## Validate

```bash
Rscript python/tests/export_fixtures.R      # (re)generate fixtures from R data
python  python/tests/validate_against_r.py  # compare Python merge to R toc_t50
```
