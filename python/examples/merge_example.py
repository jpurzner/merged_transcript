"""Worked example: recreate the merge in Python, with N datasets and adding one.

    python python/examples/merge_example.py

Uses the R-normalized fixtures in tests/fixtures/ (Scott RNA-seq + Hatten TRAP).
Demonstrates all three jobs:
  1. recreate the merge in Python   -> merge_datasets([scott, hatten])
  2. add a dataset to an existing merge -> merged.add_dataset(late)
  3. multiple datasets at once        -> merge_datasets([scott, hatten, late])
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from merged_timecourse import Dataset, merge_datasets, core  # noqa: E402

FIX = os.path.join(HERE, "..", "tests", "fixtures")
POSITIONS = {"E15": 1, "P0": 7, "P1": 8, "P7": 14, "P12": 19,
             "P14": 21, "P18": 25, "P21": 28, "P56": 63}
DAY_LABELS = {v: k for k, v in POSITIONS.items()}

scott = Dataset.from_files("scott_GNP", f"{FIX}/scott_normalized.tsv",
                           f"{FIX}/scott_meta.tsv", POSITIONS)
hatten = Dataset.from_files("hatten_GNP", f"{FIX}/hatten_normalized.tsv",
                            f"{FIX}/hatten_meta.tsv", POSITIONS)

# --- Job 1: recreate the merge --------------------------------------------
merged = merge_datasets([scott, hatten])
print("Merged:", merged.dataset_names, "| genes:", len(merged.genes),
      "| axis:", int(merged.grid[0]), "-", int(merged.grid[-1]))

for g in ["ENSMUST00000114623", "ENSMUST00000155218"]:   # Ezh2, Gabra6
    if g in merged.kinetics.index:
        k = merged.kinetics.loc[g]
        print(f"  {g}: t50={k.avg_t50:.1f} ({DAY_LABELS.get(round(k.avg_t50),'?')}) "
              f"fc={k.avg_fc:.2f} dir={k.avg_dir:+d} {k.max_t50_cat}")

# --- Job 2: add a dataset to an existing merge ----------------------------
# A synthetic "late" array: Hatten's late timepoints as a separate dataset,
# overlapping the existing merge at P18/P21/P56 (its anchor).
late_cols = [s for s, t in zip(hatten.expr.columns, hatten.sample_times)
             if t in ("P18", "P21", "P56")]
late = Dataset("late_GNP", hatten.expr[late_cols],
               [t for t in hatten.sample_times if t in ("P18", "P21", "P56")],
               POSITIONS)
merged3 = merged.add_dataset(late)
print("\nAfter add_dataset(late):", merged3.dataset_names,
      "| genes:", len(merged3.genes))

# --- Job 3: multiple datasets at once (same result) -----------------------
merged3b = merge_datasets([scott, hatten, late])
print("merge_datasets([scott, hatten, late]):", merged3b.dataset_names)

# --- gene_program overlay (needs expression; here: 'changing' = fc>=1) ----
expressed = merged.kinetics["avg_fc"].to_numpy() >= 0   # placeholder: all kept
kin = core.apply_gene_program(merged.kinetics,
                              expressed=merged.kinetics["avg_fc"] >= 1,
                              fc_cutoff=1.0)
print("\ngene_program counts:")
print(kin["gene_program"].value_counts(dropna=False).to_string())
