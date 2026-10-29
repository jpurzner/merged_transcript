"""Validate the Python merge against the R pipeline.

Feeds the SAME R-normalized Scott + Hatten matrices into the Python
merge_datasets() and compares the kinetics to the canonical R toc_t50. The
alignment + kinetics should reproduce R (the normalization is identical because
we start from the R-normalized matrices in tests/fixtures/).

    python python/tests/validate_against_r.py
"""
import os
import sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from merged_timecourse import Dataset, merge_datasets, add_dataset  # noqa: E402

FIX = os.path.join(HERE, "fixtures")
# R time_ref day axis: E15=1 ... P0=7, P1=8, P7=14, P12=19, P14=21, P18=25, P21=28, P56=63
POSITIONS = {"E15": 1, "P0": 7, "P1": 8, "P7": 14, "P12": 19,
             "P14": 21, "P18": 25, "P21": 28, "P56": 63}


def load():
    scott = Dataset.from_files("scott_GNP", f"{FIX}/scott_normalized.tsv",
                               f"{FIX}/scott_meta.tsv", POSITIONS)
    hatten = Dataset.from_files("hatten_GNP", f"{FIX}/hatten_normalized.tsv",
                                f"{FIX}/hatten_meta.tsv", POSITIONS)
    return scott, hatten


def compare(merged, label):
    exp = pd.read_csv(f"{FIX}/expected_kinetics.tsv", sep="\t").set_index("ensembl_transcript_id")
    got = merged.kinetics
    common = got.index.intersection(exp.index)
    g, e = got.loc[common], exp.loc[common]
    print(f"\n=== {label}  (n={len(common)} common genes) ===")
    for col in ["avg_t50", "avg_fc", "avg_max", "avg_min"]:
        r = np.corrcoef(g[col].astype(float), e[col].astype(float))[0, 1]
        maxd = np.nanmax(np.abs(g[col].astype(float).values - e[col].astype(float).values))
        print(f"  {col:8s} Pearson r = {r:.5f}   max|diff| = {maxd:.3g}")
    print(f"  avg_dir  agreement = {(g['avg_dir'].values == e['avg_dir'].values).mean():.5f}")
    print(f"  t50_cat  agreement = {(g['t50_cat'].astype(str).values == e['t50_cat'].astype(str).values).mean():.5f}")
    # max_cat: compare only on genes the R table calls changing+expressed
    # (E15/P0_7/P14_P56); R additionally overlays unchanged/unexpressed, which
    # needs expression + p-values (see core.apply_gene_program).
    chg = e["max_cat"].astype(str).isin(["E15", "P0_7", "P14_P56"]).values
    mc = (g.loc[chg, "max_cat"].astype(str).values == e.loc[chg, "max_cat"].astype(str).values).mean()
    print(f"  max_cat  agreement (on {chg.sum()} changing genes) = {mc:.5f}")
    return np.corrcoef(g['avg_t50'].astype(float), e['avg_t50'].astype(float))[0, 1]


if __name__ == "__main__":
    scott, hatten = load()

    # Job 3: merge two (N) datasets at once
    merged = merge_datasets([scott, hatten])
    r_t50 = compare(merged, "merge_datasets([scott, hatten])")

    # Job 2: build from one, then ADD the second -> must match the joint merge
    start = merge_datasets([scott])
    added = add_dataset(start, hatten)
    same_kin = np.allclose(added.kinetics[["avg_t50", "avg_fc"]].values,
                           merged.kinetics[["avg_t50", "avg_fc"]].values, equal_nan=True)
    print(f"\nadd_dataset(merge([scott]), hatten) == merge([scott, hatten]) : {same_kin}")

    ok = r_t50 > 0.999 and same_kin
    print("\nRESULT:", "PASS" if ok else "CHECK", "(t50 r =", round(r_t50, 5), ")")
    sys.exit(0 if ok else 1)
