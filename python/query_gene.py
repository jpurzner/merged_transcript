#!/usr/bin/env python3
"""Python query CLI over a merged timecourse — the counterpart of R/query_gene.R.

Builds the merge from normalized dataset files, then reports a gene's kinetics
and (optionally) its per-position aligned timecourse.

USAGE
  python python/query_gene.py <GENE> [<GENE> ...] [options]

OPTIONS
  --fixtures DIR     directory with scott/hatten normalized + meta TSVs
                     (default: python/tests/fixtures)
  --format json|table            output format (default: json)
  --values           also emit the aligned timecourse values per gene
  --help

EXAMPLES
  python python/query_gene.py ENSMUST00000114623            # Ezh2
  python python/query_gene.py ENSMUST00000155218 --values   # Gabra6
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from merged_timecourse import Dataset, merge_datasets  # noqa: E402

POSITIONS = {"E15": 1, "P0": 7, "P1": 8, "P7": 14, "P12": 19,
             "P14": 21, "P18": 25, "P21": 28, "P56": 63}
DAY_AXIS = ["E15", "E16", "E17", "E18", "E19", "E20"] + [f"P{i}" for i in range(0, 57)]


def day_label(idx):
    if idx != idx:   # NaN
        return None
    i = min(max(int(round(idx)), 1), len(DAY_AXIS))
    return DAY_AXIS[i - 1]


def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument("genes", nargs="*")
    ap.add_argument("--fixtures", default=os.path.join(HERE, "tests", "fixtures"))
    ap.add_argument("--format", default="json", choices=["json", "table"])
    ap.add_argument("--values", action="store_true")
    args = ap.parse_args()
    if not args.genes:
        print("No gene supplied. e.g. python python/query_gene.py ENSMUST00000114623",
              file=sys.stderr)
        sys.exit(2)

    f = args.fixtures
    scott = Dataset.from_files("scott_GNP", f"{f}/scott_normalized.tsv",
                               f"{f}/scott_meta.tsv", POSITIONS)
    hatten = Dataset.from_files("hatten_GNP", f"{f}/hatten_normalized.tsv",
                                f"{f}/hatten_meta.tsv", POSITIONS)
    merged = merge_datasets([scott, hatten])

    found, missing = [], []
    for g in args.genes:
        if g not in merged.kinetics.index:
            missing.append(g)
            continue
        k = merged.kinetics.loc[g]
        rec = dict(gene=g, avg_t50=round(float(k.avg_t50), 3),
                   t50_day=day_label(k.avg_t50), avg_fc=round(float(k.avg_fc), 3),
                   direction=int(k.avg_dir),
                   time_of_max_day=day_label(k.avg_max),
                   max_t50_cat=k.max_t50_cat)
        if args.values:
            tc = merged.gene_timecourse(g)
            tc["day"] = [day_label(p) for p in tc["position"]]
            rec["timecourse"] = tc[["dataset", "day", "position", "log2fc"]] \
                .round({"log2fc": 4}).to_dict(orient="records")
        found.append(rec)

    if args.format == "json":
        out = {"found": bool(found), "n": len(found)}
        if missing:
            out["missing"] = missing
        out["results"] = found
        print(json.dumps(out, indent=2))
    else:
        import pandas as pd
        if found:
            cols = ["gene", "t50_day", "avg_fc", "direction", "time_of_max_day", "max_t50_cat"]
            print(pd.DataFrame(found)[cols].to_string(index=False))
        if missing:
            print("missing:", ", ".join(missing), file=sys.stderr)
    sys.exit(0 if found else 2)


if __name__ == "__main__":
    main()
