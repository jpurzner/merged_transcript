"""Dataset: one developmental time-course on the shared pseudo-time axis."""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np
import pandas as pd


@dataclass
class Dataset:
    """A single time-course dataset.

    name          : identifier, used as the column label in merged output.
    expr          : genes (index) x samples (columns), NORMALIZED expression
                    (e.g. DESeq2-normalized counts, or 2**rma array intensity).
                    See normalize.py for helpers.
    sample_times  : timepoint label for each sample column (len == n columns).
    positions     : {timepoint_label: numeric position on the shared axis}, e.g.
                    developmental day. Defines where this dataset sits in time
                    and drives interpolation/alignment against other datasets.
    """
    name: str
    expr: pd.DataFrame
    sample_times: list
    positions: dict
    time_order: list = field(default_factory=list)

    def __post_init__(self):
        self.expr = self.expr.copy()
        self.sample_times = list(self.sample_times)
        if len(self.sample_times) != self.expr.shape[1]:
            raise ValueError(
                f"{self.name}: {len(self.sample_times)} sample_times for "
                f"{self.expr.shape[1]} expression columns")
        missing = set(self.sample_times) - set(self.positions)
        if missing:
            raise ValueError(f"{self.name}: no axis position for timepoint(s) {sorted(missing)}")
        # order timepoints by their numeric position
        uniq = list(dict.fromkeys(self.sample_times))
        self.time_order = sorted(uniq, key=lambda t: self.positions[t])

    @property
    def genes(self) -> pd.Index:
        return self.expr.index

    @property
    def timepoint_positions(self) -> list:
        return [self.positions[t] for t in self.time_order]

    @classmethod
    def from_files(cls, name, expr_tsv, meta_tsv, positions, gene_col=None):
        """Load from a (genes x samples) TSV + a sample->timepoint TSV
        (columns `sample`, `timepoint`)."""
        expr = pd.read_csv(expr_tsv, sep="\t")
        gene_col = gene_col or expr.columns[0]
        expr = expr.set_index(gene_col)
        meta = pd.read_csv(meta_tsv, sep="\t")
        meta = meta.set_index("sample").loc[expr.columns]
        return cls(name=name, expr=expr, sample_times=meta["timepoint"].tolist(),
                   positions=positions)
