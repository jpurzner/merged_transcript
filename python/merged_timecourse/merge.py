"""Merge N time-course datasets, and add datasets to an existing merge."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd

from .dataset import Dataset
from . import core


@dataclass
class MergedTimecourse:
    """The result of merging one or more datasets onto a shared pseudo-time axis.

    datasets : the input Dataset objects (in merge order; the first is the
               alignment reference).
    genes    : gene ids common to all datasets.
    grid     : the dense numeric pseudo-time axis (interpolation grid).
    aligned  : one (genes x len(grid)) DataFrame of aligned log2 fold change per
               dataset, keyed by dataset name.
    average  : (genes x len(grid)) cross-dataset mean (the aligned consensus).
    kinetics : per-gene t50 / fold change / direction / categories.
    shifts   : (genes x datasets) per-gene alignment shift applied.
    """
    datasets: list
    genes: pd.Index
    grid: np.ndarray
    aligned: dict
    average: pd.DataFrame
    kinetics: pd.DataFrame
    shifts: pd.DataFrame

    @property
    def dataset_names(self):
        return [d.name for d in self.datasets]

    def gene_timecourse(self, gene):
        """Aligned log2FC trajectory for one gene, long form (dataset, position,
        log2fc) across all datasets — the numbers behind a plot."""
        rows = []
        for name, mat in self.aligned.items():
            vals = mat.loc[gene]
            for pos, v in zip(self.grid, vals):
                if not np.isnan(v):
                    rows.append((name, float(pos), float(v)))
        return pd.DataFrame(rows, columns=["dataset", "position", "log2fc"])

    def add_dataset(self, new_dataset, grid=None):
        """Return a new MergedTimecourse with `new_dataset` added. The new
        dataset is anchored onto the existing scale (the original first dataset
        stays the reference)."""
        return merge_datasets(list(self.datasets) + [new_dataset], grid=grid)


def _build_grid(datasets):
    positions = [p for d in datasets for p in d.timepoint_positions]
    lo, hi = int(np.floor(min(positions))), int(np.ceil(max(positions)))
    return np.arange(lo, hi + 1, dtype=float)


def merge_datasets(datasets, grid=None):
    """Merge >=1 datasets onto a shared pseudo-time axis.

    The first dataset is the alignment reference; each subsequent dataset is
    interpolated onto `grid` and shifted to agree with the datasets already
    placed (see core.align). Genes are intersected across datasets. Returns a
    MergedTimecourse.
    """
    if len(datasets) == 0:
        raise ValueError("need at least one dataset")

    # common genes, in the first dataset's order
    genes = datasets[0].genes
    for d in datasets[1:]:
        genes = genes[genes.isin(d.genes)]
    if len(genes) == 0:
        raise ValueError("datasets share no genes")

    if grid is None:
        grid = _build_grid(datasets)

    # 1. row-normalize to log2 fold change (global zero-floor across all datasets)
    mats = [d.expr.loc[genes].to_numpy(float) for d in datasets]
    log2 = core.row_log2fc(mats)

    # 2. average per timepoint, 3. interpolate onto the grid
    interp = []
    for d, l in zip(datasets, log2):
        avg = core.average_by_timepoint(l, d.sample_times, d.time_order)
        interp.append(core.interpolate(avg, d.timepoint_positions, grid))

    # 4. align (sequential anchoring + centring), 5. cross-dataset average
    aligned_arrs, shifts = core.align(interp)
    iaa = core.average_across_datasets(aligned_arrs)

    # 6. kinetics
    kin = core.kinetics(iaa, grid)
    kin.index = genes

    aligned = {d.name: pd.DataFrame(a, index=genes, columns=grid)
               for d, a in zip(datasets, aligned_arrs)}
    return MergedTimecourse(
        datasets=list(datasets), genes=genes, grid=grid, aligned=aligned,
        average=pd.DataFrame(iaa, index=genes, columns=grid),
        kinetics=kin,
        shifts=pd.DataFrame(shifts, index=genes, columns=[d.name for d in datasets]),
    )


def add_dataset(merged, new_dataset, grid=None):
    """Add a dataset to an existing MergedTimecourse (functional form of
    MergedTimecourse.add_dataset)."""
    return merged.add_dataset(new_dataset, grid=grid)
