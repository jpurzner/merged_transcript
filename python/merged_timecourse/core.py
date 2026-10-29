"""Core algorithm: normalize -> interpolate -> align (N datasets) -> kinetics.

A Python recreation of the R merged_transcript pipeline, generalized to an
arbitrary number of developmental time-course datasets. The alignment uses
sequential anchoring: the first dataset is the reference, and each subsequent
dataset is shifted (per gene) to agree with the datasets already placed on the
shared pseudo-time axis. This is what makes both merging many datasets and
incrementally *adding* a dataset to an existing merge one and the same routine
(see merge.py).

The numeric recipe mirrors R/lib/{time_expr_interp,time_expr_align,
time_expr_avg,expr_half_max}.R; for two datasets it reproduces the R output
(validated in tests/).
"""
from __future__ import annotations

import warnings
import numpy as np
import pandas as pd


# --- normalization --------------------------------------------------------

def row_log2fc(matrices):
    """Row-normalize each (genes x samples) matrix to log2 fold change vs that
    gene's own mean. Zeros across all matrices are first floored to the global
    smallest positive value (matching the R pipeline). Returns a list of arrays.
    """
    pos_min = min(np.nanmin(m[m > 0]) for m in matrices)
    out = []
    for m in matrices:
        m = np.where(m == 0, pos_min, m)
        row_mean = m.mean(axis=1, keepdims=True)
        out.append(np.log2(m / row_mean))
    return out


def average_by_timepoint(log2fc, sample_times, time_order):
    """Mean across replicate columns within each timepoint.
    log2fc: (genes x samples) array; sample_times: length-samples labels;
    time_order: ordered unique timepoint labels. Returns (genes x timepoints).
    """
    sample_times = np.asarray(sample_times)
    cols = [log2fc[:, sample_times == t].mean(axis=1) for t in time_order]
    return np.column_stack(cols)


# --- interpolation --------------------------------------------------------

def interpolate(avg, positions, grid):
    """Linear-interpolate each gene's per-timepoint values onto `grid`.
    avg: (genes x timepoints); positions: numeric axis position of each
    timepoint; grid: dense numeric axis. Points outside a dataset's own
    sampled range become NaN (matching R approx). Returns (genes x len(grid)).
    """
    positions = np.asarray(positions, float)
    order = np.argsort(positions)
    positions = positions[order]
    avg = avg[:, order]
    lo, hi = positions[0], positions[-1]
    out = np.empty((avg.shape[0], len(grid)), float)
    for i in range(avg.shape[0]):
        out[i] = np.interp(grid, positions, avg[i], left=np.nan, right=np.nan)
    out[:, (np.asarray(grid) < lo) | (np.asarray(grid) > hi)] = np.nan
    return out


# --- alignment (N datasets) ----------------------------------------------

def align(interp_list):
    """Align interpolated datasets onto a common scale.

    Sequential anchoring: dataset 0 is the reference; each later dataset is
    shifted (per gene) by the median difference, over the time points where it
    overlaps the mean of the already-placed datasets, then everything is
    centred on zero per gene. Returns (aligned_list, shifts) where shifts is
    (genes x datasets) the per-gene shift applied before centring.
    """
    n_ds = len(interp_list)
    n_genes = interp_list[0].shape[0]
    aligned = [interp_list[0].copy()]
    shifts = np.zeros((n_genes, n_ds))
    warnings.filterwarnings("ignore", r"(Mean|All-NaN)", RuntimeWarning)
    for k in range(1, n_ds):
        placed = np.stack(aligned, axis=0)             # placed x genes x axis
        ref_mean = np.nanmean(placed, axis=0)          # genes x axis
        cur = interp_list[k]
        overlap = ~np.isnan(ref_mean) & ~np.isnan(cur)  # genes x axis
        diff = np.where(overlap, ref_mean - cur, np.nan)
        shift_k = np.nanmedian(diff, axis=1)           # genes
        shift_k = np.where(np.isnan(shift_k), 0.0, shift_k)
        shifts[:, k] = shift_k
        aligned.append(cur + shift_k[:, None])

    stacked = np.stack(aligned, axis=0)                # ds x genes x axis
    all_min = np.nanmin(stacked, axis=(0, 2))          # genes
    all_max = np.nanmax(stacked, axis=(0, 2))
    centre = ((all_max - all_min) / 2.0 + all_min)[:, None]
    aligned = [a - centre for a in aligned]
    shifts = shifts - centre                           # fold the centring in
    return aligned, shifts


def average_across_datasets(aligned_list):
    """Per gene, per axis point, mean across datasets (ignoring NaN)."""
    return np.nanmean(np.stack(aligned_list, axis=0), axis=0)


# --- kinetics -------------------------------------------------------------

def _find_y0(j):
    """Fractional axis index at 50% of the max change (ports expr_half_max's
    find_y0): NaN if the gene never exceeds 0; else interpolate the midpoint on
    the monotonic segment between the min and the max."""
    j = j[~np.isnan(j)]
    if j.size == 0 or np.max(j) < 0:
        return np.nan
    imin, imax = int(np.argmin(j)), int(np.argmax(j))   # 0-based
    start, end = min(imin, imax), max(imin, imax)
    jtrim = j[start:end + 1]
    mid = (j.max() + j.min()) / 2.0
    yidx = np.arange(1, jtrim.size + 1, dtype=float)    # 1-based, like R
    order = np.argsort(jtrim)
    maxi = np.interp(mid, jtrim[order], yidx[order])
    return maxi + (start + 1)                           # + 1-based start (R)


def half_max_t50(iaa):
    """t50 per gene: _find_y0 clamped to the last non-NaN axis index (ports
    expr_half_max for the cross-dataset average)."""
    out = np.empty(iaa.shape[0])
    for i in range(iaa.shape[0]):
        row = iaa[i]
        nn = np.where(~np.isnan(row))[0]
        last = (nn[-1] + 1) if nn.size else np.nan      # 1-based last index
        y0 = _find_y0(row)
        out[i] = np.nanmin([y0, last]) if not (np.isnan(y0) and np.isnan(last)) else np.nan
    return out


def kinetics(iaa, grid):
    """Per-gene kinetic metrics from the cross-dataset aligned average.
    Returns a DataFrame with avg_fc, avg_dir, avg_max, avg_min, avg_t50 and the
    max_cat / t50_cat / max_t50_cat categories (indices are 1-based positions on
    `grid`, matching the R day axis)."""
    grid = np.asarray(grid)
    fc = np.nanmax(iaa, axis=1) - np.nanmin(iaa, axis=1)
    amax = np.nanargmax(iaa, axis=1)
    amin = np.nanargmin(iaa, axis=1)
    direction = np.where(amax > amin, 1, -1)
    # 1-based index on the grid (R which.max). For grid == 1..N this is grid[idx].
    avg_max = (amax + 1).astype(float)
    avg_min = (amin + 1).astype(float)
    t50 = half_max_t50(iaa)
    df = pd.DataFrame(dict(avg_fc=fc, avg_dir=direction,
                           avg_max=avg_max, avg_min=avg_min, avg_t50=t50))
    df["max_cat"] = pd.cut(df["avg_max"], [0, 2, 20, 64], right=False,
                           labels=["E15", "P0_7", "P14_P56"]).astype("object")
    df["t50_cat"] = pd.cut(df["avg_t50"], [0, 6, 18, 64], right=False,
                           labels=["Early", "Mid", "Late"]).astype("object")
    df["max_t50_cat"] = "Max_" + df["max_cat"].astype(str) + " t50_" + df["t50_cat"].astype(str)
    return df


_PROGRAM = {"E15": "embryonic", "P0_7": "proliferation", "P14_P56": "differentiation",
            "unchanged": "unchanged", "unexpressed": "unexpressed"}


def apply_gene_program(kin, expressed, significant=None, fc_cutoff=1.0):
    """Overlay the `unchanged` / `unexpressed` calls onto the peak-timing
    categories and add a `gene_program` column, mirroring the R pipeline
    (Purzner et al. 2025 developmental programs).

    expressed   : boolean per gene (called expressed in >=1 dataset).
    significant : optional boolean per gene (significant change in >=1 dataset);
                  if omitted, only the fold-change cutoff gates `unchanged`.
    """
    kin = kin.copy()
    expressed = np.asarray(expressed, bool)
    unchanged = kin["avg_fc"].to_numpy() < fc_cutoff
    if significant is not None:
        unchanged = unchanged | ~np.asarray(significant, bool)
    max_cat = kin["max_cat"].astype(str).to_numpy().copy()
    max_t50 = kin["max_t50_cat"].astype(str).to_numpy().copy()
    max_cat[unchanged] = "unchanged"; max_t50[unchanged] = "unchanged"
    max_cat[~expressed] = "unexpressed"; max_t50[~expressed] = "unexpressed"
    kin["max_cat"] = max_cat
    kin["max_t50_cat"] = max_t50
    kin["gene_program"] = [_PROGRAM.get(c) for c in max_cat]
    return kin
