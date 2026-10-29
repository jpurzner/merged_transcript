"""Timecourse plotting (matplotlib) — the Python counterpart of the R figures."""
from __future__ import annotations

import numpy as np


def plot_gene(merged, gene, ax=None, labels=None):
    """Plot a gene's aligned log2FC trajectory across all datasets.

    merged : MergedTimecourse; gene : gene id in merged.genes.
    labels : optional {position: tick_label} for the x-axis.
    Returns the matplotlib Axes.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(5, 4))
    for name, mat in merged.aligned.items():
        y = mat.loc[gene].to_numpy(float)
        x = np.asarray(merged.grid, float)
        m = ~np.isnan(y)
        ax.plot(x[m], y[m], marker="o", ms=3, label=name)
    if gene in merged.kinetics.index:
        t50 = merged.kinetics.loc[gene, "avg_t50"]
        if np.isfinite(t50):
            ax.axvline(t50, color="k", lw=0.8, ls="-")
            ax.text(t50, ax.get_ylim()[0], " t50", va="bottom", fontsize=9)
    ax.axhline(0, color="grey", lw=0.5)
    ax.set_xlabel("developmental time (axis position)")
    ax.set_ylabel("log2 fold change")
    kc = merged.kinetics.loc[gene, "max_t50_cat"] if gene in merged.kinetics.index else ""
    ax.set_title(f"{gene}  {kc}")
    if labels:
        pos = sorted(labels)
        ax.set_xticks(pos); ax.set_xticklabels([labels[p] for p in pos])
    ax.legend(frameon=False, fontsize=8)
    ax.figure.patch.set_facecolor("white")
    return ax
