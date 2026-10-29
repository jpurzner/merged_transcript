"""Optional normalization helpers (apply BEFORE building a Dataset).

The merge operates on normalized expression. For RNA-seq counts use
`size_factor_normalize` (DESeq2 median-of-ratios); for log2 microarray values
use `unlog2`. You can also supply already-normalized matrices and skip these.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def size_factor_normalize(counts: pd.DataFrame) -> pd.DataFrame:
    """DESeq2-style median-of-ratios normalization of a counts matrix
    (genes x samples). Divides each sample by its size factor."""
    x = counts.to_numpy(float)
    log_geomean = np.mean(np.log(np.where(x > 0, x, np.nan)), axis=1)
    finite = np.isfinite(log_geomean)
    size_factors = np.empty(x.shape[1])
    for j in range(x.shape[1]):
        ratios = np.log(x[finite, j]) - log_geomean[finite]
        ratios = ratios[np.isfinite(ratios)]
        size_factors[j] = np.exp(np.median(ratios))
    return pd.DataFrame(x / size_factors, index=counts.index, columns=counts.columns)


def unlog2(values: pd.DataFrame) -> pd.DataFrame:
    """2 ** values — put log2 (e.g. RMA microarray) data on a linear scale."""
    return 2.0 ** values
