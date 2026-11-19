"""Expression transforms and gene-set scores shared by the steps."""

import numpy as np


def log2_cp10k(x, lib):
    """log2(CP10k + 1) of a cells x genes count matrix, given each cell's library size."""
    return np.log2(x / lib[:, None] * 1e4 + 1)


def pctile_score(x, pctile):
    """Gene-set score per cell: each gene divided by its `pctile` percentile over the cells, capped at 1,
    mean over the genes (cells x genes in). 0 = no expression; a gene zero in more than (100 - pctile)%
    of cells keeps its raw value, capped at 1."""
    hi = np.percentile(x, pctile, axis=0)
    hi = np.where(hi > 0, hi, 1.0)
    return np.clip(x / hi, 0.0, 1.0).mean(axis=1)


def minmax_score(x, lo_pctile, hi_pctile):
    """Gene-set score per cell: each gene min-max scaled between its `lo_pctile` and `hi_pctile`
    percentiles over the cells, clipped to [0, 1], mean over the genes (cells x genes in)."""
    lo = np.percentile(x, lo_pctile, axis=0)
    hi = np.percentile(x, hi_pctile, axis=0)
    rng = np.where(hi > lo, hi - lo, 1.0)
    return np.clip((x - lo) / rng, 0, 1).mean(axis=1)
