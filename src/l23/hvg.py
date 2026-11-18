"""Highly variable gene selection."""

import numpy as np


def select_hvg(x_raw, depths_raw, n_top_genes):
    """Select top highly variable genes via normalize→log1p→variance ranking."""
    median_depth = np.median(depths_raw)
    x_norm = x_raw / depths_raw[:, None] * median_depth
    gene_var = np.log1p(x_norm).var(axis=0)
    top_idx = np.argsort(gene_var)[::-1][:n_top_genes]
    mask = np.zeros(x_raw.shape[1], dtype=bool)
    mask[top_idx] = True
    return mask
