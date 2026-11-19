"""Shared encoding of the regulon x archetype enrichment figures (3_regulon/05, 06).

  colour  = log2_enr on a blue (below 0) / YlOrRd (above 0) ramp with a hard break at 0,
            the value where the observed overlap equals the expression-matched expectation
  area    = overlap / n_markers (FRAC_REF drawn at SIZE_REF)
  gray    = overlap < MASK_MIN_OVERLAP (too few shared genes to trust)
  outline = star rule (fdr_strat < STAR_FDR, log2_enr > STAR_LOG2ENR, overlap >= MASK_MIN_OVERLAP)
"""

import re

import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from plotly.colors import sample_colorscale

from .constants import ARCH_COLORS  # noqa: F401  (archetype colours of the figures)

SIGN = '+/+'                 # activating regulons only
STAR_FDR = 0.05
STAR_LOG2ENR = 1.0
MASK_MIN_OVERLAP = 5
MASK_COLOR = '#d9d9d9'
COLOR_MIN, COLOR_MAX = -1.5, 3.5
FRAC_REF = 0.40
SIZE_REF = 260.0
SIZE_LEGEND = [0.05, 0.15, 0.25, 0.35]
BOX_LW = 1.4                 # outline width for significant cells


def log2enr_cmap():
    """log2_enr colormap: YlOrRd above zero, a blue tail below it, the break exactly at zero.

    The plotly scale repeats a position (the hard step at 0); matplotlib requires strictly
    increasing stops, so the duplicate is nudged by one part in a million.
    """
    zero = (0.0 - COLOR_MIN) / (COLOR_MAX - COLOR_MIN)
    neg = sample_colorscale('Blues', np.linspace(0.62, 0.04, 12))
    pos = sample_colorscale('YlOrRd', np.linspace(0, 0.74, 21))
    scale = ([[zero * t, c] for t, c in zip(np.linspace(0, 1, 12), neg)]
             + [[zero + (1 - zero) * t, c] for t, c in zip(np.linspace(0, 1, 21), pos)])
    stops, eps = [], 1e-6
    for pos_, css in scale:
        rgb = tuple(int(v) / 255 for v in re.findall(r'\d+', css)[:3])
        if stops and pos_ <= stops[-1][0]:
            pos_ = stops[-1][0] + eps
        stops.append((pos_, rgb))
    lo, hi = stops[0][0], stops[-1][0]
    stops = [((p - lo) / (hi - lo), c) for p, c in stops]
    return LinearSegmentedColormap.from_list('log2enr', stops)
