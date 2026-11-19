"""NR vs DR over postnatal day: the per-sample OLS test and the time-course panel.

Test (2_development/04-05, 4_dr/01): over the P12-P21 samples (TEST_DAYS; NR and DR at each),
  y ~ age (factor) + DR + DR * (day - D0)
D0 = 16 = the mean day of the DR samples, so b_DR is the average DR effect (the effect at day 16) and
b_int is the change in the DR effect per day; t-test p-values. The same design as the per-gene DESeq2
of 2_development/03.
"""

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

from .constants import TEST_DAYS, D0, STARS, COND_COLORS


def dr_ols(t, y, suffix=''):
    """OLS of column `y` of the P12-P21 samples `t` (columns day, condition). Returns {DR_d16, se_DR_d16,
    pvalue_DR_d16, slope_per_day, se_slope, pvalue_slope}, each key + suffix."""
    t = t.assign(age=pd.Categorical(t['day'].map(lambda d: f'P{d}'), categories=[f'P{d}' for d in TEST_DAYS]),
                 dr=(t['condition'] == 'DR').astype(float))
    if t.loc[t['dr'] == 1, 'day'].mean() != D0:
        raise ValueError(f'mean DR day is {t.loc[t["dr"] == 1, "day"].mean()}, D0 = {D0}')
    t['dr_dayc'] = t['dr'] * (t['day'] - D0)
    fit = smf.ols(f'{y} ~ age + dr + dr_dayc', data=t).fit()
    return {f'DR_d16{suffix}': fit.params['dr'], f'se_DR_d16{suffix}': fit.bse['dr'],
            f'pvalue_DR_d16{suffix}': fit.pvalues['dr'], f'slope_per_day{suffix}': fit.params['dr_dayc'],
            f'se_slope{suffix}': fit.bse['dr_dayc'], f'pvalue_slope{suffix}': fit.pvalues['dr_dayc']}


def bh(p):
    """Benjamini-Hochberg adjusted p-values."""
    return multipletests(p, method='fdr_bh')[1]


def stars(padj):
    return next((s for cut, s in STARS if padj < cut), '')   # NaN padj -> ''


def bracket(ax, x, y_nr, y_dr, labels, fontsize=7, lw=1):
    """Bracket right of the last day, from x: short ticks from the NR and DR line ends out to a vertical
    bar, labelled with the significant tests."""
    ax.plot([x, x + 0.3, x + 0.3, x], [y_nr, y_nr, y_dr, y_dr], color='black', lw=lw)
    ax.text(x + 0.6, (y_nr + y_dr) / 2, '\n'.join(labels), fontsize=fontsize, va='center', ha='left')


def draw_nr_dr(ax, s, y, days, labels):
    """One time-course panel: dashed line at 0; per condition, a line through the per-day mean of column
    `y` and the samples as dots; a bracket right of P21 when `labels` (significant tests) is not empty."""
    ax.axhline(0, color='gray', lw=0.8, ls='--')
    ends = {}
    for c, color in COND_COLORS.items():
        sc = s[s['condition'] == c]
        t = sc.groupby('day')[y].mean()
        ax.plot(t.index, t.values, color=color, lw=1.5, label=c)
        ax.scatter(sc['day'], sc[y], color=color, s=10, zorder=3)
        ends[c] = t.loc[t.index.max()]
    if labels:
        bracket(ax, 21.8, ends['NR'], ends['DR'], labels)
    ax.set_xlim(5, 25.5)
    ax.set_xticks(days)
    ax.tick_params(axis='x', labelsize=6.5)
    ax.spines[['top', 'right']].set_visible(False)
