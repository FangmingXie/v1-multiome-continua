"""Dot matrix of 03's L2/3 enrichment: enriched +/+ TFs x A'/B'/C' (developmental archetype genes).

Rows are 03's selected TFs (star rule in >= 1 archetype), grouped by peak archetype and ordered
by log2_enr within each group. Encoding (src/l23/regulon_plot.py):
  colour  = log2_enr
  area    = overlap / n_markers
  gray    = overlap < 5 (too few shared genes to trust)
  outline = star rule (fdr_strat < 0.05, log2_enr > 1, overlap >= 5)

Reads:
  results/3_regulon/03.l23_stratified_enrichment.tsv
  results/3_regulon/03.enriched_regulon_selection.tsv
Outputs:
  figures/3_regulon/07.regulon_dotplot.pdf
"""

import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

from l23.regulon_plot import (SIGN, STAR_FDR, STAR_LOG2ENR, MASK_MIN_OVERLAP, MASK_COLOR, COLOR_MIN,
                              COLOR_MAX, FRAC_REF, SIZE_REF, SIZE_LEGEND, BOX_LW, log2enr_cmap)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_ENRICHMENT = os.path.join(REPO_ROOT, 'results', '3_regulon', '03.l23_stratified_enrichment.tsv')
INPUT_SELECTION = os.path.join(REPO_ROOT, 'results', '3_regulon', '03.enriched_regulon_selection.tsv')
OUT_PDF = os.path.join(REPO_ROOT, 'figures', '3_regulon', '07.regulon_dotplot.pdf')

COLS = ["A'", "B'", "C'"]
CELL_W, CELL_H = 0.50, 0.21      # inches per matrix cell

os.makedirs(os.path.dirname(OUT_PDF), exist_ok=True)


def main():
    enr = pd.read_csv(INPUT_ENRICHMENT, sep='\t')
    enr = enr[enr['regulation_direction'] == SIGN]
    sel = pd.read_csv(INPUT_SELECTION, sep='\t')
    rows = list(sel['TF'])           # already grouped by peak archetype, log2_enr descending
    if set(enr['arch_letter']) != set(COLS):
        raise ValueError(f'unexpected archetypes: {sorted(set(enr["arch_letter"]))}')

    mats = {k: enr.pivot_table(index='TF', columns='arch_letter', values=k, aggfunc='first')
                  .reindex(index=rows, columns=COLS).values
            for k in ['log2_enr', 'fdr_strat', 'overlap', 'n_markers']}
    log2, overlap = mats['log2_enr'], mats['overlap']
    if not np.isfinite(log2).all():
        raise ValueError('a selected TF is missing an archetype in the enrichment table')
    thin = overlap < MASK_MIN_OVERLAP
    sig = ((mats['fdr_strat'] < STAR_FDR) & (log2 > STAR_LOG2ENR)
           & (overlap >= MASK_MIN_OVERLAP))
    star_rows = {tf for tf, s in zip(rows, sig) if s.any()}
    if star_rows != set(rows):
        raise ValueError(f'selection and star rule disagree: {sorted(set(rows) ^ star_rows)}')
    frac = overlap / mats['n_markers']
    if frac.max() > FRAC_REF + 1e-9:
        raise ValueError(f'coverage {frac.max():.3f} exceeds the size reference {FRAC_REF}')

    cmap = log2enr_cmap()
    norm = Normalize(vmin=COLOR_MIN, vmax=COLOR_MAX)
    fig, ax = plt.subplots(layout='constrained',
                           figsize=(CELL_W * len(COLS) + 3.0, CELL_H * len(rows) + 1.6))

    yy, xx = np.indices(log2.shape).reshape(2, -1)
    sizes = frac[yy, xx] / FRAC_REF * SIZE_REF
    mappable = None
    for keep, kw in [(thin[yy, xx], dict(color=MASK_COLOR)),
                     (~thin[yy, xx], dict(c=log2[yy, xx][~thin[yy, xx]], cmap=cmap, norm=norm))]:
        if not keep.any():
            continue
        s = sig[yy, xx][keep]
        sc = ax.scatter(xx[keep], yy[keep], s=sizes[keep], zorder=3,
                        edgecolors=np.where(s, 'black', 'none'),
                        linewidths=np.where(s, BOX_LW, 0.0), **kw)
        if 'c' in kw:
            mappable = sc

    # separators between peak-archetype groups
    peak = sel['arch_letter'].values
    for i in np.nonzero(peak[1:] != peak[:-1])[0]:
        ax.axhline(i + 0.5, color='black', lw=0.8, zorder=2)

    ax.set_xticks(range(len(COLS)))
    ax.set_xticklabels(COLS, fontsize=9)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontsize=7)
    ax.set_xlim(-0.6, len(COLS) - 0.4)
    ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_axisbelow(True)
    ax.grid(True, color='0.90', lw=0.5, zorder=0)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_edgecolor('0.6')

    cbar = fig.colorbar(mappable, ax=ax, fraction=0.06, pad=0.04, aspect=30)
    cbar.set_label('log2 enrichment\n(obs / expression-matched exp)', fontsize=7)
    cbar.ax.tick_params(labelsize=6)
    cbar.ax.axhline(0.0, color='black', lw=0.8)

    handles = [Line2D([], [], marker='o', linestyle='none', markerfacecolor='0.55',
                      markeredgecolor='none',
                      markersize=np.sqrt(f / FRAC_REF * SIZE_REF), label=f'{f:.2f}')
               for f in SIZE_LEGEND]
    handles += [Line2D([], [], marker='o', linestyle='none', markerfacecolor=MASK_COLOR,
                       markeredgecolor='none', markersize=7,
                       label=f'overlap<{MASK_MIN_OVERLAP}'),
                Line2D([], [], marker='o', linestyle='none', markerfacecolor='white',
                       markeredgecolor='black', markeredgewidth=BOX_LW, markersize=7,
                       label='significant')]
    fig.legend(handles=handles, title='overlap / n_markers', loc='outside right lower',
               frameon=False, fontsize=6, title_fontsize=7, labelspacing=1.0)

    n_mk = dict(zip(COLS, mats['n_markers'][0].astype(int)))
    ax.set_title(f'Yoo25 L2/3 {SIGN} regulons enriched in L2/3 archetypes\n'
                 f'NR-union archetype genes ({", ".join(f"{a} {n}" for a, n in n_mk.items())})\n'
                 f'outlined: FDR<{STAR_FDR:g}, log2 enr>{STAR_LOG2ENR:g}, '
                 f'overlap>={MASK_MIN_OVERLAP}', fontsize=7)

    fig.savefig(OUT_PDF, bbox_inches='tight')
    plt.close(fig)
    print(f'{len(rows)} TFs x {len(COLS)} archetypes, {int(sig.sum())} significant cells, '
          f'{int(thin.sum())} masked; saved {OUT_PDF}')


if __name__ == '__main__':
    main()
