"""Archetype-gene overlap of 03's enriched L2/3 regulons: one hbar panel per archetype (A', B', C').

Each enriched +/+ regulon of 03's selection (NR-union archetype genes, star rule) is drawn in the
panel of the archetype it is enriched in. Its bar stacks the number of A', B' and C' marker genes
among the regulon's targets (03's `overlap` per archetype; the marker sets are disjoint), with the
enriched archetype's segment first so it starts at 0. Bars are sorted by the enriched archetype's overlap, largest on top. Every
selected TF is starred in exactly one archetype, so each appears once.

Reads:
  results/3_regulon/03.enriched_regulon_selection.tsv
  results/3_regulon/03.l23_stratified_enrichment.tsv   (overlap with every archetype)
Outputs:
  figures/3_regulon/08.regulon_overlap_bar.pdf
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from l23.regulon_plot import ARCH_COLORS

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# --- file paths ---
INPUT_SELECTION = os.path.join(REPO_ROOT, 'results', '3_regulon', '03.enriched_regulon_selection.tsv')
INPUT_ENRICHMENT = os.path.join(REPO_ROOT, 'results', '3_regulon', '03.l23_stratified_enrichment.tsv')
OUT_PDF = os.path.join(REPO_ROOT, 'figures', '3_regulon', '08.regulon_overlap_bar.pdf')

ROW_H = 0.22                                         # inches per bar

os.makedirs(os.path.dirname(OUT_PDF), exist_ok=True)

sel = pd.read_csv(INPUT_SELECTION, sep='\t')
if set(sel['arch_letter']) != set(ARCH_COLORS):
    raise ValueError(f'unexpected archetypes: {sorted(set(sel["arch_letter"]))}')
if sel['TF'].duplicated().any():
    raise ValueError(f'TF selected in >1 archetype: {sorted(sel.loc[sel["TF"].duplicated(), "TF"])}')
enr = pd.read_csv(INPUT_ENRICHMENT, sep='\t')
overlap = (enr[enr['regulon'].isin(sel['regulon'])]
           .pivot(index='regulon', columns='arch_letter', values='overlap')
           .reindex(columns=list(ARCH_COLORS)))
if overlap.isna().any().any() or len(overlap) != len(sel):
    raise ValueError('a selected regulon is missing an archetype in the enrichment table')
sel = sel.join(overlap.add_prefix('ov_'), on='regulon')
if not (sel['overlap'] == sel.apply(lambda r: r[f"ov_{r['arch_letter']}"], axis=1)).all():
    raise ValueError('selection overlap disagrees with the enrichment table')
sel['ov_total'] = sel[[f'ov_{a}' for a in ARCH_COLORS]].sum(axis=1)

n_max = sel.groupby('arch_letter').size().max()
plt.rcParams['pdf.fonttype'] = 42
fig, axes = plt.subplots(1, len(ARCH_COLORS), figsize=(3.0 * len(ARCH_COLORS), ROW_H * n_max + 1.0),
                         layout='constrained')
for ax, (arch, color) in zip(axes, ARCH_COLORS.items()):
    sub = sel[sel['arch_letter'] == arch].sort_values('overlap', ascending=True)
    left = 0
    for a in [arch] + [a for a in ARCH_COLORS if a != arch]:       # enriched archetype first
        ax.barh(sub['TF'], sub[f'ov_{a}'], left=left, color=ARCH_COLORS[a], height=0.75,
                label=f'{a} genes')
        left = left + sub[f'ov_{a}'].values
    ax.set_ylim(len(sub) - n_max - 0.6, len(sub) - 0.4)   # same bar height in every panel; top-aligned
    ax.set_xlim(0, sub['ov_total'].max() * 1.05)
    ax.set_title(f"{arch} regulons (n = {len(sub)}; {sub['n_markers'].iloc[0]} {arch} genes)",
                 fontsize=9, color=color)
    ax.tick_params(axis='y', labelsize=8)
    sns.despine(ax=ax)
fig.supxlabel('# L2/3 archetype genes among regulon targets', fontsize=9)
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in ARCH_COLORS.values()]
axes[-1].legend(handles, [f'{a} genes' for a in ARCH_COLORS], loc='lower right', fontsize=8,
                frameon=False)                       # the short C' panel has room below its bars
fig.suptitle('Enriched L2/3 regulons: archetype genes among targets', fontsize=10)
fig.savefig(OUT_PDF, bbox_inches='tight', dpi=300)
plt.close(fig)
print(sel[['TF', 'arch_letter'] + [f'ov_{a}' for a in ARCH_COLORS] + ['ov_total', 'n_targets']].to_string(index=False))
print(f'Saved {OUT_PDF}')
