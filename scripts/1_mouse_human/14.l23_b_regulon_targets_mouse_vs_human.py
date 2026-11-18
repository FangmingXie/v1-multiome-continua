"""B'-enriched L2/3 regulons: target-gene activity over normal development, mouse vs human.

Mouse counterpart of 13 (human Wang25 AP-1 / EGR regulon targets over stage), same score and
style, shown side by side with it. NR vs DR in mouse is 4_dr/02.

TFs: the AP-1 + EGR members of 3_regulon/04's B' TFs (developmental archetype genes, star rule)
that have a +/+ regulon in both species (human has no FOSB regulon). One colour per TF in both
panels.

Mouse score (as 13): Yoo25 L2/3 cells of the normally reared (NR) ages P6-P21; each +/+ target
(3_regulon/01) log2(CP10k+1) divided by its NORM_MAX_PCTILE percentile over those cells, capped
at 1, averaged over the regulon's targets (l23.regulon.load_expression / regulon_score). Per age:
mean +/- SEM over cells, divided by the largest per-age mean (1 = peak). Human: 13's table,
unchanged.

Reads:
  data/raw/superdupermegaRNA_hasraw.h5ad
  results/3_regulon/01.yoo25_l23_regulon_targets.tsv
  results/3_regulon/04.enriched_regulon_selection.tsv
  results/1_mouse_human/13.wang25_l23_ap1_egr_regulon_targets_over_stage.tsv
Outputs:
  figures/1_mouse_human/14.l23_b_regulon_targets_mouse_vs_human.pdf
  results/1_mouse_human/14.l23_b_regulon_targets_nr_mouse.tsv
"""

import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt

from l23.regulon import load_expression, regulon_score

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
INPUT_SOURCE = os.path.join(REPO_ROOT, 'data', 'raw', 'superdupermegaRNA_hasraw.h5ad')
INPUT_REGULONS = os.path.join(REPO_ROOT, 'results', '3_regulon', '01.yoo25_l23_regulon_targets.tsv')
INPUT_SELECTION = os.path.join(REPO_ROOT, 'results', '3_regulon', '04.enriched_regulon_selection.tsv')
INPUT_HUMAN = os.path.join(REPO_ROOT, 'results', '1_mouse_human', '13.wang25_l23_ap1_egr_regulon_targets_over_stage.tsv')
OUT_PDF = os.path.join(REPO_ROOT, 'figures', '1_mouse_human', '14.l23_b_regulon_targets_mouse_vs_human.pdf')
OUT_TSV = os.path.join(REPO_ROOT, 'results', '1_mouse_human', '14.l23_b_regulon_targets_nr_mouse.tsv')

ARCH = "B'"
SIGN = '+/+'
# B' TFs by family; every B' TF of 3_regulon/04 must land in exactly one
FAMILIES = [('AP-1', ['Fos', 'Fosb', 'Fosl2', 'Junb']),
            ('EGR', ['Egr1', 'Egr2', 'Egr3', 'Egr4']),
            ('other', ['Smad3', 'Irf2', 'Mxi1', 'Etv5', 'Tfdp2'])]
HUMAN_FAMILIES = ['AP-1', 'EGR']
HUMAN_STAGES = {'Second_trimester': '2nd trim.', 'Third_trimester': '3rd trim.',
                'Infancy': 'Infancy', 'Adolescence': 'Adolescence'}
# colour: TF i (AP-1 then EGR order) = tab20 colour i; marker = member order within its family
CMAP = 'tab20'
LINE_MARKERS = ['o', 's', '^', 'D']

os.makedirs(os.path.dirname(OUT_PDF), exist_ok=True)
os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)


def load_targets():
    """{TF: sorted +/+ targets} for every B' TF of 3_regulon/04, in FAMILIES order."""
    sel = pd.read_csv(INPUT_SELECTION, sep='\t')
    tfs = set(sel.loc[sel['arch_letter'] == ARCH, 'TF'])
    fam_tfs = [tf for _, fam in FAMILIES for tf in fam]
    if sorted(fam_tfs) != sorted(tfs):
        raise ValueError(f'FAMILIES {sorted(fam_tfs)} != {ARCH} TFs of {INPUT_SELECTION} {sorted(tfs)}')
    reg = pd.read_csv(INPUT_REGULONS, sep='\t')
    reg = reg[reg['regulation_direction'] == SIGN]
    targets = {tf: sorted(set(reg.loc[reg['TF'] == tf, 'Gene'])) for tf in fam_tfs}
    empty = [tf for tf, t in targets.items() if not t]
    if empty:
        raise ValueError(f'no {SIGN} regulon in {INPUT_REGULONS} for {empty}')
    return targets


def draw(ax, s, color, marker, label):
    ax.errorbar(s['x'], s['mean_norm'], yerr=s['sem_norm'], color=color, marker=marker,
                markersize=4.5, capsize=2.0, linewidth=1.3, elinewidth=0.9, label=label)


def main():
    hum = pd.read_csv(INPUT_HUMAN, sep='\t')
    if set(hum['stage']) != set(HUMAN_STAGES):
        raise ValueError(f'unexpected human stages: {sorted(set(hum["stage"]))}')
    targets = load_targets()
    mouse_tfs = [tf for fam, tfs in FAMILIES if fam in HUMAN_FAMILIES for tf in tfs]
    tfs = [tf for tf in mouse_tfs if tf.upper() in set(hum['TF'])]
    print(f'mouse vs human: {tfs}; mouse-only (no human regulon): '
          f'{[tf for tf in mouse_tfs if tf not in tfs]}')

    L, cells = load_expression(INPUT_SOURCE, sorted({g for tf in tfs for g in targets[tf]}))
    nr = (cells['condition'] == 'NR').to_numpy()
    L, cells = L[nr].reset_index(drop=True), cells[nr].reset_index(drop=True)

    rows = []
    for tf in tfs:
        score, n_used = regulon_score(tf, targets[tf], L)
        st = (cells.assign(score=score).groupby(['Age', 'day'])['score']
                   .agg(n_cells='size', mean='mean', std='std').reset_index().sort_values('day'))
        st['sem'] = st['std'] / np.sqrt(st['n_cells'])
        peak = st['mean'].max()
        if not peak > 0:
            raise ValueError(f'{tf}: peak mean is {peak}; cannot max-scale')
        st['mean_norm'], st['sem_norm'] = st['mean'] / peak, st['sem'] / peak
        st.insert(0, 'TF', tf)
        st.insert(1, 'n_targets_used', n_used)
        rows.append(st)
    mouse = pd.concat(rows, ignore_index=True)
    mouse.to_csv(OUT_TSV, sep='\t', index=False)
    print(f'wrote {OUT_TSV}')

    days = sorted(mouse['day'].unique())
    stage_x = {s: i for i, s in enumerate(HUMAN_STAGES)}
    fig, (axm, axh) = plt.subplots(1, 2, figsize=(6.0, 3.6), sharey=True)
    style = {}
    for fam, fam_tfs in FAMILIES:
        if fam in HUMAN_FAMILIES:
            for j, tf in enumerate(t for t in fam_tfs if t in tfs):
                style[tf] = (plt.get_cmap(CMAP)(len(style)), LINE_MARKERS[j])
    for i, tf in enumerate(tfs):
        color, marker = style[tf]
        m = mouse[mouse['TF'] == tf].assign(x=lambda d: d['day'])
        draw(axm, m, color, marker, None)
        h = hum[hum['TF'] == tf.upper()].assign(x=lambda d: d['stage'].map(stage_x)).sort_values('x')
        draw(axh, h, color, marker,
             f'{tf} ({m["n_targets_used"].iloc[0]} / {h["n_targets_used"].iloc[0]})')
    axm.set_xticks(days)
    axm.set_xticklabels([f'P{d}' for d in days], fontsize=7)
    axm.set_xlim(days[0] - 1, days[-1] + 1)
    axm.set_xlabel('postnatal age')
    axm.set_title('mouse V1 L2/3 (Yoo25, NR)', fontsize=9)
    axm.set_ylabel('regulon target activity\n(per-gene max-scaled mean, /peak)', fontsize=8)
    axh.set_xticks(range(len(HUMAN_STAGES)))
    axh.set_xticklabels(list(HUMAN_STAGES.values()), rotation=45, ha='right', fontsize=7)
    axh.set_xlim(-0.4, len(HUMAN_STAGES) - 0.6)
    axh.set_xlabel('developmental stage')
    axh.set_title('human V1 L2/3 (Wang25)', fontsize=9)
    for ax in (axm, axh):
        ax.set_ylim(bottom=0)
        ax.axhline(1.0, color='0.8', lw=0.6, zorder=0)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    fig.suptitle(f'AP-1 & EGR {SIGN} regulon target activity over development, mouse vs human\n'
                 'mean ± SEM over cells', fontsize=10)
    axh.legend(fontsize=7, frameon=False, title='TF (n targets:\nmouse / human)', title_fontsize=7,
               loc='center left', bbox_to_anchor=(1.02, 0.5))
    fig.tight_layout()
    fig.savefig(OUT_PDF, bbox_inches='tight')
    plt.close(fig)
    print(f'wrote {OUT_PDF}')


if __name__ == '__main__':
    main()
