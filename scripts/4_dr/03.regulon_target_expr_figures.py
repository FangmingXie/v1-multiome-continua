"""Figures of 01: archetype-enriched Yoo25 L2/3 regulons, target-gene activity over postnatal age, NR vs DR.

  lines:   the top TOP_N TFs per archetype by log2_enr (one row per archetype; C has only 5), one
           panel per TF, NR (C1) vs DR (black), own y-axis; dots = replicates, line = replicate mean;
           title: OLS DR effect at day 16 and DR slope per day (padj); bracket right of P21: padj < FDR.
  barplot: every TF, height = rel_change_P21 = (P21 DR - P21 NR) / P21 NR of the replicate-mean score;
           dots = each DR replicate against the P21 NR mean. Bars coloured by archetype, ordered by value
           (largest decrease first); stars = DR effect padj.

Reads:
  results/4_dr/01.regulon_target_expr_nr_vs_dr_replicates.tsv
  results/4_dr/01.regulon_target_expr_nr_vs_dr.tsv
  results/4_dr/01.regulon_target_expr_nr_vs_dr_stats.tsv
Outputs:
  figures/4_dr/03.regulon_activity_nr_vs_dr.pdf
  figures/4_dr/03.regulon_activity_p21_dr_change.pdf
"""

import os

import pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt

from l23.constants import ARCH_COLORS, COND_COLORS, FDR, TEST_DAYS
from l23.timecourse import stars, bracket

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES_DIR = os.path.join(REPO_ROOT, 'results', '4_dr')
FIG_DIR = os.path.join(REPO_ROOT, 'figures', '4_dr')
IN_REP = os.path.join(RES_DIR, '01.regulon_target_expr_nr_vs_dr_replicates.tsv')
IN_AGE = os.path.join(RES_DIR, '01.regulon_target_expr_nr_vs_dr.tsv')
IN_STATS = os.path.join(RES_DIR, '01.regulon_target_expr_nr_vs_dr_stats.tsv')
OUT_PDF = os.path.join(FIG_DIR, '03.regulon_activity_nr_vs_dr.pdf')
OUT_PDF_BAR = os.path.join(FIG_DIR, '03.regulon_activity_p21_dr_change.pdf')

ARCHS = list(ARCH_COLORS)
TOP_N = 8                      # TFs per archetype row = panels per row
SIGN = '+/+'
CONDITIONS = ['NR', 'DR']
REP_DOT_SIZE = 9
BAR_DAY = 21

os.makedirs(FIG_DIR, exist_ok=True)


def mark(ax, age, r):
    """Bracket right of the last day, labelled with the significant terms."""
    lab = [f'{name} {stars(r[f"padj_{key}"])}' for name, key in [('DR', 'DR_d16'), ('slope', 'slope')]
           if stars(r[f'padj_{key}'])]
    if not lab:
        return
    last = max(TEST_DAYS)
    ends = {c: age.loc[(age['condition'] == c) & (age['day'] == last), 'mean_norm'].iloc[0]
            for c in CONDITIONS}
    bracket(ax, last + 0.8, ends['NR'], ends['DR'], lab, fontsize=6, lw=0.8)


def draw(ax, rep, age, color, label):
    """Replicate dots + line through the replicate means."""
    ax.plot(age['day'], age['mean_norm'], color=color, linewidth=1.3, label=label)
    ax.scatter(rep['day'], rep['mean_norm'], s=REP_DOT_SIZE, color=color, linewidths=0, zorder=3)


def plot(rows, rep, age, st, days):
    fig, axes = plt.subplots(len(ARCHS), TOP_N, figsize=(1.9 * TOP_N, 2.4 * len(ARCHS)),
                             sharex=True, squeeze=False)
    rows = {a: tfs[:TOP_N] for a, tfs in rows.items()}
    for r, a in enumerate(ARCHS):
        tfs = rows[a]
        for c in range(TOP_N):
            ax = axes[r, c]
            if c >= len(tfs):
                ax.set_visible(False)
                continue
            tf = tfs[c]
            rt, at = rep[rep['TF'] == tf], age[age['TF'] == tf]
            for cond in CONDITIONS:
                draw(ax, rt[rt['condition'] == cond], at[at['condition'] == cond],
                     COND_COLORS[cond], cond)
            r_st = st.loc[tf]
            mark(ax, at, r_st)
            ax.set_title(f'{tf} ({at["n_targets_used"].iloc[0]})\n'
                         f'DR {r_st["DR_d16"]:+.2f} ({r_st["padj_DR_d16"]:.1g})\n'
                         f'slope {r_st["slope_per_day"]:+.3f} ({r_st["padj_slope"]:.1g})',
                         fontsize=6.5)
            ax.set_xticks(days)
            ax.set_xticklabels([f'P{d}' for d in days], fontsize=6, rotation=90)
            ax.set_xlim(days[0] - 1, days[-1] + 5)
            ax.set_ylim(bottom=0)
            ax.tick_params(labelsize=6)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            # x labels on the lowest visible panel of every column
            if r == len(ARCHS) - 1 or c >= len(rows[ARCHS[r + 1]]):
                ax.tick_params(labelbottom=True)
        axes[r, 0].set_ylabel(f'{a}\ntarget activity (/NR peak)', fontsize=7)
    axes[0, 0].legend(fontsize=6, frameon=False)
    fig.suptitle(f'Top {TOP_N} archetype-enriched Yoo25 L2/3 {SIGN} regulons per archetype '
                 '(rows A\', B\', C\'; by log2 enrichment): target activity in L2/3, NR vs DR\n'
                 'dots = replicates, line = replicate mean; title: OLS DR effect at day 16 and '
                 f'DR slope per day (padj, BH over TFs); bracket: padj < {FDR}, '
                 '* < 0.05, ** < 0.01, *** < 0.001', fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT_PDF, bbox_inches='tight')
    plt.close(fig)
    print(f'wrote {OUT_PDF}')


def plot_bar(rep, age, st):
    """Bars of rel. P21 DR change, largest decrease first, coloured by archetype; stars in one
    row at the top where the DR effect (OLS, P12-P21) has padj < FDR."""
    tf_change = st['rel_change_P21']
    nr_mean = age[(age['day'] == BAR_DAY) & (age['condition'] == 'NR')].set_index('TF')['mean']
    rep_change = rep[(rep['day'] == BAR_DAY) & (rep['condition'] == 'DR')]
    rep_change = rep_change.assign(rel_change=(rep_change['mean'] / rep_change['TF'].map(nr_mean) - 1).values)
    order = list(tf_change.sort_values().index)
    x = {tf: i for i, tf in enumerate(order)}
    y = tf_change[order] * 100
    fig, ax = plt.subplots(figsize=(0.24 * len(order) + 1.2, 3.0))
    ax.bar(range(len(order)), y, width=0.75, color=[ARCH_COLORS[st.loc[tf, 'arch']] for tf in order])
    ax.scatter(rep_change['TF'].map(x), rep_change['rel_change'] * 100, s=6, color='black',
               linewidths=0, zorder=3)
    lo = rep_change.groupby('TF')['rel_change'].min() * 100
    hi = rep_change.groupby('TF')['rel_change'].max() * 100
    y_star = max(hi.max(), y.max()) + 4      # one row of stars above everything
    for tf in order:
        ax.text(x[tf], y_star, stars(st.loc[tf, 'padj_DR_d16']), rotation=90, ha='center',
                va='bottom', fontsize=7)
    ax.axhline(0, color='black', lw=0.8)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90, fontsize=7)
    ax.set_ylim(min(lo.min(), y.min()) - 3, y_star + 9)   # room for the star row
    ax.set_xlim(-0.7, len(order) - 0.3)
    ax.set_ylabel(f'(P{BAR_DAY} DR - P{BAR_DAY} NR) / P{BAR_DAY} NR (%)\nregulon target activity',
                  fontsize=8)
    ax.tick_params(axis='y', labelsize=7)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend([plt.Rectangle((0, 0), 1, 1, color=ARCH_COLORS[a]) for a in ARCHS],
              [f'L2/3 {a}' for a in ARCHS], fontsize=7, frameon=False, loc='lower right')
    ax.set_title(f'Archetype-enriched Yoo25 L2/3 {SIGN} regulons: P{BAR_DAY} DR vs NR\n'
                 'bar = replicate-mean change; dots = DR replicates vs the NR mean; stars: DR effect '
                 f'(OLS, P12-P21) padj < {FDR}, BH over {len(st)} TFs; '
                 '* < 0.05, ** < 0.01, *** < 0.001', fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_PDF_BAR, bbox_inches='tight')
    plt.close(fig)
    print(f'wrote {OUT_PDF_BAR}')


rep = pd.read_csv(IN_REP, sep='\t')
age = pd.read_csv(IN_AGE, sep='\t')
st = pd.read_csv(IN_STATS, sep='\t').set_index('TF')
rows = {a: list(dict.fromkeys(age.loc[age['arch'] == a, 'TF'])) for a in ARCHS}   # by descending log2_enr
plot(rows, rep, age, st, sorted(age['day'].unique()))
plot_bar(rep, age, st)
