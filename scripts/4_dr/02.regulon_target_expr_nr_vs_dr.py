"""Top archetype-enriched Yoo25 L2/3 regulons: target-gene activity over postnatal age, NR vs DR.

Regulons: every +/+ TF of 3_regulon/04's selection (developmental archetype genes, star rule;
A' 12, B' 13, C' 5), +/+ targets from 3_regulon/01. The line figure shows, per archetype row, the
top TOP_N by log2_enr (C' has only 5, so its row is shorter); statistics and the barplot cover all
of them. The expression loader and the per-cell score are l23.regulon's (as 1_mouse_human/13-14):
each target's log2(CP10k+1) divided by its NORM_MAX_PCTILE percentile
over ALL Yoo25 L2/3 cells (every age, NR and DR pooled, so both share one scale), capped at 1,
averaged over the regulon's targets.

Replicates: per (condition, Age, Sample) the mean score over the sample's cells (one dot per
replicate); the line is the mean over replicates. Each TF is divided by its largest NR
replicate-mean, so 1 = the NR peak (DR can exceed 1).

Statistics (as 2_development/04): per TF, OLS of the replicate scores (mean_norm) over the P12-P21 samples,
score ~ age (factor) + DR + DR*(day - D0), D0 = 16 = mean DR day: b_DR = average DR effect at
day 16 (in NR-peak units), b_slope = change in the DR effect per day; t-test p-values, BH padj
over all tested TFs. Marked (bracket right of P21, as 2_development/05): padj < FDR.

Figures
  lines:   one panel per TF, NR (C1) vs DR (black), own y-axis; one row per archetype.
  barplot: every TF, height = (P21 DR - P21 NR) / P21 NR of the replicate-mean score
           (rel_change_P21); dots = each DR replicate against the P21 NR mean. Bars coloured by
           archetype, ordered by value (largest decrease first); stars = DR effect padj (OLS above).

Reads:
  data/raw/superdupermegaRNA_hasraw.h5ad
  results/3_regulon/01.yoo25_l23_regulon_targets.tsv
  results/3_regulon/04.enriched_regulon_selection.tsv
Outputs:
  results/4_dr/02.regulon_target_expr_nr_vs_dr_replicates.tsv
  results/4_dr/02.regulon_target_expr_nr_vs_dr.tsv
  results/4_dr/02.regulon_target_expr_nr_vs_dr_stats.tsv
  figures/4_dr/02.regulon_target_expr_nr_vs_dr.pdf
  figures/4_dr/02.regulon_target_expr_p21_dr_change_barplot.pdf
"""

import os

import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt

from l23.regulon import load_expression, regulon_score

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES_DIR = os.path.join(REPO_ROOT, 'results', '4_dr')
FIG_DIR = os.path.join(REPO_ROOT, 'figures', '4_dr')
INPUT_SOURCE = os.path.join(REPO_ROOT, 'data', 'raw', 'superdupermegaRNA_hasraw.h5ad')
INPUT_REGULONS = os.path.join(REPO_ROOT, 'results', '3_regulon', '01.yoo25_l23_regulon_targets.tsv')
INPUT_SELECTION = os.path.join(REPO_ROOT, 'results', '3_regulon', '04.enriched_regulon_selection.tsv')
OUT_REP_TSV = os.path.join(RES_DIR, '02.regulon_target_expr_nr_vs_dr_replicates.tsv')
OUT_TSV = os.path.join(RES_DIR, '02.regulon_target_expr_nr_vs_dr.tsv')
OUT_STATS = os.path.join(RES_DIR, '02.regulon_target_expr_nr_vs_dr_stats.tsv')
OUT_PDF = os.path.join(FIG_DIR, '02.regulon_target_expr_nr_vs_dr.pdf')
OUT_PDF_BAR = os.path.join(FIG_DIR, '02.regulon_target_expr_p21_dr_change_barplot.pdf')

ARCHS = ["A'", "B'", "C'"]
TOP_N = 8                      # TFs per archetype row = panels per row
SIGN = '+/+'
CONDITIONS = ['NR', 'DR']
COND_COLORS = {'NR': 'C1', 'DR': 'black'}
REP_DOT_SIZE = 9
MIN_REPLICATES = 2
# NR vs DR test, as 2_development/04
TEST_DAYS = [12, 14, 17, 21]   # Ages with both NR and DR
D0 = 16                        # centering day: mean day of the DR samples (checked)
FDR = 0.05
STARS = [(0.001, '***'), (0.01, '**'), (0.05, '*')]
BAR_DAY = 21
ARCH_COLORS = {"A'": 'C0', "B'": 'C1', "C'": 'C2'}

os.makedirs(RES_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)


def load_targets():
    """{arch: [TF, ...]} (all, descending log2_enr) and {TF: sorted +/+ targets}."""
    sel = pd.read_csv(INPUT_SELECTION, sep='\t')
    if set(sel['arch_letter']) != set(ARCHS):
        raise ValueError(f'unexpected archetypes in {INPUT_SELECTION}: {sorted(set(sel["arch_letter"]))}')
    rows = {a: list(sel[sel['arch_letter'] == a].sort_values('log2_enr', ascending=False)['TF'])
            for a in ARCHS}
    reg = pd.read_csv(INPUT_REGULONS, sep='\t')
    reg = reg[reg['regulation_direction'] == SIGN]
    targets = {tf: sorted(set(reg.loc[reg['TF'] == tf, 'Gene'])) for tfs in rows.values() for tf in tfs}
    empty = [tf for tf, t in targets.items() if not t]
    if empty:
        raise ValueError(f'no {SIGN} regulon in {INPUT_REGULONS} for {empty}')
    for a, tfs in rows.items():
        print(f'  {a}: {", ".join(tfs)}')
    return rows, targets


def tf_stats(tf, arch, genes, L, cells):
    """(per-replicate, per-age) tables of the regulon score, both scaled to the NR peak."""
    score, n_used = regulon_score(tf, genes, L)
    rep = (cells.assign(score=score).groupby(['condition', 'Age', 'day', 'Sample'])['score']
                .agg(n_cells='size', mean='mean').reset_index())
    age = (rep.groupby(['condition', 'Age', 'day'])['mean']
              .agg(n_replicates='size', mean='mean', sd='std').reset_index())
    few = age[age['n_replicates'] < MIN_REPLICATES]
    if len(few):
        raise ValueError(f'{tf}: fewer than {MIN_REPLICATES} replicates at {few["Age"].tolist()}')
    peak = age.loc[age['condition'] == 'NR', 'mean'].max()
    if not peak > 0:
        raise ValueError(f'{tf}: NR peak mean is {peak}; cannot max-scale')
    rep['mean_norm'] = rep['mean'] / peak
    age['mean_norm'] = age['mean'] / peak
    age['sd_norm'] = age['sd'] / peak
    for df in (rep, age):
        df.insert(0, 'TF', tf)
        df.insert(1, 'arch', arch)
        df.insert(2, 'n_targets_used', n_used)
    return rep, age.sort_values(['condition', 'day'])


def dr_test(rep):
    """2_development/04's OLS per TF on the P12-P21 replicate scores; BH over TFs."""
    out = []
    for tf, t in rep[rep['day'].isin(TEST_DAYS)].groupby('TF', sort=False):
        t = t.assign(age=pd.Categorical(t['day'].map(lambda d: f'P{d}'),
                                        categories=[f'P{d}' for d in TEST_DAYS]),
                     dr=(t['condition'] == 'DR').astype(float))
        if t.loc[t['dr'] == 1, 'day'].mean() != D0:
            raise ValueError(f'{tf}: mean DR day is {t.loc[t["dr"] == 1, "day"].mean()}, D0 = {D0}')
        t['dr_dayc'] = t['dr'] * (t['day'] - D0)
        fit = smf.ols('mean_norm ~ age + dr + dr_dayc', data=t).fit()
        out.append({'TF': tf, 'arch': t['arch'].iloc[0],
                    'n_NR': int((t['dr'] == 0).sum()), 'n_DR': int((t['dr'] == 1).sum()),
                    'DR_d16': fit.params['dr'], 'se_DR_d16': fit.bse['dr'],
                    'pvalue_DR_d16': fit.pvalues['dr'],
                    'slope_per_day': fit.params['dr_dayc'], 'se_slope': fit.bse['dr_dayc'],
                    'pvalue_slope': fit.pvalues['dr_dayc']})
    st = pd.DataFrame(out)
    for key in ['DR_d16', 'slope']:
        st[f'padj_{key}'] = multipletests(st[f'pvalue_{key}'], method='fdr_bh')[1]
    return st


def stars(padj):
    return next((s for cut, s in STARS if padj < cut), '')


def mark(ax, age, r):
    """Bracket right of the last day (as 2_development/05), labelled with the significant terms."""
    lab = [f'{name} {stars(r[f"padj_{key}"])}' for name, key in [('DR', 'DR_d16'), ('slope', 'slope')]
           if stars(r[f'padj_{key}'])]
    if not lab:
        return
    last = max(TEST_DAYS)
    ends = {c: age.loc[(age['condition'] == c) & (age['day'] == last), 'mean_norm'].iloc[0]
            for c in CONDITIONS}
    ax.plot([last + 0.8, last + 1.1, last + 1.1, last + 0.8],
            [ends['NR'], ends['NR'], ends['DR'], ends['DR']], color='black', lw=0.8)
    ax.text(last + 1.4, (ends['NR'] + ends['DR']) / 2, '\n'.join(lab), fontsize=6,
            va='center', ha='left')


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
            r_st = st.set_index('TF').loc[tf]
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


def p21_change(rep, age):
    """(DR - NR) / NR at BAR_DAY of the replicate-mean score, per TF and per DR replicate."""
    a = age[age['day'] == BAR_DAY].pivot(index='TF', columns='condition', values='mean')
    tf_change = (a['DR'] - a['NR']) / a['NR']
    r = rep[(rep['day'] == BAR_DAY) & (rep['condition'] == 'DR')]
    rep_change = r['mean'] / r['TF'].map(a['NR']) - 1
    return tf_change, r.assign(rel_change=rep_change.values)[['TF', 'Sample', 'rel_change']]


def plot_bar(tf_change, rep_change, st):
    """Bars of rel. P21 DR change, largest decrease first, coloured by archetype; stars in one
    row at the top where the DR effect (OLS, P12-P21) has padj < FDR."""
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


def main():
    rows, targets = load_targets()
    L, cells = load_expression(INPUT_SOURCE, sorted({g for t in targets.values() for g in t}))
    print(f'replicates per Age: {cells.groupby("Age")["Sample"].nunique().to_dict()}')

    out = [tf_stats(tf, a, targets[tf], L, cells) for a in ARCHS for tf in rows[a]]
    rep = pd.concat([o[0] for o in out], ignore_index=True)
    age = pd.concat([o[1] for o in out], ignore_index=True)
    rep.to_csv(OUT_REP_TSV, sep='\t', index=False)
    age.to_csv(OUT_TSV, sep='\t', index=False)
    print(f'wrote {OUT_REP_TSV}\nwrote {OUT_TSV}')

    tf_change, rep_change = p21_change(rep, age)
    st = dr_test(rep)
    st['rel_change_P21'] = st['TF'].map(tf_change)
    st.to_csv(OUT_STATS, sep='\t', index=False)
    print(st.to_string(index=False, float_format=lambda v: f'{v:.3g}'))
    print(f'wrote {OUT_STATS}')

    plot(rows, rep, age, st, sorted(age['day'].unique()))
    plot_bar(tf_change, rep_change, st.set_index('TF'))


if __name__ == '__main__':
    main()
