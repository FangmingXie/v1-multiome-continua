"""Archetype-enriched Yoo25 L2/3 regulons: target-gene activity over postnatal age, NR vs DR.

Regulons: every +/+ TF of 3_regulon/03's selection (developmental archetype genes, star rule;
A' 12, B' 13, C' 5), +/+ targets from 3_regulon/01, in archetype order and by descending log2_enr. The
expression loader and the per-cell score are l23.regulon's (as 3_regulon/06): each target's
log2(CP10k+1) divided by its NORM_MAX_PCTILE percentile over ALL Yoo25 L2/3 cells (every age, NR and DR
pooled, so both share one scale), capped at 1, averaged over the regulon's targets.

Replicates: per (condition, Age, Sample) the mean score over the sample's cells; per age, the mean over
replicates. Each TF is divided by its largest NR replicate-mean, so 1 = the NR peak (DR can exceed 1).

Statistics (l23.timecourse.dr_ols, as 2_development/04): per TF, OLS of the replicate scores (mean_norm)
over the P12-P21 samples, score ~ age (factor) + DR + DR*(day - D0), D0 = 16 = mean DR day: b_DR =
average DR effect at day 16 (in NR-peak units), b_slope = change in the DR effect per day; t-test
p-values, BH padj over all tested TFs. rel_change_P21 = (P21 DR - P21 NR) / P21 NR of the replicate-mean
score. Figures: 03.

Reads:
  results/0_prep/01.yoo25_l23.h5ad
  results/3_regulon/01.yoo25_l23_regulon_targets.tsv
  results/3_regulon/03.enriched_regulon_selection.tsv
Outputs:
  results/4_dr/01.regulon_target_expr_nr_vs_dr_replicates.tsv
  results/4_dr/01.regulon_target_expr_nr_vs_dr.tsv
  results/4_dr/01.regulon_target_expr_nr_vs_dr_stats.tsv
"""

import os

import pandas as pd

from l23.constants import ARCH_COLORS, TEST_DAYS
from l23.regulon import load_expression, regulon_score
from l23.timecourse import dr_ols, bh

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES_DIR = os.path.join(REPO_ROOT, 'results', '4_dr')
INPUT_YOO25 = os.path.join(REPO_ROOT, 'results', '0_prep', '01.yoo25_l23.h5ad')
INPUT_REGULONS = os.path.join(REPO_ROOT, 'results', '3_regulon', '01.yoo25_l23_regulon_targets.tsv')
INPUT_SELECTION = os.path.join(REPO_ROOT, 'results', '3_regulon', '03.enriched_regulon_selection.tsv')
OUT_REP_TSV = os.path.join(RES_DIR, '01.regulon_target_expr_nr_vs_dr_replicates.tsv')
OUT_TSV = os.path.join(RES_DIR, '01.regulon_target_expr_nr_vs_dr.tsv')
OUT_STATS = os.path.join(RES_DIR, '01.regulon_target_expr_nr_vs_dr_stats.tsv')

ARCHS = list(ARCH_COLORS)
SIGN = '+/+'
MIN_REPLICATES = 2
BAR_DAY = 21

os.makedirs(RES_DIR, exist_ok=True)


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
    """dr_ols per TF on the P12-P21 replicate scores; BH over TFs."""
    st = pd.DataFrame([{'TF': tf, 'arch': t['arch'].iloc[0],
                        'n_NR': int((t['condition'] == 'NR').sum()), 'n_DR': int((t['condition'] == 'DR').sum()),
                        **dr_ols(t, 'mean_norm')}
                       for tf, t in rep[rep['day'].isin(TEST_DAYS)].groupby('TF', sort=False)])
    for key in ['DR_d16', 'slope']:
        st[f'padj_{key}'] = bh(st[f'pvalue_{key}'])
    return st


def main():
    rows, targets = load_targets()
    L, cells = load_expression(INPUT_YOO25, sorted({g for t in targets.values() for g in t}))
    print(f'replicates per Age: {cells.groupby("Age")["Sample"].nunique().to_dict()}')

    out = [tf_stats(tf, a, targets[tf], L, cells) for a in ARCHS for tf in rows[a]]
    rep = pd.concat([o[0] for o in out], ignore_index=True)
    age = pd.concat([o[1] for o in out], ignore_index=True)
    rep.to_csv(OUT_REP_TSV, sep='\t', index=False)
    age.to_csv(OUT_TSV, sep='\t', index=False)
    print(f'wrote {OUT_REP_TSV}\nwrote {OUT_TSV}')

    st = dr_test(rep)
    a = age[age['day'] == BAR_DAY].pivot(index='TF', columns='condition', values='mean')
    st['rel_change_P21'] = st['TF'].map((a['DR'] - a['NR']) / a['NR'])
    st.to_csv(OUT_STATS, sep='\t', index=False)
    print(st.to_string(index=False, float_format=lambda v: f'{v:.3g}'))
    print(f'wrote {OUT_STATS}')


if __name__ == '__main__':
    main()
