"""A/B/C archetype gene modules (M1-M5) over postnatal day in Yoo25 L2/3, NR vs DR, per cell group of
2_development/01 (all L2/3, A', B', C' cells).

Modules: each archetype gene set (gene_sets_l23{a,b,c}_250423.csv; = the union of that archetype's
Table S4 genes (0_prep/02) over all conditions, checked) is split by its column clst 0-4, named M1-M5 here. A gene
in two archetype sets is in a module of each.
Per gene and sample (2_development/01 pseudobulk): log2(CPM+1) = log2(counts / lib * 1e6 + 1) minus the gene's
reference REF (mean over the P21 NR samples, all L2/3 cells), as in 4_dr/02. Module score per sample =
mean of that value over the module's genes that pass the filter in the cell group (>= MIN_PB_COUNTS
counts in >= MIN_PB_SAMPLES of the group's samples, all Ages).
Statistics per module and cell group (l23.timecourse.dr_ols): OLS of the module score over the P12-P21
samples (9 NR, 8 DR) with the 2_development/02 model, score ~ age (factor) + DR + DR*(day - D0), D0 = 16:
b_DR = average DR effect (at day 16), b_int = change in the DR effect per day; t-test p-values, BH padj
over the 15 modules of the cell group. Context: fraction of the module's genes with padj < FDR in
2_development/02 (DR, slope; columns *_14).
Figures: 2_development/06.

Reads:
  results/2_development/01.yoo25_l23_pseudobulk_counts.tsv
  results/2_development/01.yoo25_l23_pseudobulk_samples.tsv
  results/2_development/02.yoo25_l23_time_course_dr_by_day_slope_deseq2.tsv
  data/raw/gene_sets_l23{a,b,c}_250423.csv
  results/0_prep/02.yoo25_l23_abc_genes.tsv
Outputs:
  results/2_development/03.yoo25_l23_gene_modules_genes.tsv     (gene -> module, pass filter per cell group)
  results/2_development/03.yoo25_l23_gene_modules_samples.tsv   (module score per sample)
  results/2_development/03.yoo25_l23_gene_modules_stats.tsv     (OLS DR and slope per module)
"""

import os
import numpy as np
import pandas as pd

from l23.constants import TEST_DAYS, FDR
from l23.timecourse import dr_ols, bh

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_PB         = os.path.join(REPO_ROOT, 'results', '2_development', '01.yoo25_l23_pseudobulk_counts.tsv')
INPUT_PB_SAMPLES = os.path.join(REPO_ROOT, 'results', '2_development', '01.yoo25_l23_pseudobulk_samples.tsv')
INPUT_STATS      = os.path.join(REPO_ROOT, 'results', '2_development', '02.yoo25_l23_time_course_dr_by_day_slope_deseq2.tsv')
INPUT_SETS       = {k: os.path.join(REPO_ROOT, 'data', 'raw', f'gene_sets_l23{k.lower()}_250423.csv')
                    for k in ['A', 'B', 'C']}
INPUT_MARKERS    = os.path.join(REPO_ROOT, 'results', '0_prep', '02.yoo25_l23_abc_genes.tsv')
OUT_GENES        = os.path.join(REPO_ROOT, 'results', '2_development', '03.yoo25_l23_gene_modules_genes.tsv')
OUT_SAMPLES      = os.path.join(REPO_ROOT, 'results', '2_development', '03.yoo25_l23_gene_modules_samples.tsv')
OUT_STATS        = os.path.join(REPO_ROOT, 'results', '2_development', '03.yoo25_l23_gene_modules_stats.tsv')

GROUPS         = ['all L2/3', "A' cells", "B' cells", "C' cells"]
ARCHETYPES     = ['A', 'B', 'C']
MODULES        = ['M1', 'M2', 'M3', 'M4', 'M5']   # = clst 0-4
REF            = ('all L2/3', 'NR', 21)   # (cells, condition, day) of the per-gene reference
MIN_PB_COUNTS  = 10
MIN_PB_SAMPLES = 2

os.makedirs(os.path.dirname(OUT_STATS), exist_ok=True)

# modules
mod = pd.concat([pd.read_csv(p, index_col=0).assign(archetype=k) for k, p in INPUT_SETS.items()],
                ignore_index=True)
if sorted(mod['clst'].unique()) != list(range(len(MODULES))):
    raise ValueError(f'expected clst 0-{len(MODULES) - 1}, got {sorted(mod["clst"].unique())}')
if mod.duplicated(['archetype', 'gene']).any():
    raise ValueError('a gene is listed twice within one archetype gene set')
mod['module'] = mod['clst'].map(lambda c: MODULES[c])
mk = pd.read_csv(INPUT_MARKERS, sep='\t', usecols=['archetype', 'gene'])
for k in ARCHETYPES:
    if set(mod.loc[mod['archetype'] == k, 'gene']) != set(mk.loc[mk['archetype'] == k, 'gene']):
        raise ValueError(f'{INPUT_SETS[k]} genes differ from the {k} genes of {INPUT_MARKERS}')
print('genes per module:', mod.groupby(['archetype', 'module']).size().unstack().to_string(), sep='\n')

# per-gene log2(CPM+1) minus the P21 NR reference
genes = sorted(mod['gene'].unique())
absent = sorted(set(genes) - set(pd.read_csv(INPUT_PB, sep='\t', nrows=0).columns))
if absent:
    raise ValueError(f'{len(absent)} module genes not in {INPUT_PB}: {absent}')
pb = pd.read_csv(INPUT_PB, sep='\t', index_col=['cells', 'Sample'], usecols=['cells', 'Sample', *genes])
meta = pd.read_csv(INPUT_PB_SAMPLES, sep='\t', index_col=['cells', 'Sample'])
if not pb.index.equals(meta.index):
    raise ValueError(f'rows of {INPUT_PB} and {INPUT_PB_SAMPLES} differ')
if sorted(meta.index.get_level_values('cells').unique()) != sorted(GROUPS):
    raise ValueError(f'expected cell groups {GROUPS} in {INPUT_PB_SAMPLES}')
logcpm = np.log2(pb.div(meta['lib'], axis=0) * 1e6 + 1)
is_ref = ((meta.index.get_level_values('cells') == REF[0]) & (meta['condition'] == REF[1])
          & (meta['day'] == REF[2])).to_numpy()
if not is_ref.any():
    raise ValueError(f'no reference samples {REF} in {INPUT_PB_SAMPLES}')
rel = logcpm - logcpm[is_ref].mean(axis=0)

# gene filter per cell group, module scores per sample
passes = pd.DataFrame({g: (pb.loc[g] >= MIN_PB_COUNTS).sum(axis=0) >= MIN_PB_SAMPLES for g in GROUPS})
gene_tab = mod[['archetype', 'module', 'gene']].join(passes, on='gene')
gene_tab.to_csv(OUT_GENES, sep='\t', index=False)
print(f'wrote {OUT_GENES}')
rows = []
for grp in GROUPS:
    for (k, m), t in mod.groupby(['archetype', 'module']):
        keep = [g for g in t['gene'] if passes.loc[g, grp]]
        if not keep:
            raise ValueError(f'{grp}: no gene of {k} {m} passes the filter')
        rows.append(meta.loc[grp, ['Age', 'condition', 'day']]
                    .assign(score=rel.loc[grp, keep].mean(axis=1), n_genes=len(keep), n_genes_set=len(t),
                            cells=grp, archetype=k, module=m).reset_index())
samp = pd.concat(rows, ignore_index=True)[['cells', 'archetype', 'module', 'Sample', 'Age', 'condition', 'day',
                                           'n_genes', 'n_genes_set', 'score']]
samp.to_csv(OUT_SAMPLES, sep='\t', index=False)
print(f'wrote {OUT_SAMPLES}')

# OLS per module and cell group (2_development/02 model), P12-P21
deg = pd.read_csv(INPUT_STATS, sep='\t').set_index(['cells', 'gene'])
stats = []
for (grp, k, m), t in samp[samp['day'].isin(TEST_DAYS)].groupby(['cells', 'archetype', 'module']):
    keep = [g for g in mod.loc[(mod['archetype'] == k) & (mod['module'] == m), 'gene'] if passes.loc[g, grp]]
    d = deg.loc[grp].reindex(keep)
    stats.append({'cells': grp, 'archetype': k, 'module': m, 'n_genes': len(keep),
                  'n_NR': int((t['condition'] == 'NR').sum()), 'n_DR': int((t['condition'] == 'DR').sum()),
                  **dr_ols(t, 'score'),
                  'n_genes_in_14': int(d['padj_DR_d16'].notna().sum()),
                  'frac_genes_DR_14': (d['padj_DR_d16'] < FDR).sum() / d['padj_DR_d16'].notna().sum(),
                  'frac_genes_slope_14': (d['padj_slope'] < FDR).sum() / d['padj_slope'].notna().sum()})
stats = pd.DataFrame(stats)
for key in ['DR_d16', 'slope']:
    stats[f'padj_{key}'] = stats.groupby('cells')[f'pvalue_{key}'].transform(bh)
stats.to_csv(OUT_STATS, sep='\t', index=False)
print(stats[['cells', 'archetype', 'module', 'n_genes', 'DR_d16', 'padj_DR_d16', 'slope_per_day', 'padj_slope',
             'frac_genes_DR_14', 'frac_genes_slope_14']].to_string(index=False, float_format=lambda v: f'{v:.3g}'))
print(f'wrote {OUT_STATS}')
