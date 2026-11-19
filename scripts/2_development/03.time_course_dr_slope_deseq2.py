"""Does the dark-rearing (DR) effect change steadily with postnatal day? Yoo25 L2/3, P12-P21, per cell
group of 2_development/02 (all L2/3, A', B', C' cells), pseudobulk DESeq2.

Pseudobulk counts from 2_development/02 (raw counts summed per sample and cell group), restricted to the Ages with
both NR and DR (P12, P14, P17, P21; 9 NR and 8 DR samples). Sex-linked genes (EXCLUDE_GENES) are
dropped; genes need >= MIN_PB_COUNTS counts in >= MIN_PB_SAMPLES samples. Per gene,
negative-binomial GLM:
  log2 q = b0 + b_P14 [P14] + b_P17 [P17] + b_P21 [P21]      age as a factor (NR P12 = reference)
         + b_DR * DR                                         DR effect at day D0
         + b_int * DR * (day - D0)                           change in the DR effect per day
D0 = 16 = the mean day of the DR samples (2 each at 12, 14, 17, 21), so DR * (day - D0) sums to 0 over
the DR samples: b_DR is the average DR effect over the DR samples and b_int and b_DR are
nearly uncorrelated. The DR x (day - D0) column is built explicitly (dr_dayc): the formula term
condition:dayc expands to one slope per condition, collinear with the age factor (singular design).
Wald tests of b_int (slope) and b_DR, each with BH padj over all genes tested in the cell group
(log2FC = unshrunken MLE estimate).
Implied DR log2FC at P12 and P21 = b_DR + b_int * (day - D0) (reported, not tested).

Reads:
  results/2_development/02.yoo25_l23_pseudobulk_counts.tsv
  results/2_development/02.yoo25_l23_pseudobulk_samples.tsv
Outputs:
  results/2_development/03.yoo25_l23_time_course_dr_by_day_slope_deseq2.tsv
"""

import os
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from l23.constants import TEST_DAYS, D0, FDR

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_PB         = os.path.join(REPO_ROOT, 'results', '2_development', '02.yoo25_l23_pseudobulk_counts.tsv')
INPUT_PB_SAMPLES = os.path.join(REPO_ROOT, 'results', '2_development', '02.yoo25_l23_pseudobulk_samples.tsv')
OUT_TSV          = os.path.join(REPO_ROOT, 'results', '2_development', '03.yoo25_l23_time_course_dr_by_day_slope_deseq2.tsv')

GENES          = ['Cdh13', 'Robo1', 'Pcdh9', 'Baz1a']
GROUPS         = ['all L2/3', "A' cells", "B' cells", "C' cells"]
DAYS           = TEST_DAYS   # Ages with both NR and DR; D0 = their mean DR day (checked below)
MIN_PB_COUNTS  = 10
MIN_PB_SAMPLES = 2
# Sex-linked genes.
EXCLUDE_GENES = ['Ddx3y', 'Eif2s3y', 'Kdm5d', 'Uty', 'Usp9y', 'Zfy1', 'Zfy2', 'Sry',
                 'Ube1y1', 'Rbmy', 'Tspy1', 'Xist', 'Tsix']

os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)

pb = pd.read_csv(INPUT_PB, sep='\t', index_col=['cells', 'Sample'])
meta_all = pd.read_csv(INPUT_PB_SAMPLES, sep='\t', index_col=['cells', 'Sample'])
if not pb.index.equals(meta_all.index):
    raise ValueError(f'rows of {INPUT_PB} and {INPUT_PB_SAMPLES} differ')
if sorted(pb.index.get_level_values('cells').unique()) != sorted(GROUPS):
    raise ValueError(f'expected cell groups {GROUPS} in {INPUT_PB}')
keep = meta_all['day'].isin(DAYS).to_numpy()
if sorted(meta_all[keep].groupby('day')['condition'].nunique().unique()) != [2]:
    raise ValueError(f'expected NR and DR at every day in {DAYS}')
pb, meta_all = pb[keep], meta_all[keep]
pb = pb.loc[:, ~pb.columns.isin(EXCLUDE_GENES)]

out = []
for grp in GROUPS:
    counts = pb.loc[grp]
    counts = counts.loc[:, (counts >= MIN_PB_COUNTS).sum(axis=0) >= MIN_PB_SAMPLES]
    missing = [g for g in GENES if g not in counts.columns]
    if missing:
        raise ValueError(f'{grp}: {missing} fail the pseudobulk count filter')
    smp = meta_all.loc[grp].loc[counts.index]
    if sorted(smp['day'].unique()) != DAYS:
        raise ValueError(f'{grp}: expected days {DAYS}, got {sorted(smp["day"].unique())}')
    dr = (smp['condition'] == 'DR').astype(float)
    if smp.loc[dr == 1, 'day'].mean() != D0:
        raise ValueError(f'{grp}: mean DR day is {smp.loc[dr == 1, "day"].mean()}, D0 = {D0}')
    meta = pd.DataFrame({'age': pd.Categorical([f'P{d}' for d in smp['day']], categories=[f'P{d}' for d in DAYS]),
                         'condition': pd.Categorical(smp['condition'], categories=['NR', 'DR']),
                         'dr_dayc': dr * (smp['day'] - D0)}, index=counts.index)
    print(f'=== {grp}: {len(counts)} samples ({meta["condition"].value_counts().to_dict()}), '
          f'{counts.shape[1]} genes tested', flush=True)
    dds = DeseqDataSet(counts=counts, metadata=meta, design='~age + condition + dr_dayc', quiet=True)
    dds.deseq2()
    cols = list(dds.obsm['design_matrix'].columns)
    expected = ['Intercept', 'age[T.P14]', 'age[T.P17]', 'age[T.P21]', 'condition[T.DR]', 'dr_dayc']
    if cols != expected:
        raise ValueError(f'unexpected design columns {cols}; expected {expected}')

    res = {}
    for key, coef in [('DR_d16', 'condition[T.DR]'), ('slope', 'dr_dayc')]:
        st = DeseqStats(dds, contrast=np.array([float(c == coef) for c in cols]), alpha=FDR, quiet=True)
        st.summary()
        res[key] = st.results_df
    r = pd.DataFrame({'baseMean': res['DR_d16']['baseMean'],
                      'log2FC_DR_d16': res['DR_d16']['log2FoldChange'], 'lfcSE_DR_d16': res['DR_d16']['lfcSE'],
                      'pvalue_DR_d16': res['DR_d16']['pvalue'], 'padj_DR_d16': res['DR_d16']['padj'],
                      'slope_per_day': res['slope']['log2FoldChange'], 'lfcSE_slope': res['slope']['lfcSE'],
                      'pvalue_slope': res['slope']['pvalue'], 'padj_slope': res['slope']['padj']})
    for d in (DAYS[0], DAYS[-1]):
        r[f'log2FC_DR_P{d}_implied'] = r['log2FC_DR_d16'] + r['slope_per_day'] * (d - D0)
    r.insert(0, 'n_genes_tested', counts.shape[1])
    r.insert(0, 'n_DR', int((meta['condition'] == 'DR').sum()))
    r.insert(0, 'n_NR', int((meta['condition'] == 'NR').sum()))
    r.insert(0, 'cells', grp)
    out.append(r.rename_axis('gene').reset_index())

out = pd.concat(out)
out.to_csv(OUT_TSV, sep='\t', index=False)
show = ['gene', 'cells', 'log2FC_DR_d16', 'padj_DR_d16', 'slope_per_day', 'lfcSE_slope', 'pvalue_slope', 'padj_slope',
        'log2FC_DR_P12_implied', 'log2FC_DR_P21_implied']
print(out.loc[out['gene'].isin(GENES), show].to_string(index=False, float_format=lambda v: f'{v:.3g}'))
print(f'genes with padj < {FDR} per cell group:',
      {k: {'DR_d16': int((t['padj_DR_d16'] < FDR).sum()), 'slope': int((t['padj_slope'] < FDR).sum())}
       for k, t in out.groupby('cells')})
print(f'wrote {OUT_TSV}')
