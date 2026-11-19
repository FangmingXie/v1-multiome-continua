"""GENES (rows) x cell groups (all L2/3, A', B', C' cells; columns) of 2_development/02 over postnatal day, NR vs DR,
with the 2_development/03 statistics. A'/B'/C' cells = top 20% per Age by A/B/C scores from all
Table S4 conditions' genes merged, cells in more than one group dropped (see 2_development/02).
y = per-sample pseudobulk log2(CPM+1) = log2(counts / lib * 1e6 + 1) (2_development/02: counts and library sizes
summed over the cells) minus the gene's reference REF (mean over the P21 NR samples, all L2/3 cells;
dashed line at 0), labeled 'log2(CPM+1) - P21 NR mean' (= a log2FC vs P21 NR L2/3 where CPM >> 1);
the scale DESeq2 tests in 2_development/03, so NR-DR gaps match the log2FCs in the titles; archetype panels show
enrichment relative to all L2/3. Row label = gene; one shared y-axis label on the left.
Lines pass through the mean across samples per age x condition; dots are samples. y is shared within
a row.

Statistics (2_development/03: pseudobulk DESeq2 over P12-P21, 9 NR vs 8 DR samples, ~ age + condition + DR*(day - 16)):
  DR     b_DR, the average DR vs NR log2FC over the DR samples (the DR effect at day 16)
  slope  b_int, change in the DR log2FC per day
Each panel title gives both (estimate, padj; padj NA = removed by DESeq2 independent filtering).
Significant: DR = padj < FDR and |log2FC| > MIN_ABS_LOG2FC; slope = padj < FDR. For a significant
test, a bracket right of P21 spans the NR and DR line ends, labeled 'DR' and / or 'slope' with
* < 0.05, ** < 0.01, *** < 0.001 (padj).

Reads:
  results/2_development/02.yoo25_l23_pseudobulk_counts.tsv
  results/2_development/02.yoo25_l23_pseudobulk_samples.tsv
  results/2_development/03.yoo25_l23_time_course_dr_by_day_slope_deseq2.tsv
Outputs:
  figures/4_dr/02.yoo25_l23_genes_time_course_pseudobulk_log2fc.pdf
"""

import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

from l23.constants import FDR
from l23.timecourse import stars, draw_nr_dr

matplotlib.rcParams['pdf.fonttype'] = 42

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_PB         = os.path.join(REPO_ROOT, 'results', '2_development', '02.yoo25_l23_pseudobulk_counts.tsv')
INPUT_PB_SAMPLES = os.path.join(REPO_ROOT, 'results', '2_development', '02.yoo25_l23_pseudobulk_samples.tsv')
INPUT_STATS      = os.path.join(REPO_ROOT, 'results', '2_development', '03.yoo25_l23_time_course_dr_by_day_slope_deseq2.tsv')
OUT_FIG          = os.path.join(REPO_ROOT, 'figures', '4_dr', '02.yoo25_l23_genes_time_course_pseudobulk_log2fc.pdf')

GENES       = ['Cdh13', 'Robo1', 'Pcdh9', 'Baz1a']
GROUPS      = ['all L2/3', "A' cells", "B' cells", "C' cells"]
REF         = ('all L2/3', 'NR', 21)   # (cells, condition, day) of the pseudobulk log2FC reference
MIN_ABS_LOG2FC = 0.6   # DR effect also needs |log2FC| above this (~1.5-fold) to be marked

os.makedirs(os.path.dirname(OUT_FIG), exist_ok=True)

absent = sorted(set(GENES) - set(pd.read_csv(INPUT_PB, sep='\t', nrows=0).columns))
if absent:
    raise ValueError(f'genes {absent} not in {INPUT_PB}')
pb = pd.read_csv(INPUT_PB, sep='\t', index_col=['cells', 'Sample'], usecols=['cells', 'Sample', *GENES])
meta = pd.read_csv(INPUT_PB_SAMPLES, sep='\t', index_col=['cells', 'Sample'])
if not pb.index.equals(meta.index):
    raise ValueError(f'rows of {INPUT_PB} and {INPUT_PB_SAMPLES} differ')
if sorted(meta.index.get_level_values('cells').unique()) != sorted(GROUPS):
    raise ValueError(f'expected cell groups {GROUPS} in {INPUT_PB_SAMPLES}')
samp = (np.log2(pb.div(meta['lib'], axis=0) * 1e6 + 1).join(meta[['condition', 'day']])
        .melt(id_vars=['condition', 'day'], value_vars=GENES, var_name='gene', value_name='pb_log2cpm',
              ignore_index=False).reset_index())
ref = samp[(samp['cells'] == REF[0]) & (samp['condition'] == REF[1]) & (samp['day'] == REF[2])]
ref_val = ref.groupby('gene')['pb_log2cpm'].mean()
if sorted(ref_val.index) != sorted(GENES):
    raise ValueError(f'no reference samples {REF} for some of {GENES} in {INPUT_PB}')
samp['pb_log2fc'] = samp['pb_log2cpm'] - samp['gene'].map(ref_val)
print(f'reference {REF}, pseudobulk log2(CPM+1) per gene:', ref_val.round(3).to_dict())

st = pd.read_csv(INPUT_STATS, sep='\t').set_index(['gene', 'cells'])
absent = [(g, c) for g in GENES for c in GROUPS if (g, c) not in st.index]
if absent:
    raise ValueError(f'{absent} not in {INPUT_STATS}')


def marks(k):
    """Labels of the significant tests of one panel: DR and slope (2_development/03)."""
    out = []
    if abs(k['log2FC_DR_d16']) > MIN_ABS_LOG2FC and stars(k['padj_DR_d16']):
        out.append(f"DR {stars(k['padj_DR_d16'])}")
    if stars(k['padj_slope']):
        out.append(f"slope {stars(k['padj_slope'])}")
    return out


def fmt_p(padj):
    return 'NA' if np.isnan(padj) else f'{padj:.2g}'


def draw(ycol, ylabel, out):
    fig, axes = plt.subplots(len(GENES), len(GROUPS), figsize=(2.5 * len(GROUPS), 3.2 * len(GENES)),
                             sharey='row', squeeze=False)
    for i, (row, gene) in enumerate(zip(axes, GENES)):
        for ax, grp in zip(row, GROUPS):
            k = st.loc[(gene, grp)]
            draw_nr_dr(ax, samp[(samp['gene'] == gene) & (samp['cells'] == grp)], ycol, sorted(samp['day'].unique()),
                       marks(k))
            ax.set_title(f'DR {k["log2FC_DR_d16"]:.2f} ({fmt_p(k["padj_DR_d16"])})\n'
                         f'slope {k["slope_per_day"]:+.3f} ({fmt_p(k["padj_slope"])})', fontsize=6.5)
            if i == 0:
                ax.text(0.5, 1.22, grp, transform=ax.transAxes, ha='center', va='bottom', fontsize=10)   # column header
        row[0].set_ylabel(gene, fontsize=11, fontstyle='italic')
    for ax in axes[-1]:
        ax.set_xlabel('postnatal day')
    axes[0, 0].legend(fontsize=8, frameon=False, loc='lower left')
    fig.supylabel(ylabel)
    fig.suptitle('NR vs DR, pseudobulk; reference = P21 NR mean, all L2/3 (dashed line)\n'
                 'DESeq2 over P12-P21, ~ age + condition + DR·(day - 16): '
                 'DR = average DR log2FC (at day 16); slope = change in DR log2FC per day. Panel: estimate (padj)\n'
                 f'marked (bracket): DR padj < {FDR} and |log$_2$FC| > {MIN_ABS_LOG2FC:g}; slope padj < {FDR}. '
                 '* padj < 0.05, ** < 0.01, *** < 0.001', fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out)
    print(f'wrote {out}')


draw('pb_log2fc', 'log$_2$(CPM+1) − P21 NR mean', OUT_FIG)
