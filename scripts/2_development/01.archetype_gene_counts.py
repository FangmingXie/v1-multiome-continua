"""Number of Yoo25 L2/3 A/B/C genes per condition (from 0_prep/02), NR vs DR; checked against Table S4.

Reads:
  results/0_prep/02.yoo25_l23_abc_genes.tsv
  data/TableS4_L23_archetypal_genes.tsv              (published Table S4; exact-match check only)
Outputs:
  figures/2_development/01.yoo25_l23_abc_gene_counts.pdf
"""

import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams['pdf.fonttype'] = 42

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_GENES = os.path.join(REPO_ROOT, 'results', '0_prep', '02.yoo25_l23_abc_genes.tsv')
INPUT_S4    = os.path.join(REPO_ROOT, 'data', 'TableS4_L23_archetypal_genes.tsv')
OUT_PDF     = os.path.join(REPO_ROOT, 'figures', '2_development', '01.yoo25_l23_abc_gene_counts.pdf')

NR = {'P6': 6, 'P8': 8, 'P10': 10, 'P12': 12, 'P14': 14, 'P17': 17, 'P21': 21}
DR = {'P12DR': 12, 'P14DR': 14, 'P17DR': 17, 'P21DR': 21}
ARCH_COLORS = {'A': 'C0', 'B': 'C1', 'C': 'C2'}
PANELS = ['AC', 'B']

res = pd.read_csv(INPUT_GENES, sep='\t')
s4 = pd.read_csv(INPUT_S4, sep='\t', usecols=['Condition', 'Archetype', 'Gene']).apply(lambda c: c.str.strip())
ours = res.groupby(['cond', 'archetype'])['gene'].apply(set)
theirs = s4.groupby(['Condition', 'Archetype'])['Gene'].apply(set)
diff = [k for k in theirs.index.union(ours.index) if ours.get(k, set()) != theirs.get(k, set())]
if diff:
    raise ValueError(f'{INPUT_GENES} differs from Table S4 in {len(diff)} cond x archetype: {diff}')
print('matches Table S4 exactly')

n = res.groupby(['archetype', 'cond']).size()
fig, axs = plt.subplots(1, len(PANELS), figsize=(4 * len(PANELS), 4))
for ax, ks in zip(axs, PANELS):
    for k in ks:
        ax.plot(list(NR.values()), [n.get((k, c), 0) for c in NR], '-o', fillstyle='none', c=ARCH_COLORS[k], label=f'{k} NR')
        ax.plot(list(DR.values()), [n.get((k, c), 0) for c in DR], '--s', fillstyle='none', c=ARCH_COLORS[k], label=f'{k} DR')
    ax.set(xticks=list(NR.values()), xlabel='time (P)', ylabel='num. of genes' if ks == PANELS[0] else None)
    ax.set_ylim(bottom=0)
    ax.legend()
    ax.spines[['top', 'right']].set_visible(False)
fig.tight_layout()

os.makedirs(os.path.dirname(OUT_PDF), exist_ok=True)
fig.savefig(OUT_PDF, bbox_inches='tight')
print(f'wrote {OUT_PDF}')
