"""Per-cell L2/3 A'-C' axis for Yoo25 P21 NR (condition NR of the 0_prep/01 extract).

Per-cell archetype score = mean over the MARKER_CONDITION A' or C' archetype genes (= Table S4;
0_prep/02) of log2(CP10k+1) scaled per gene to [0, 1]: 0 = no expression, 1 = the gene's SCORE_MAX_PCTILE percentile across
cells (values above capped at 1). EXCLUDE_GENES (Cdh13 is an A' gene, Chrm2 a C' gene, Trpc6 is B')
are removed from the sets so that these marker genes, when plotted against the axis, do not help
define it.
A'-C' axis = score_C' - score_A' (negative = A' end, positive = C' end).

Reads:
  results/0_prep/01.yoo25_l23_p21nrdr.h5ad
  results/0_prep/02.yoo25_l23_abc_genes.tsv
Outputs:
  results/5_easifish/02.yoo25_p21nr_ac_axis.tsv  (cell, replicate, Type, score_A, score_C, ac_axis)
"""

import os
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_YOO25   = os.path.join(REPO_ROOT, 'results', '0_prep', '01.yoo25_l23_p21nrdr.h5ad')
INPUT_MARKERS = os.path.join(REPO_ROOT, 'results', '0_prep', '02.yoo25_l23_abc_genes.tsv')
OUT_TSV       = os.path.join(REPO_ROOT, 'results', '5_easifish', '02.yoo25_p21nr_ac_axis.tsv')

MARKER_CONDITION = 'P21'
EXCLUDE_GENES    = ['Cdh13', 'Trpc6', 'Chrm2']
SCORE_MAX_PCTILE = 98    # per-gene scale: 0 = no expression, this percentile = 1

os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)

a = ad.read_h5ad(INPUT_YOO25)
a = a[a.obs['condition'] == 'NR'].copy()
mk = pd.read_csv(INPUT_MARKERS, sep='\t', usecols=['cond', 'archetype', 'gene'])
mk = mk[mk['cond'] == MARKER_CONDITION]
sets = {k: sorted(set(g['gene']) - set(EXCLUDE_GENES)) for k, g in mk.groupby('archetype')}
for k in ['A', 'C']:
    miss = sorted(set(sets[k]) - set(a.var_names))
    if miss:
        raise ValueError(f"{k}' genes not in {INPUT_YOO25}: {miss}")

counts = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
lib = np.asarray(counts.sum(axis=1)).ravel()


def score(genes):
    idx = [a.var_names.get_loc(g) for g in genes]
    l2 = np.log2(counts[:, idx].toarray() / lib[:, None] * 1e4 + 1)
    hi = np.percentile(l2, SCORE_MAX_PCTILE, axis=0)
    hi = np.where(hi > 0, hi, 1.0)   # gene zero in >= 98% of cells: raw value, capped at 1
    return np.clip(l2 / hi, 0.0, 1.0).mean(axis=1)


out = pd.DataFrame({'replicate': a.obs['replicate'].astype(str).to_numpy(),
                    'Type': a.obs['Type'].astype(str).to_numpy(),
                    'score_A': score(sets['A']), 'score_C': score(sets['C'])},
                   index=pd.Index(a.obs_names, name='cell'))
out['ac_axis'] = out['score_C'] - out['score_A']
out.to_csv(OUT_TSV, sep='\t')
print(f"Yoo25 P21 NR: {a.n_obs} cells; A' {len(sets['A'])} and C' {len(sets['C'])} genes "
      f"(excluded {EXCLUDE_GENES})")
print(out.groupby('Type')['ac_axis'].describe().round(3).to_string())
print(f'wrote {OUT_TSV}')
