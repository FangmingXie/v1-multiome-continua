"""Per-cell L2/3 A'-C' axis for Yoo25 P21 NR (Age P21 of the 0_prep/01 extract).

Per-cell archetype score = mean over the MARKER_CONDITION A' or C' archetype genes (= Table S4;
2_development/01) of log2(CP10k+1) scaled per gene to [0, 1]: 0 = no expression, 1 = the gene's SCORE_MAX_PCTILE percentile across
cells (values above capped at 1). EXCLUDE_GENES (Cdh13 is an A' gene, Chrm2 a C' gene, Trpc6 is B')
are removed from the sets so that these marker genes, when plotted against the axis, do not help
define it.
A'-C' axis = score_C' - score_A' (negative = A' end, positive = C' end).

Reads:
  results/0_prep/01.yoo25_l23.h5ad
  results/2_development/01.yoo25_l23_abc_genes.tsv
Outputs:
  results/5_easifish/01.yoo25_p21nr_ac_axis.tsv  (cell, replicate, Type, score_A, score_C, ac_axis)
"""

import os
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

from l23.expr import log2_cp10k, pctile_score

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_YOO25   = os.path.join(REPO_ROOT, 'results', '0_prep', '01.yoo25_l23.h5ad')
INPUT_MARKERS = os.path.join(REPO_ROOT, 'results', '2_development', '01.yoo25_l23_abc_genes.tsv')
OUT_TSV       = os.path.join(REPO_ROOT, 'results', '5_easifish', '01.yoo25_p21nr_ac_axis.tsv')

AGE              = 'P21'   # Yoo25 P21 NR
MARKER_CONDITION = 'P21'
EXCLUDE_GENES    = ['Cdh13', 'Trpc6', 'Chrm2']
SCORE_MAX_PCTILE = 98    # per-gene scale: 0 = no expression, this percentile = 1

os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)

a = ad.read_h5ad(INPUT_YOO25)
a = a[a.obs['Age'] == AGE].copy()
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
    return pctile_score(log2_cp10k(counts[:, idx].toarray(), lib), SCORE_MAX_PCTILE)


out = pd.DataFrame({'replicate': a.obs['Sample'].astype(str).to_numpy(),
                    'Type': a.obs['Type'].astype(str).to_numpy(),
                    'score_A': score(sets['A']), 'score_C': score(sets['C'])},
                   index=pd.Index(a.obs_names, name='cell'))
out['ac_axis'] = out['score_C'] - out['score_A']
out.to_csv(OUT_TSV, sep='\t')
print(f"Yoo25 P21 NR: {a.n_obs} cells; A' {len(sets['A'])} and C' {len(sets['C'])} genes "
      f"(excluded {EXCLUDE_GENES})")
print(out.groupby('Type')['ac_axis'].describe().round(3).to_string())
print(f'wrote {OUT_TSV}')
