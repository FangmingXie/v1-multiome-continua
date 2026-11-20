"""Mouse L2/3 archetype scores on human L2/3 IT (Jorstad23) cells.

Each mouse archetype's marker set (script 02) is mapped through 1-to-1 human/mouse orthologs and
scored on the human cells with the same recipe as the mouse scores: per gene min-max clipped at the
SCORE_PCTILE_LO / SCORE_PCTILE_HI percentiles, mean over the archetype's ortholog list.
Human .X is already log-normalized (ln(CPM+1)).

Reads:
  results/1_mouse_human/02.mouse_l23_archetype_markers.tsv
  data/input/human_mouse_orthologs.tsv
  data/input/jorstad23_human_WithinArea_L23IT.h5ad
Outputs:
  results/1_mouse_human/03.mouse_l23_archetype_scores_on_human_cells.tsv
"""

import os
import numpy as np
import pandas as pd
import anndata as ad

from l23.crossspecies import load_orthologs
from l23.expr import minmax_score

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

RES_DIR       = os.path.join(REPO_ROOT, 'results', '1_mouse_human')
IN_MARKERS    = os.path.join(RES_DIR, '02.mouse_l23_archetype_markers.tsv')
IN_ORTHOLOGS  = os.path.join(REPO_ROOT, 'data', 'input', 'human_mouse_orthologs.tsv')
INPUT_HUMAN   = os.path.join(REPO_ROOT, 'data', 'input', 'jorstad23_human_WithinArea_L23IT.h5ad')
OUT_SCORES    = os.path.join(RES_DIR, '03.mouse_l23_archetype_scores_on_human_cells.tsv')

NOC             = 3          # as script 02
GENE_NAME_COL   = 'feature_name'
SCORE_PCTILE_LO = 2
SCORE_PCTILE_HI = 98
ALPHABET        = ['A', 'B', 'C', 'D', 'E', 'F']

ortho = load_orthologs(IN_ORTHOLOGS)
mouse_to_human = dict(zip(ortho['mouse_symbol'], ortho['human_symbol']))
print(f'{len(mouse_to_human)} 1-to-1 orthologs')
markers_df = pd.read_csv(IN_MARKERS, sep='\t')

adata = ad.read_h5ad(INPUT_HUMAN)
gene_names  = adata.var[GENE_NAME_COL].values
gene_to_idx = {g: i for i, g in enumerate(gene_names)}

archetype_names = ALPHABET[:NOC]
scores = np.zeros((adata.n_obs, NOC), dtype=np.float32)
for k, name in enumerate(archetype_names):
    mouse_genes = markers_df.loc[markers_df['archetype'] == f'archetype_{k+1}', 'gene'].values
    human_genes = [mouse_to_human[g] for g in mouse_genes
                   if g in mouse_to_human and mouse_to_human[g] in gene_to_idx]
    if not human_genes:
        raise ValueError(f'no human orthologs of the archetype_{k+1} markers')
    print(f'Score {name}: {len(human_genes)}/{len(mouse_genes)} mouse markers mapped')
    mat = adata.X[:, [gene_to_idx[g] for g in human_genes]].toarray().astype(np.float32)
    scores[:, k] = minmax_score(mat, SCORE_PCTILE_LO, SCORE_PCTILE_HI)

pd.DataFrame(scores, index=adata.obs_names.values,
             columns=[f'score_{n}' for n in archetype_names]).to_csv(OUT_SCORES, sep='\t')
print(f'Saved {OUT_SCORES}')
