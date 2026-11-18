"""Human L2/3 IT (Jorstad23): PCA + varimax rotation + variance partitioning, top-2000 HVGs.

.X is already log-normalized (ln(CPM+1)) and is used directly. Gene identifiers are ENSG in
var_names; symbols (var['feature_name']) index the loadings, because the ortholog table is
symbol-based. HVGs = top N_HVG by raw variance; StandardScaler -> PCA(N_PCS) -> Kaiser varimax; VX
components ordered by descending score variance. Variance partition over four factors (cell_type /
donor / source / library_size; development_stage is excluded: it is 1:1 with donor_id); the
post-varimax table is what the Gate-A VX selection (VX dominated by cell type) of script 05 is read from.

Reads:
  data/raw/jorstad23_human_WithinArea_L23IT.h5ad
Outputs:
  results/1_mouse_human/02.human_l23_varimax_coords.tsv
  results/1_mouse_human/02.human_l23_varimax_loadings.tsv
  results/1_mouse_human/02.human_l23_pc_variance_partition.tsv
  results/1_mouse_human/02.human_l23_vx_variance_partition.tsv
"""

import os
import numpy as np
import pandas as pd
import anndata as ad

from l23.crossspecies import pca_varimax, variance_partition, dummies, library_size_factor

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_HUMAN     = os.path.join(REPO_ROOT, 'data', 'raw', 'jorstad23_human_WithinArea_L23IT.h5ad')
OUT_DIR         = os.path.join(REPO_ROOT, 'results', '1_mouse_human')
OUT_VX_COORDS   = os.path.join(OUT_DIR, '02.human_l23_varimax_coords.tsv')
OUT_VX_LOADINGS = os.path.join(OUT_DIR, '02.human_l23_varimax_loadings.tsv')
OUT_PC_VAR_TSV  = os.path.join(OUT_DIR, '02.human_l23_pc_variance_partition.tsv')
OUT_VX_VAR_TSV  = os.path.join(OUT_DIR, '02.human_l23_vx_variance_partition.tsv')

N_HVG         = 2000
N_PCS         = 10
CLUSTER_COL   = 'WithinArea_cluster'
DONOR_COL     = 'donor_id'
SOURCE_COL    = 'Source'
LIBSIZE_COL   = 'nCount_RNA'
GENE_NAME_COL = 'feature_name'
META_COLS     = [CLUSTER_COL, DONOR_COL, SOURCE_COL, 'development_stage', LIBSIZE_COL]

os.makedirs(OUT_DIR, exist_ok=True)

print(f'Loading {INPUT_HUMAN}...')
adata = ad.read_h5ad(INPUT_HUMAN)
print(f'  {adata.n_obs} cells, {adata.n_vars} genes')
obs = adata.obs
factors = {'cell_type': dummies(obs[CLUSTER_COL], 'cell_type'), 'donor': dummies(obs[DONOR_COL], 'donor'),
           'source': dummies(obs[SOURCE_COL], 'source'),
           'library_size': library_size_factor(obs[LIBSIZE_COL].values.astype(float))}

X_norm = adata.X.toarray().astype(np.float32)
print(f'Top {N_HVG} HVGs -> PCA -> varimax...')
hvg_names, scores, vx_scores, vx_loadings = pca_varimax(X_norm, adata.var[GENE_NAME_COL].values, N_HVG, N_PCS)
del X_norm

pc_cols = [f'PC{i+1}' for i in range(N_PCS)]
vx_cols = [f'VX{i+1}' for i in range(N_PCS)]
variance_partition(scores, pc_cols, factors, 'pre-varimax PCs').to_csv(OUT_PC_VAR_TSV, sep='\t')
variance_partition(vx_scores, vx_cols, factors, 'post-varimax VX').to_csv(OUT_VX_VAR_TSV, sep='\t')

vx_df = pd.DataFrame(vx_scores, index=adata.obs_names, columns=vx_cols)
for col in META_COLS:
    vx_df[col] = obs[col].values
vx_df.to_csv(OUT_VX_COORDS, sep='\t')
pd.DataFrame(vx_loadings, index=hvg_names, columns=vx_cols).to_csv(OUT_VX_LOADINGS, sep='\t')
for p in (OUT_PC_VAR_TSV, OUT_VX_VAR_TSV, OUT_VX_COORDS, OUT_VX_LOADINGS):
    print(f'Saved {p}')
