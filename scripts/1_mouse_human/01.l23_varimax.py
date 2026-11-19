"""Mouse (Cheng22, P28 NR) and human (Jorstad23) L2/3 IT: PCA + varimax rotation + variance partitioning, top-2000 HVGs.

Per species: HVGs = top N_HVG by raw variance of the log-normalized matrix; StandardScaler -> PCA(N_PCS)
-> Kaiser varimax; VX components ordered by descending score variance (l23.crossspecies.pca_varimax).
The post-varimax variance partition is what the Gate-A VX selection (VX dominated by cell type) of
scripts 02 and 04 is read from.
  mouse  the integer raw counts of the 0_prep/01 Cheng22 L2/3 extract are normalized with
         log2(CP10k + 1). Factors: cell_type / sample / library_size.
  human  .X is already log-normalized (ln(CPM+1)) and is used directly. Gene identifiers are ENSG in
         var_names; symbols (var['feature_name']) index the loadings, because the ortholog table is
         symbol-based. Factors: cell_type / donor / source / library_size (development_stage is
         excluded: it is 1:1 with donor_id).

Reads:
  results/0_prep/01.cheng22_p28nr_l23.h5ad
  data/raw/jorstad23_human_WithinArea_L23IT.h5ad
Outputs (species = mouse, human):
  results/1_mouse_human/01.{species}_l23_varimax_coords.tsv
  results/1_mouse_human/01.{species}_l23_varimax_loadings.tsv
  results/1_mouse_human/01.{species}_l23_pc_variance_partition.tsv
  results/1_mouse_human/01.{species}_l23_vx_variance_partition.tsv
"""

import os
import numpy as np
import pandas as pd
import anndata as ad

from l23.crossspecies import pca_varimax, variance_partition, dummies, library_size_factor, mouse_log2_cp10k

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_MOUSE = os.path.join(REPO_ROOT, 'results', '0_prep', '01.cheng22_p28nr_l23.h5ad')
INPUT_HUMAN = os.path.join(REPO_ROOT, 'data', 'raw', 'jorstad23_human_WithinArea_L23IT.h5ad')
OUT_DIR     = os.path.join(REPO_ROOT, 'results', '1_mouse_human')
OUT         = {k: os.path.join(OUT_DIR, f'01.{{}}_l23_{k}.tsv') for k in
               ['varimax_coords', 'varimax_loadings', 'pc_variance_partition', 'vx_variance_partition']}

N_HVG = 2000
N_PCS = 10
# mouse
M_CLUSTER_COL  = 'Type'
M_SAMPLE_COL   = 'Sample'
M_LIBSIZE_COL  = 'n_counts'
# human
H_CLUSTER_COL   = 'WithinArea_cluster'
H_DONOR_COL     = 'donor_id'
H_SOURCE_COL    = 'Source'
H_LIBSIZE_COL   = 'nCount_RNA'
H_GENE_NAME_COL = 'feature_name'
H_META_COLS     = [H_CLUSTER_COL, H_DONOR_COL, H_SOURCE_COL, 'development_stage', H_LIBSIZE_COL]

os.makedirs(OUT_DIR, exist_ok=True)


def fit_and_save(species, X_norm, gene_names, factors, cells, meta):
    """pca_varimax + variance partitions; coords (VX scores + meta columns) and loadings written."""
    print(f'{species}: top {N_HVG} HVGs -> PCA -> varimax...')
    hvg_names, scores, vx_scores, vx_loadings = pca_varimax(X_norm, gene_names, N_HVG, N_PCS)
    pc_cols = [f'PC{i+1}' for i in range(N_PCS)]
    vx_cols = [f'VX{i+1}' for i in range(N_PCS)]
    out = {k: p.format(species) for k, p in OUT.items()}
    variance_partition(scores, pc_cols, factors, 'pre-varimax PCs').to_csv(out['pc_variance_partition'], sep='\t')
    variance_partition(vx_scores, vx_cols, factors, 'post-varimax VX').to_csv(out['vx_variance_partition'], sep='\t')
    vx_df = pd.DataFrame(vx_scores, index=cells, columns=vx_cols)
    for col, v in meta.items():
        vx_df[col] = v
    vx_df.to_csv(out['varimax_coords'], sep='\t')
    pd.DataFrame(vx_loadings, index=hvg_names, columns=vx_cols).to_csv(out['varimax_loadings'], sep='\t')
    for p in out.values():
        print(f'Saved {p}')


print('Loading mouse cheng22 data...')
m_adata = ad.read_h5ad(INPUT_MOUSE)
print(f'  {m_adata.n_obs} cells, {m_adata.n_vars} genes')
types, samples = m_adata.obs[M_CLUSTER_COL].values, m_adata.obs[M_SAMPLE_COL].values
fit_and_save('mouse', mouse_log2_cp10k(m_adata), m_adata.var_names.values,
             {'cell_type': dummies(types, 'cell_type'), 'sample': dummies(samples, 'sample'),
              'library_size': library_size_factor(m_adata.obs[M_LIBSIZE_COL].values)},
             m_adata.obs_names, {M_CLUSTER_COL: types, M_SAMPLE_COL: samples})
del m_adata

print(f'Loading {INPUT_HUMAN}...')
adata = ad.read_h5ad(INPUT_HUMAN)
print(f'  {adata.n_obs} cells, {adata.n_vars} genes')
obs = adata.obs
fit_and_save('human', adata.X.toarray().astype(np.float32), adata.var[H_GENE_NAME_COL].values,
             {'cell_type': dummies(obs[H_CLUSTER_COL], 'cell_type'), 'donor': dummies(obs[H_DONOR_COL], 'donor'),
              'source': dummies(obs[H_SOURCE_COL], 'source'),
              'library_size': library_size_factor(obs[H_LIBSIZE_COL].values.astype(float))},
             adata.obs_names, {col: obs[col].values for col in H_META_COLS})
