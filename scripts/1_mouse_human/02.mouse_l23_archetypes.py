"""Mouse L2/3 (Cheng22, P28 NR) archetypes: seeded PCHA on the Gate-A varimax subspace, markers and scores.

PCHA (NOC archetypes) on the Gate-A mouse VX coordinates of script 01 (VX dominated by cell type in 01's
post-varimax variance partition). Top N_TOP_CELLS cells per archetype by distance in VX space;
conservative one-vs-each Wilcoxon markers (worst case over the other archetypes: min log2FC, max p;
BH; kept if frac_in >= FRAC_IN_THRESH, fdr < FDR_THRESH, log2FC > LOG2FC_THRESH); per-cell [0, 1]
archetype score = mean over the archetype's markers of expression min-max clipped at the
SCORE_PCTILE_LO / SCORE_PCTILE_HI percentiles.
Expression is log2(CP10k + 1) from the .raw integer counts, restricted to 01's HVGs (the genes the
loadings were fit on). log2FC is a difference of log2 values (no division by ln 2).
Letters: archetype_k -> A, B, C (internal PCHA order); downstream relabel to the published primed
letters: C -> A', B -> B', A -> C' (archetype_3 -> A', archetype_2 -> B', archetype_1 -> C').

Reads:
  results/1_mouse_human/01.mouse_l23_varimax_coords.tsv
  results/1_mouse_human/01.mouse_l23_varimax_loadings.tsv
  data/raw/superdupermegaRNA_cheng22_IT_P28NR.h5ad
Outputs:
  results/1_mouse_human/02.mouse_l23_pcha_xp.tsv
  results/1_mouse_human/02.mouse_l23_pcha_aa.tsv
  results/1_mouse_human/02.mouse_l23_pcha_aa_vx.tsv      (the archetypes in the Gate-A VX space)
  results/1_mouse_human/02.mouse_l23_top_cells.tsv
  results/1_mouse_human/02.mouse_l23_archetype_markers.tsv
  results/1_mouse_human/02.mouse_l23_archetype_scores.tsv
"""

import os
import numpy as np
import pandas as pd
import anndata as ad
import scipy.stats
from statsmodels.stats.multitest import multipletests

from l23.archetype import SCA
from l23.crossspecies import mouse_log2_cp10k
from l23.expr import minmax_score

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

RES_DIR        = os.path.join(REPO_ROOT, 'results', '1_mouse_human')
IN_VX_COORDS   = os.path.join(RES_DIR, '01.mouse_l23_varimax_coords.tsv')
IN_VX_LOADINGS = os.path.join(RES_DIR, '01.mouse_l23_varimax_loadings.tsv')
INPUT_MOUSE    = os.path.join(REPO_ROOT, 'data', 'raw', 'superdupermegaRNA_cheng22_IT_P28NR.h5ad')
OUT_PCHA_XP    = os.path.join(RES_DIR, '02.mouse_l23_pcha_xp.tsv')
OUT_PCHA_AA    = os.path.join(RES_DIR, '02.mouse_l23_pcha_aa.tsv')
OUT_PCHA_AA_VX = os.path.join(RES_DIR, '02.mouse_l23_pcha_aa_vx.tsv')
OUT_TOP_CELLS  = os.path.join(RES_DIR, '02.mouse_l23_top_cells.tsv')
OUT_MARKERS    = os.path.join(RES_DIR, '02.mouse_l23_archetype_markers.tsv')
OUT_SCORES     = os.path.join(RES_DIR, '02.mouse_l23_archetype_scores.tsv')

MOUSE_SUBCLASS   = 'L2/3'
VX_COLS          = ['VX1', 'VX2', 'VX5', 'VX7', 'VX8', 'VX9']   # Gate A (01's vx variance partition)
NOC              = 3
N_TOP_CELLS      = 300
SUBCLASS_COL     = 'Subclass'
CLUSTER_COL      = 'Type'
PCHA_SEED        = 0
FRAC_IN_THRESH   = 0.25
FDR_THRESH       = 0.001
LOG2FC_THRESH    = np.log2(1.5)
SCORE_PCTILE_LO  = 2
SCORE_PCTILE_HI  = 98
TOP_CELL_FRAC_MAX = 0.25             # NOC * N_TOP_CELLS must stay under this fraction of cells
ALPHABET         = ['A', 'B', 'C', 'D', 'E', 'F']

os.makedirs(RES_DIR, exist_ok=True)
ndim = len(VX_COLS) - 1
archetype_names = ALPHABET[:NOC]

vx_df = pd.read_csv(IN_VX_COORDS, sep='\t', index_col=0)
xn    = vx_df[VX_COLS].values
n_cells = len(vx_df)
if NOC * N_TOP_CELLS > TOP_CELL_FRAC_MAX * n_cells:
    raise ValueError(f'NOC*N_TOP_CELLS = {NOC * N_TOP_CELLS} exceeds {TOP_CELL_FRAC_MAX:.0%} of {n_cells} cells')

# --- PCHA; py_pcha draws its start from the global NumPy RNG and takes no seed argument ---
print(f'Fitting PCHA: NOC={NOC}, NDIM={ndim}, VX={VX_COLS}...')
sca = SCA(xn, vx_df[CLUSTER_COL].values)
sca.setup_feature_matrix(method='data')
np.random.seed(PCHA_SEED)
sca.proj_and_pcha(ndim, NOC)

pc_names = [f'PC{i+1}' for i in range(ndim)]
pd.DataFrame(sca.xp, index=vx_df.index, columns=pc_names).to_csv(OUT_PCHA_XP, sep='\t')
pd.DataFrame(sca.aa.T, columns=pc_names,
             index=[f'archetype_{k+1}' for k in range(NOC)]).to_csv(OUT_PCHA_AA, sep='\t')

# --- top cells per archetype, by distance in VX space ---
inner_comps = sca.pca_.components_[:ndim]                    # VX axes -> PCHA space
aa_vx = sca.aa.T @ inner_comps + sca.pca_.mean_              # (NOC, len(VX_COLS))
pd.DataFrame(aa_vx, columns=VX_COLS,
             index=[f'archetype_{k+1}' for k in range(NOC)]).to_csv(OUT_PCHA_AA_VX, sep='\t')
dists = np.stack([np.linalg.norm(xn - aa_vx[k], axis=1) for k in range(NOC)], axis=1)
top_cells = [np.argsort(dists[:, k])[:N_TOP_CELLS] for k in range(NOC)]
for k in range(NOC):
    for j in range(k + 1, NOC):
        if len(np.intersect1d(top_cells[k], top_cells[j])):
            raise ValueError(f'top-cell sets of archetype_{k+1} and archetype_{j+1} overlap; N_TOP_CELLS too large')
cell_index = vx_df.index.values
pd.DataFrame([{'archetype': f'archetype_{k+1}', 'rank': r, 'cell': cell_index[i]}
              for k in range(NOC) for r, i in enumerate(top_cells[k])]).to_csv(OUT_TOP_CELLS, sep='\t', index=False)

# --- expression on 01's HVGs ---
print('Loading expression...')
m_adata = ad.read_h5ad(INPUT_MOUSE)
m_adata = m_adata[m_adata.obs[SUBCLASS_COL] == MOUSE_SUBCLASS]
m_adata = m_adata[vx_df.index]        # align to the varimax coord order
X_norm = mouse_log2_cp10k(m_adata)
hvg_genes  = pd.read_csv(IN_VX_LOADINGS, sep='\t', index_col=0).index.values
gene_names = m_adata.var_names.values
hvg_idx    = np.where(np.isin(gene_names, hvg_genes))[0]
X_norm     = X_norm[:, hvg_idx]
hvg_names  = gene_names[hvg_idx]
n_hvg = X_norm.shape[1]
print(f'X_norm shape: {X_norm.shape}')

# --- conservative one-vs-each Wilcoxon per archetype ---
all_markers = []
for k in range(NOC):
    X_in = X_norm[top_cells[k]]
    log2fc_mat, pval_mat, frac_out_list = [], [], []
    for j in range(NOC):
        if j == k:
            continue
        X_out = X_norm[top_cells[j]]
        log2fc_mat.append(X_in.mean(axis=0) - X_out.mean(axis=0))
        pval_mat.append(np.array([scipy.stats.ranksums(X_in[:, g], X_out[:, g])[1] for g in range(n_hvg)]))
        frac_out_list.append((X_out > 0).mean(axis=0))
    pvals = np.stack(pval_mat).max(axis=0)
    df = pd.DataFrame({'gene': hvg_names, 'archetype': f'archetype_{k+1}',
                       'log2FC': np.stack(log2fc_mat).min(axis=0), 'pval': pvals,
                       'fdr': multipletests(pvals, method='fdr_bh')[1],
                       'frac_in': (X_in > 0).mean(axis=0), 'frac_out': np.stack(frac_out_list).mean(axis=0)})
    df = df[(df['frac_in'] >= FRAC_IN_THRESH) & (df['fdr'] < FDR_THRESH) & (df['log2FC'] > LOG2FC_THRESH)]
    all_markers.append(df.sort_values('log2FC', ascending=False))
    print(f'  archetype_{k+1}: {len(df)} markers')
    if df.empty:
        raise ValueError(f'archetype_{k+1} has no markers')
markers_df = pd.concat(all_markers, ignore_index=True)
markers_df.to_csv(OUT_MARKERS, sep='\t', index=False)

# --- archetype scores ---
scores = np.zeros((n_cells, NOC), dtype=np.float32)
for k, name in enumerate(archetype_names):
    genes = markers_df.loc[markers_df['archetype'] == f'archetype_{k+1}', 'gene'].values
    mat = X_norm[:, [np.where(hvg_names == g)[0][0] for g in genes]]
    scores[:, k] = minmax_score(mat, SCORE_PCTILE_LO, SCORE_PCTILE_HI)
    print(f'  Score {name}: {len(genes)} genes')
pd.DataFrame(scores, index=m_adata.obs_names.values,
             columns=[f'score_{n}' for n in archetype_names]).to_csv(OUT_SCORES, sep='\t')
for p in (OUT_PCHA_XP, OUT_PCHA_AA, OUT_PCHA_AA_VX, OUT_TOP_CELLS, OUT_MARKERS, OUT_SCORES):
    print(f'Saved {p}')
