"""Mouse vs human L2/3 axes: closed-form CCA on ortholog-paired Gate-A varimax loadings (top-2000 HVGs).

Genes are observations, Gate-A VX components the features: the human block X = 02's loadings on
HUMAN_VX, the mouse block Y = 01's loadings on MOUSE_VX, rows paired by 1-to-1 orthology over the
genes that are HVGs in both species. Statistics (spectrum, permutation null, module-blocked CV,
bootstrap weight stability) are described in l23.crossspecies. Gate A = the VX components whose
largest variance-partition factor is cell type (02 / 01 post-varimax tables); the mouse set is the
PCHA subspace of script 03.

Reads:
  results/1_mouse_human/02.human_l23_varimax_loadings.tsv
  results/1_mouse_human/01.mouse_l23_varimax_loadings.tsv
  data/raw/human_mouse_orthologs.tsv
Outputs:
  results/1_mouse_human/05.l23_axis_cca_spectrum.tsv
  results/1_mouse_human/05.l23_axis_cca_weights_human.tsv
  results/1_mouse_human/05.l23_axis_cca_weights_mouse.tsv
"""

import os
import pandas as pd

from l23.crossspecies import load_orthologs, run_cca, write_spectrum, write_weights

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

RES_DIR       = os.path.join(REPO_ROOT, 'results', '1_mouse_human')
IN_H_LOADINGS = os.path.join(RES_DIR, '02.human_l23_varimax_loadings.tsv')
IN_M_LOADINGS = os.path.join(RES_DIR, '01.mouse_l23_varimax_loadings.tsv')
IN_ORTHOLOGS  = os.path.join(REPO_ROOT, 'data', 'raw', 'human_mouse_orthologs.tsv')
OUT_SPECTRUM  = os.path.join(RES_DIR, '05.l23_axis_cca_spectrum.tsv')
OUT_W_HUMAN   = os.path.join(RES_DIR, '05.l23_axis_cca_weights_human.tsv')
OUT_W_MOUSE   = os.path.join(RES_DIR, '05.l23_axis_cca_weights_mouse.tsv')

HUMAN_VX = ['VX2', 'VX6', 'VX7', 'VX8', 'VX9', 'VX10']   # Gate A, 02's vx variance partition
MOUSE_VX = ['VX1', 'VX2', 'VX5', 'VX7', 'VX8', 'VX9']    # Gate A, 01's vx variance partition (= 03)

human = pd.read_csv(IN_H_LOADINGS, sep='\t', index_col=0)
mouse = pd.read_csv(IN_M_LOADINGS, sep='\t', index_col=0)
ortho = load_orthologs(IN_ORTHOLOGS)
shared = ortho[ortho['human_symbol'].isin(human.index) & ortho['mouse_symbol'].isin(mouse.index)].reset_index(drop=True)
X = human.loc[shared['human_symbol'].values, HUMAN_VX].values
Y = mouse.loc[shared['mouse_symbol'].values, MOUSE_VX].values

print(f'L2/3 CCA, HVGs in both species ({len(shared)} ortholog pairs)')
res, A, B, boot_h, boot_m = run_cca(X, Y)
write_spectrum(res, OUT_SPECTRUM)
write_weights(res, A, B, boot_h, boot_m, HUMAN_VX, MOUSE_VX, OUT_W_HUMAN, OUT_W_MOUSE)
