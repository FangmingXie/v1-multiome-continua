"""Mouse vs human L2/3 axes at top-4000 HVGs: closed-form CCA on ortholog-paired Gate-A loadings.

Script 05 on the top-4000-HVG record (07 human, 08 mouse, extended to every expressed gene by 09),
over the gene universe `hvg_union`: 1-to-1 ortholog pairs expressed in both datasets that are an
HVG in either species (n = 5485). Statistics as in 05 (see l23.crossspecies). The canonical r falls
as the universe grows, so it is not comparable with 05's r — always quote it with n.
Gate A is re-derived, not inherited: varimax orders components by score variance and that order
changes with the gene set. HUMAN_VX / MOUSE_VX = the components whose largest variance-partition
factor is cell type in 07 / 08's post-varimax tables (human: VX3 .725, VX6 .366, VX8 .246, VX9 .287,
VX10 .196; dropped VX1 library size, VX2/VX4/VX5 donor, VX7 source. Mouse: VX1 .611, VX3 .576,
VX4 .180, VX6 .170, VX9 .323, VX10 .156; dropped VX2/VX5/VX8 sample, VX7 library size).

Reads:
  results/1_mouse_human/07.human_l23_varimax_loadings_hvg4000.tsv       (HVG membership)
  results/1_mouse_human/08.mouse_l23_varimax_loadings_hvg4000.tsv       (HVG membership)
  results/1_mouse_human/09.{human,mouse}_l23_varimax_loadings_full_hvg4000.tsv
  data/human_mouse_orthologs.tsv
Outputs:
  results/1_mouse_human/10.l23_axis_cca_spectrum_hvg4000_union.tsv
  results/1_mouse_human/10.l23_axis_cca_weights_human_hvg4000_union.tsv
  results/1_mouse_human/10.l23_axis_cca_weights_mouse_hvg4000_union.tsv
"""

import os
import pandas as pd

from l23.crossspecies import load_orthologs, run_cca, write_spectrum, write_weights

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

RES_DIR      = os.path.join(REPO_ROOT, 'results', '1_mouse_human')
IN_H_HVG     = os.path.join(RES_DIR, '07.human_l23_varimax_loadings_hvg4000.tsv')
IN_M_HVG     = os.path.join(RES_DIR, '08.mouse_l23_varimax_loadings_hvg4000.tsv')
IN_H_FULL    = os.path.join(RES_DIR, '09.human_l23_varimax_loadings_full_hvg4000.tsv')
IN_M_FULL    = os.path.join(RES_DIR, '09.mouse_l23_varimax_loadings_full_hvg4000.tsv')
IN_ORTHOLOGS = os.path.join(REPO_ROOT, 'data', 'human_mouse_orthologs.tsv')
OUT_SPECTRUM = os.path.join(RES_DIR, '10.l23_axis_cca_spectrum_hvg4000_union.tsv')
OUT_W_HUMAN  = os.path.join(RES_DIR, '10.l23_axis_cca_weights_human_hvg4000_union.tsv')
OUT_W_MOUSE  = os.path.join(RES_DIR, '10.l23_axis_cca_weights_mouse_hvg4000_union.tsv')

HUMAN_VX = ['VX3', 'VX6', 'VX8', 'VX9', 'VX10']
MOUSE_VX = ['VX1', 'VX3', 'VX4', 'VX6', 'VX9', 'VX10']

hvg_h = pd.read_csv(IN_H_HVG, sep='\t', index_col=0)
hvg_m = pd.read_csv(IN_M_HVG, sep='\t', index_col=0)
human = pd.read_csv(IN_H_FULL, sep='\t', index_col=0)
mouse = pd.read_csv(IN_M_FULL, sep='\t', index_col=0)
ortho = load_orthologs(IN_ORTHOLOGS)
shared = ortho[ortho['human_symbol'].isin(human.index) & ortho['mouse_symbol'].isin(mouse.index)]
shared = shared[shared['human_symbol'].isin(hvg_h.index)
                | shared['mouse_symbol'].isin(hvg_m.index)].reset_index(drop=True)
X = human.loc[shared['human_symbol'].values, HUMAN_VX].values
Y = mouse.loc[shared['mouse_symbol'].values, MOUSE_VX].values

print(f'L2/3 CCA at N_HVG=4000, universe hvg_union ({len(shared)} ortholog pairs)')
res, A, B, boot_h, boot_m = run_cca(X, Y)
write_spectrum(res, OUT_SPECTRUM)
write_weights(res, A, B, boot_h, boot_m, HUMAN_VX, MOUSE_VX, OUT_W_HUMAN, OUT_W_MOUSE)
