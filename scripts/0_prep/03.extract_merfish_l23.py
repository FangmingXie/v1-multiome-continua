"""Extract the L2/3 MERFISH cells, NR P28 (1_mouse_human/06) and NR P21 (5_easifish/02).

Every later step that needs MERFISH L2/3 cells reads these extracts instead of the raw objects.
L2/3 = gated_pred_subclass == 'L2/3' of the MERFISH cells (modality == 'merfish') of the combined objects.
No depth filter here (the analysis scripts keep their own). uns and obsm are not copied (no later use).
  P28 NR  the L2/3 cells of the 4 NR sections, labelled in the P28 NR/DR object, with the expression of the
          section's own 500-gene file (X and layers jnorm, norm; obs incl. depth and width): sections in
          P28_SECTIONS order, cells in section-file order. obs_names = the combined object's (sample + section
          cell id); columns sample and gated_pred_subclass added. depth must equal the combined object's.
  P21 NR  the L2/3 cells of the P21 NR sections (sample prefix P21NR) of the P8/P14/P21 NR object, in the
          source order; X, every layer, obs and var (361 genes).

Reads:
  data/raw/P28NRDR_v1glut_rna_merfish_240729.h5ad                               (L2/3 labels)
  data/raw/{P28NRa_ant,P28NRa_pos,P28NRb_ant,P28NRb_pos}_l2_v1_240723.h5ad      (section expression)
  data/raw/P8P14P21NR_v1l23glut_rna_merfish_250411.h5ad
Outputs:
  results/0_prep/03.merfish_p28nr_l23.h5ad
  results/0_prep/03.merfish_p21nr_l23.h5ad
"""

import os
import numpy as np
import anndata as ad

from l23.constants import L23_SUBCLASS

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_P28        = os.path.join(REPO_ROOT, 'data', 'raw', 'P28NRDR_v1glut_rna_merfish_240729.h5ad')
INPUT_P28_SECTION = os.path.join(REPO_ROOT, 'data', 'raw', '{}_l2_v1_240723.h5ad')
INPUT_P21        = os.path.join(REPO_ROOT, 'data', 'raw', 'P8P14P21NR_v1l23glut_rna_merfish_250411.h5ad')
OUT_P28          = os.path.join(REPO_ROOT, 'results', '0_prep', '03.merfish_p28nr_l23.h5ad')
OUT_P21          = os.path.join(REPO_ROOT, 'results', '0_prep', '03.merfish_p21nr_l23.h5ad')

P28_SECTIONS = ['P28NRa_ant', 'P28NRa_pos', 'P28NRb_ant', 'P28NRb_pos']   # every P28 NR section
P21_PREFIX   = 'P21NR'                                                   # sample-name prefix of the P21 NR sections
SUBCLASS_COL = 'gated_pred_subclass'

os.makedirs(os.path.dirname(OUT_P28), exist_ok=True)


def l23_merfish(obs):
    return (obs['modality'] == 'merfish') & (obs[SUBCLASS_COL] == L23_SUBCLASS)


# --- P28 NR: labels from the combined object, expression from each section's file ---
mobs = ad.read_h5ad(INPUT_P28, backed='r').obs
parts = []
for sec in P28_SECTIONS:
    lab = mobs[l23_merfish(mobs) & (mobs['sample'].astype(str) == sec)]
    if lab.empty:
        raise ValueError(f'no {L23_SUBCLASS} MERFISH cells of {sec} in {INPUT_P28}')
    path = INPUT_P28_SECTION.format(sec)
    sd = ad.read_h5ad(path, backed='r')
    rows = sd.obs_names.get_indexer(lab.index.str.slice(len(sec)))
    if (rows < 0).any():
        raise ValueError(f'{(rows < 0).sum()} {sec} {L23_SUBCLASS} cells not in {path}')
    order = np.argsort(rows)
    rows, lab = rows[order], lab.iloc[order]
    part = ad.AnnData(X=np.asarray(sd.X[rows]), obs=sd.obs.iloc[rows].copy(), var=sd.var.copy(),
                      layers={k: np.asarray(sd.layers[k][rows]) for k in sd.layers})
    if not np.array_equal(part.obs['depth'].to_numpy(), lab['depth'].to_numpy()):
        raise ValueError(f'{sec}: section depth differs from {INPUT_P28}')
    part.obs_names = lab.index
    part.obs['sample'] = sec
    part.obs[SUBCLASS_COL] = lab[SUBCLASS_COL].astype(str).to_numpy()
    print(f'P28 {sec}: {part.n_obs} {L23_SUBCLASS} cells of {sd.n_obs}', flush=True)
    parts.append(part)
if not all(p.var_names.equals(parts[0].var_names) for p in parts):
    raise ValueError('the P28 section files have different genes')
p28 = ad.concat(parts, merge='same')
p28.write_h5ad(OUT_P28)
print(f'wrote {OUT_P28} {p28.shape}', flush=True)

# --- P21 NR ---
src = ad.read_h5ad(INPUT_P21, backed='r')
o = src.obs
idx = np.flatnonzero((l23_merfish(o) & o['sample'].astype(str).str.startswith(P21_PREFIX)).to_numpy())
if len(idx) == 0:
    raise ValueError(f'no {L23_SUBCLASS} MERFISH cells with sample prefix {P21_PREFIX} in {INPUT_P21}')
p21 = src[idx].to_memory()
p21 = ad.AnnData(X=p21.X, obs=p21.obs, var=p21.var, layers=dict(p21.layers))
print(f'P21: {p21.n_obs} {L23_SUBCLASS} cells, samples {sorted(p21.obs["sample"].astype(str).unique())}')
p21.write_h5ad(OUT_P21)
print(f'wrote {OUT_P21} {p21.shape}')
