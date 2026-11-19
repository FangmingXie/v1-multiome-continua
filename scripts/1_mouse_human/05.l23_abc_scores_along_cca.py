"""Mouse L2/3 archetype scores (A', B', C') along the conserved CCA1 and CCA2 axes, human vs mouse: binned profiles.

A cell's CCA coordinate is `(VX coords - mean) . w` with the Gate-A canonical weights of script 04;
human cells carry the mouse-derived scores (script 03), mouse cells the native scores (script 02).
Letters: 02's A/B/C are relabelled to the published primed letters (RELABEL). CCA axes are in their
raw orientation. Panels (PANELS): A' along CCA1, B' along CCA2, C' along CCA1.

Profiles normalized for the cross-species comparison, per species:
  - x = CCA percentile within species (bin k is the same slice of each species' distribution);
  - y = archetype score z-scored over all cells of that species (removes the level / spread gap of
    human scores imputed from mouse), mean +/- SEM in N_BINS equal-count bins;
  - per panel, Pearson r between the N_BINS paired human and mouse bin means, with a one-sided
    permutation p (scores shuffled across cells within each species, N_PERM times; one RNG stream
    over the panels in PANELS order).
Also written: every cell's CCA1, CCA2 and A', B', C' scores, and the mouse archetypes (02's PCHA
vertices in the VX space) projected the same way as the mouse cells (figure 09).
Figures: 07, 09.

Reads:
  results/1_mouse_human/01.{mouse,human}_l23_varimax_coords.tsv
  results/1_mouse_human/02.mouse_l23_archetype_scores.tsv
  results/1_mouse_human/02.mouse_l23_pcha_aa_vx.tsv
  results/1_mouse_human/03.mouse_l23_archetype_scores_on_human_cells.tsv
  results/1_mouse_human/04.l23_axis_cca_weights_{human,mouse}.tsv
Outputs:
  results/1_mouse_human/05.l23_abc_scores_along_cca_profiles.tsv   (axis, archetype, species, bin, pct, mean, sem)
  results/1_mouse_human/05.l23_abc_scores_along_cca_stats.tsv      (axis, archetype, r, p_perm, n_perm)
  results/1_mouse_human/05.l23_cca_cells.tsv                       (species, cell, CCA1, CCA2, score_A', score_B', score_C')
  results/1_mouse_human/05.l23_cca_mouse_archetypes.tsv            (archetype, CCA1, CCA2)
"""

import os
import numpy as np
import pandas as pd

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

RES_DIR      = os.path.join(REPO_ROOT, 'results', '1_mouse_human')
IN_M_COORDS  = os.path.join(RES_DIR, '01.mouse_l23_varimax_coords.tsv')
IN_H_COORDS  = os.path.join(RES_DIR, '01.human_l23_varimax_coords.tsv')
IN_M_SCORES  = os.path.join(RES_DIR, '02.mouse_l23_archetype_scores.tsv')
IN_M_AA_VX   = os.path.join(RES_DIR, '02.mouse_l23_pcha_aa_vx.tsv')
IN_H_SCORES  = os.path.join(RES_DIR, '03.mouse_l23_archetype_scores_on_human_cells.tsv')
IN_H_WEIGHTS = os.path.join(RES_DIR, '04.l23_axis_cca_weights_human.tsv')
IN_M_WEIGHTS = os.path.join(RES_DIR, '04.l23_axis_cca_weights_mouse.tsv')
OUT_PROFILES = os.path.join(RES_DIR, '05.l23_abc_scores_along_cca_profiles.tsv')
OUT_STATS    = os.path.join(RES_DIR, '05.l23_abc_scores_along_cca_stats.tsv')
OUT_CELLS    = os.path.join(RES_DIR, '05.l23_cca_cells.tsv')
OUT_M_AA     = os.path.join(RES_DIR, '05.l23_cca_mouse_archetypes.tsv')

HUMAN_VX = ['VX2', 'VX6', 'VX7', 'VX8', 'VX9', 'VX10']   # Gate A, as script 04
MOUSE_VX = ['VX1', 'VX2', 'VX5', 'VX7', 'VX8', 'VX9']
RELABEL  = {'C': "A'", 'B': "B'", 'A': "C'"}             # 02's letters -> published letters
CCA_AXES = ['CCA1', 'CCA2']
PANELS   = [('CCA1', "A'"), ('CCA1', "C'"), ('CCA2', "B'")]   # (axis, archetype); the permutation order
N_BINS   = 20                                            # equal-count bins along each CCA axis
N_PERM   = 1000                                          # permutations for the human-mouse r null
SEED     = 0

os.makedirs(RES_DIR, exist_ok=True)


def cca_coords(coords_path, vx_cols, weights_path):
    """(cell index, cells x 2 CCA coordinates, [CCA1 r, CCA2 r], projection of VX points) for one species."""
    vx_df   = pd.read_csv(coords_path, sep='\t', index_col=0)
    weights = pd.read_csv(weights_path, sep='\t', index_col=0)
    W = weights.loc[CCA_AXES, vx_cols].values.T
    C = vx_df[vx_cols].values
    project = lambda P: (P - C.mean(axis=0)) @ W
    return vx_df.index, project(C), weights.loc[CCA_AXES, 'canonical_r'].values, project


def bin_means(bins, y):
    return np.bincount(bins, weights=y, minlength=N_BINS) / np.bincount(bins, minlength=N_BINS)


m_scores = pd.read_csv(IN_M_SCORES, sep='\t', index_col=0)
h_scores = pd.read_csv(IN_H_SCORES, sep='\t', index_col=0)
if list(h_scores.columns) != list(m_scores.columns):
    raise ValueError(f'Score columns differ between human ({list(h_scores.columns)}) '
                     f'and mouse ({list(m_scores.columns)}).')
old_letter = {RELABEL[c.split('_')[1]]: c.split('_')[1] for c in m_scores.columns}

h_index, h_cca, r, _ = cca_coords(IN_H_COORDS, HUMAN_VX, IN_H_WEIGHTS)
m_index, m_cca, _, m_project = cca_coords(IN_M_COORDS, MOUSE_VX, IN_M_WEIGHTS)
if not h_scores.index.equals(h_index):
    raise ValueError(f'Human scores ({IN_H_SCORES}) and coords ({IN_H_COORDS}) cover different cells/order.')
if not m_scores.index.equals(m_index):
    raise ValueError(f'Mouse scores ({IN_M_SCORES}) and coords ({IN_M_COORDS}) cover different cells/order.')
species = [('human', h_cca, h_scores), ('mouse', m_cca, m_scores)]
cells = pd.concat([pd.DataFrame({'species': name, 'cell': scores.index, 'CCA1': cca[:, 0], 'CCA2': cca[:, 1],
                                 **{f'score_{d}': scores[f'score_{old_letter[d]}'].values for d in sorted(old_letter)}})
                   for name, cca, scores in species], ignore_index=True)
cells.to_csv(OUT_CELLS, sep='\t', index=False)
aa_vx = pd.read_csv(IN_M_AA_VX, sep='\t', index_col=0)        # archetype_k = 02's letter ALPHABET[k-1]
m_aa = pd.DataFrame(m_project(aa_vx[MOUSE_VX].values), columns=CCA_AXES)
m_aa.insert(0, 'archetype', [RELABEL['ABC'[int(k.split('_')[1]) - 1]] for k in aa_vx.index])
m_aa.to_csv(OUT_M_AA, sep='\t', index=False)
print(f'--- L2/3; human cells {len(h_cca)}, mouse cells {len(m_cca)}; CCA r = {r[0]:.3f}, {r[1]:.3f} ---')

rng = np.random.default_rng(SEED)
pct = (np.arange(N_BINS) + 0.5) * 100 / N_BINS                # bin centers on the percentile scale
profiles, stats = [], []
for axis, disp in PANELS:
    i = CCA_AXES.index(axis)
    prof, z_bins = {}, {}
    for name, cca, scores in species:
        y = scores[f'score_{old_letter[disp]}'].values
        z = (y - y.mean()) / y.std()
        bins = pd.qcut(cca[:, i], N_BINS, labels=False)
        g = pd.Series(z).groupby(bins)
        prof[name], z_bins[name] = g.mean().values, (bins, z)
        profiles.append(pd.DataFrame({'axis': axis, 'archetype': disp, 'species': name, 'bin': np.arange(N_BINS),
                                      'pct': pct, 'mean': prof[name], 'sem': g.sem().values}))
    (hb, hz), (mb, mz) = z_bins['human'], z_bins['mouse']
    r_obs = np.corrcoef(prof['human'], prof['mouse'])[0, 1]
    r_null = np.array([np.corrcoef(bin_means(hb, rng.permutation(hz)),
                                   bin_means(mb, rng.permutation(mz)))[0, 1]
                       for _ in range(N_PERM)])
    p_perm = (1 + np.sum(r_null >= r_obs)) / (1 + N_PERM)
    print(f'  {disp} along {axis}: human-mouse r = {r_obs:.3f}, permutation p = {p_perm:.4f}')
    stats.append({'axis': axis, 'archetype': disp, 'r': r_obs, 'p_perm': p_perm, 'n_perm': N_PERM})

pd.concat(profiles, ignore_index=True).to_csv(OUT_PROFILES, sep='\t', index=False)
pd.DataFrame(stats).to_csv(OUT_STATS, sep='\t', index=False)
print(f'Saved {OUT_PROFILES}\nSaved {OUT_STATS}\nSaved {OUT_CELLS}\nSaved {OUT_M_AA}')
