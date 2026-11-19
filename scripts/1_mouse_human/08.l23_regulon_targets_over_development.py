"""AP-1 and EGR regulon-target activity over normal development in L2/3 neurons, human and mouse.

Score per cell (both species; l23.regulon.regulon_score): each +/+ target's log2(CP10k+1) divided by
its NORM_MAX_PCTILE percentile over the species' cells (no baseline subtraction, so no expression stays
0), capped at 1, averaged over the regulon's targets present in the data. Per stage / age: mean +/- SEM
over cells, divided by the largest per-stage mean (1 = peak stage, 0 = no activity; relative
magnitudes across stages are preserved).

  human  Wang25 developing-cortex 10x multiome GEX, V1 L2/3 excitatory neurons (EN-L2_3-IT); the
         activating (+/+) eRegulons of 07 for HUMAN_PANELS; cells pooled by developmental stage
         (second trimester -> third trimester -> infancy -> adolescence).
  mouse  Yoo25 L2/3, normally reared (NR) ages P6-P21 (0_prep/01); the AP-1 + EGR members of
         3_regulon/03's B' TFs (developmental archetype genes, star rule) that have a +/+ regulon in
         both species (human has no FOSB regulon); +/+ targets from 3_regulon/01. Every B' TF of
         3_regulon/03 must land in exactly one of FAMILIES.
Figure (mouse vs human, side by side): 11. NR vs DR in mouse is 4_dr/01.

Reads:
  data/raw/wang25_human_multiome_gex_v1_l23.h5ad
  results/1_mouse_human/07.human_wang25_regulon_targets.tsv
  results/0_prep/01.yoo25_l23.h5ad
  results/3_regulon/01.yoo25_l23_regulon_targets.tsv
  results/3_regulon/03.enriched_regulon_selection.tsv
Outputs:
  results/1_mouse_human/08.wang25_l23_ap1_egr_regulon_targets_over_stage.tsv
  results/1_mouse_human/08.l23_b_regulon_targets_nr_mouse.tsv
"""

import os

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

from l23.expr import log2_cp10k, pctile_score
from l23.regulon import NORM_MAX_PCTILE, load_expression, regulon_score

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES_DIR   = os.path.join(REPO_ROOT, 'results', '1_mouse_human')

IN_HUMAN_H5AD    = os.path.join(REPO_ROOT, 'data', 'raw', 'wang25_human_multiome_gex_v1_l23.h5ad')
IN_HUMAN_REG     = os.path.join(RES_DIR, '07.human_wang25_regulon_targets.tsv')
IN_MOUSE_H5AD    = os.path.join(REPO_ROOT, 'results', '0_prep', '01.yoo25_l23.h5ad')
IN_MOUSE_REG     = os.path.join(REPO_ROOT, 'results', '3_regulon', '01.yoo25_l23_regulon_targets.tsv')
IN_SELECTION     = os.path.join(REPO_ROOT, 'results', '3_regulon', '03.enriched_regulon_selection.tsv')
OUT_HUMAN        = os.path.join(RES_DIR, '08.wang25_l23_ap1_egr_regulon_targets_over_stage.tsv')
OUT_MOUSE        = os.path.join(RES_DIR, '08.l23_b_regulon_targets_nr_mouse.tsv')

SIGN = '+/+'                   # activating regulons only
# human: panels (label, TFs)
HUMAN_PANELS = [('AP-1', ['FOS', 'FOSL2', 'JUN', 'JUNB', 'JUND']),
                ('EGR',  ['EGR1', 'EGR2', 'EGR3', 'EGR4'])]
KEEP_TYPE = 'EN-L2_3-IT'       # L2/3 excitatory neurons
TYPE_COL  = 'Type'
GROUP_COL = 'Group'
STAGE_ORDER = ['Second_trimester', 'Third_trimester', 'Infancy', 'Adolescence']   # temporal order
# mouse: B' TFs by family
ARCH = "B'"
FAMILIES = [('AP-1', ['Fos', 'Fosb', 'Fosl2', 'Junb']),
            ('EGR', ['Egr1', 'Egr2', 'Egr3', 'Egr4']),
            ('other', ['Smad3', 'Irf2', 'Mxi1', 'Etv5', 'Tfdp2'])]
MOUSE_FAMILIES = ['AP-1', 'EGR']

os.makedirs(RES_DIR, exist_ok=True)


def peak_scale(st, tf):
    """sem from std / n_cells; mean and sem divided by the peak (largest) mean."""
    st['sem'] = st['std'] / np.sqrt(st['n_cells'])
    peak = st['mean'].max()
    if not peak > 0:
        raise ValueError(f'{tf}: peak mean is {peak}; cannot max-scale')
    st['mean_norm'], st['sem_norm'] = st['mean'] / peak, st['sem'] / peak
    return st


def human():
    reg = pd.read_csv(IN_HUMAN_REG, sep='\t')
    reg = reg[reg['regulation_direction'] == SIGN]
    tf_targets = {tf: sorted(set(reg.loc[reg['regulon'] == f'{tf}_{SIGN}', 'Gene']))
                  for _, tfs in HUMAN_PANELS for tf in tfs}

    print(f'Loading {IN_HUMAN_H5AD}')
    adata = ad.read_h5ad(IN_HUMAN_H5AD)
    adata = adata[(adata.obs[TYPE_COL] == KEEP_TYPE).values].copy()   # the dataset is this type already
    if adata.n_obs == 0:
        raise ValueError(f'No {KEEP_TYPE} cells found')
    stage = adata.obs[GROUP_COL].astype(str).values
    unknown = set(stage) - set(STAGE_ORDER)
    if unknown:
        raise ValueError(f'Unknown developmental stage(s) not in STAGE_ORDER: {unknown}')
    stages = [s for s in STAGE_ORDER if s in set(stage)]
    print(f'{adata.n_obs} {KEEP_TYPE} cells across {len(stages)} stages: {stages}')

    # .X holds raw counts; library size = total counts over all genes
    X = sp.csr_matrix(adata.X).astype(np.float64)
    depth = np.asarray(X.sum(axis=1)).ravel()
    if np.any(depth == 0):
        raise ValueError('Found cells with zero total counts; cannot normalize')

    out = []
    for label, tfs in HUMAN_PANELS:
        for tf in tfs:
            used = [g for g in tf_targets[tf] if g in adata.var_names]
            if not used:
                print(f'  [{label}] {tf}: no {SIGN} regulon / no targets — skipped')
                continue
            x = np.asfortranarray(X[:, adata.var_names.get_indexer(used)].toarray())   # (sum order of the mean)
            score = pctile_score(log2_cp10k(x, depth), NORM_MAX_PCTILE)
            n_used = len(used)
            grp = pd.DataFrame({'stage': stage, 'score': score}).groupby('stage', observed=True)['score']
            st = (pd.DataFrame({'n_cells': grp.size(), 'mean': grp.mean(), 'std': grp.std(ddof=1)})
                  .reindex(stages).reset_index().rename(columns={'index': 'stage'}))
            st = peak_scale(st, tf)
            st.insert(0, 'TF', tf)
            st.insert(1, 'n_targets_used', n_used)
            st.insert(0, 'panel', label)
            out.append(st)
            print(f'  [{label}] {tf}: {n_used} targets used')
    out = pd.concat(out, ignore_index=True)
    out.to_csv(OUT_HUMAN, sep='\t', index=False)
    print(f'Wrote {OUT_HUMAN}')
    return out


def mouse(human_tfs):
    sel = pd.read_csv(IN_SELECTION, sep='\t')
    arch_tfs = set(sel.loc[sel['arch_letter'] == ARCH, 'TF'])
    fam_tfs = [tf for _, fam in FAMILIES for tf in fam]
    if sorted(fam_tfs) != sorted(arch_tfs):
        raise ValueError(f'FAMILIES {sorted(fam_tfs)} != {ARCH} TFs of {IN_SELECTION} {sorted(arch_tfs)}')
    family = {tf: fam for fam, tfs in FAMILIES for tf in tfs}
    reg = pd.read_csv(IN_MOUSE_REG, sep='\t')
    reg = reg[reg['regulation_direction'] == SIGN]
    candidates = [tf for fam, tfs in FAMILIES if fam in MOUSE_FAMILIES for tf in tfs]
    tfs = [tf for tf in candidates if tf.upper() in human_tfs]
    print(f'mouse vs human: {tfs}; mouse-only (no human regulon): {[tf for tf in candidates if tf not in tfs]}')
    targets = {tf: sorted(set(reg.loc[reg['TF'] == tf, 'Gene'])) for tf in tfs}
    empty = [tf for tf, t in targets.items() if not t]
    if empty:
        raise ValueError(f'no {SIGN} regulon in {IN_MOUSE_REG} for {empty}')

    L, cells = load_expression(IN_MOUSE_H5AD, sorted({g for t in targets.values() for g in t}))
    nr = (cells['condition'] == 'NR').to_numpy()
    L, cells = L[nr].reset_index(drop=True), cells[nr].reset_index(drop=True)
    rows = []
    for tf in tfs:
        score, n_used = regulon_score(tf, targets[tf], L)
        st = (cells.assign(score=score).groupby(['Age', 'day'])['score']
                   .agg(n_cells='size', mean='mean', std='std').reset_index().sort_values('day'))
        st = peak_scale(st, tf)
        st.insert(0, 'TF', tf)
        st.insert(1, 'n_targets_used', n_used)
        st['family'] = family[tf]
        rows.append(st)
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(OUT_MOUSE, sep='\t', index=False)
    print(f'wrote {OUT_MOUSE}')


if __name__ == '__main__':
    mouse(set(human()['TF']))
