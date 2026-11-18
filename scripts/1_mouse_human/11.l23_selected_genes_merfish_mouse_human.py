"""Hand-picked conserved L2/3 A/C genes in situ: A vs C score, per-gene maps and depth profiles, mouse vs human MERFISH.

Genes: SELECTED (mouse symbol -> side; A = superficial, C = deep), mapped to human by 1-to-1 orthologs
(data/raw/human_mouse_orthologs.tsv). Every gene must be on both MERFISH panels (checked). The CCA1 z per
species (09 loadings x 10 CCA1 weights, z over all CCA ortholog pairs; CCA1 oriented so mouse A'
markers load positive and C' markers negative, checked) is carried only as a label; genes outside the
CCA universe get NaN. The summary-score figure uses mouse MOUSE_SECTION; the per-gene figure shows every NR mouse section (MOUSE_SECTIONS) plus the human section in
each gene's row, each panel scaled by the gene's GENE_PCTILE percentile within that panel.
Relative scale (gene and profile figures): per section, relative L2/3 depth = (depth - p2) / (p98 - p2) over
its kept L2/3 cells (0 = top, 1 = bottom of the L2/3 band), and position along the section in the same
unit (L2/3 thicknesses), so mouse (~300 um band) and human (~1500 um band) share one relative scale. In the
in situ figure every panel has the human panel's shape: human at equal aspect, mouse stretched along the
section (its aspect is not fixed) so both panels are the same size. Genes: SELECTED feeds the summary
score; INSITU_GENES are drawn in situ; PROFILE_GENES (a row per side) get depth profiles. Profile figure: per gene, mean log expression in N_REL_BINS equal-width
relative-depth bins (cells outside [0, 1] left out), divided by its highest bin (peak = 1), for every
PROFILE_MOUSE_SECTIONS and PROFILE_HUMAN_SECTIONS section (each human section loaded, rotated,
outlier-filtered and depth-normalized as SECTION); plotted as mean +- SEM over sections per species. Sections,
not cells, are the replicates: pooling cells first would weight sections by cell count and mix panels of
different sensitivity (250 vs 4000 genes), whereas per-section peak normalization removes those offsets.
Every section (mouse and human; all figures): before outlier removal, depth is re-measured from a 2nd-order
polynomial fit to the top of L2/3 (CURVE_PCTILE depth percentile in CURVE_N_BINS equal-count bins along the
section); check figure OUT_FIG_CURVE.
Human sections are kept only if every INSITU_GENES gene is detected in > H_MIN_DETECT of their L2/3 cells
(some 4000-gene sections barely detect KCNH5); this applies to the profile and all-sections figures.
In situ figure: axes stripped; dashed lines at relative depth 0 and 1 (top and bottom of L2/3).
All-sections in situ figure (to inspect depth alignment): rows = every profile section (mouse, then human),
columns = INSITU_GENES; all panels the same size, y = relative depth over ALLSEC_YLIM, x = position stretched
to the panel; dashed lines at 0 and 1. Note CDH13 is on the A side in mouse
(CCA1 z 2.2) but near 0 in human (0.3).
Cells and expression:
  mouse  P28 NR sections: MERFISH cells (modality == 'merfish') with gated_pred_subclass == 'L2/3' and depth
         from pia in MOUSE_DEPTH_RANGE (um; drops label-transfer outliers), from the P28 NR/DR MERFISH object;
         expression log2(jnorm + 1) from the section's 500-gene file; position = width x depth from pia (um).
  human  Fang22 SECTION (250 genes): eL2/3.IT cells; counts scaled to the median total per cell, log1p;
         position = pre-expansion adjusted.x/y (um), rotated so pia is up (superficial eL2/3.IT above deep cells).
  In each species, cells whose depth lies outside Tukey's fences (Q1 - TUKEY_K IQR, Q3 + TUKEY_K IQR)
  are removed before scoring.
Score (per species): mean over the side's genes of the gene scaled by its SCORE_MAX_PCTILE
percentile over the kept L2/3 cells, capped at 1; z-scored; A vs C = C z - A z. Rows = species; columns =
A z, C z, A vs C (C0 -> grey -> C2, colour bar flipped so the A side is on top); colour ranges clipped at
CLIP_PCTILE percentiles per species (A vs C keeps grey at 0). Depth axes start at the top L2/3 cell.

Reads:
  results/1_mouse_human/{07.human,08.mouse}_l23_varimax_loadings_hvg4000.tsv     (HVG membership)
  results/1_mouse_human/09.{human,mouse}_l23_varimax_loadings_full_hvg4000.tsv
  results/1_mouse_human/10.l23_axis_cca_weights_{human,mouse}_hvg4000_union.tsv
  results/3_regulon/03.l23_adult_archetype_markers.tsv                         (orientation check)
  data/raw/human_mouse_orthologs.tsv
  data/raw/P28NRDR_v1glut_rna_merfish_240729.h5ad                               (L2/3 cells, depth)
  data/raw/{P28NRa_ant,P28NRa_pos,P28NRb_ant,P28NRb_pos}_l2_v1_240723.h5ad      (section expression)
  data/raw/fang22_human_merfish/{SECTION}.{features,genes,matrix}.csv          (every PROFILE_HUMAN_SECTIONS)
Outputs:
  results/1_mouse_human/11.l23_selected_genes.tsv
  results/1_mouse_human/11.l23_selected_genes_ac_scores_cells.tsv
  results/1_mouse_human/11.l23_selected_genes_rel_depth_profiles.tsv
  figures/1_mouse_human/11.l23_selected_genes_ac_scores_mouse_human.pdf
  figures/1_mouse_human/11.l23_selected_genes_in_situ_mouse_human.pdf   (each selected ortholog, in situ)
  figures/1_mouse_human/11.l23_selected_genes_rel_depth_profiles.pdf   (depth profiles, mouse vs human)
  figures/1_mouse_human/11.l23_selected_genes_in_situ_all_sections.pdf (INSITU_GENES, every profile section)
  figures/1_mouse_human/11.l23_top_fit_all_sections.pdf                (every profile section: fitted top of L2/3)
"""

import os
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from scipy.sparse import coo_matrix
from scipy.stats import spearmanr
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

RES_DIR    = os.path.join(REPO_ROOT, 'results', '1_mouse_human')
FIG_DIR    = os.path.join(REPO_ROOT, 'figures', '1_mouse_human')
IN_H_HVG   = os.path.join(RES_DIR, '07.human_l23_varimax_loadings_hvg4000.tsv')
IN_M_HVG   = os.path.join(RES_DIR, '08.mouse_l23_varimax_loadings_hvg4000.tsv')
IN_H_FULL  = os.path.join(RES_DIR, '09.human_l23_varimax_loadings_full_hvg4000.tsv')
IN_M_FULL  = os.path.join(RES_DIR, '09.mouse_l23_varimax_loadings_full_hvg4000.tsv')
IN_W_HUMAN = os.path.join(RES_DIR, '10.l23_axis_cca_weights_human_hvg4000_union.tsv')
IN_W_MOUSE = os.path.join(RES_DIR, '10.l23_axis_cca_weights_mouse_hvg4000_union.tsv')
IN_MOUSE_MARKERS = os.path.join(REPO_ROOT, 'results', '3_regulon', '03.l23_adult_archetype_markers.tsv')
IN_ORTHO   = os.path.join(REPO_ROOT, 'data', 'raw', 'human_mouse_orthologs.tsv')
IN_MERFISH = os.path.join(REPO_ROOT, 'data', 'raw', 'P28NRDR_v1glut_rna_merfish_240729.h5ad')
MOUSE_SECTION  = 'P28NRb_pos'                          # summary-score figure
MOUSE_SECTIONS = ['P28NRb_pos']   # gene figure
MOUSE_SECTION_FILE = os.path.join(REPO_ROOT, 'data', 'raw', '{}_l2_v1_240723.h5ad')
IN_MOUSE_SECTION = MOUSE_SECTION_FILE.format(MOUSE_SECTION)
HUMAN_DIR  = os.path.join(REPO_ROOT, 'data', 'raw', 'fang22_human_merfish')
SECTION    = 'H18.06.006.MTG.250.expand.rep1'
IN_H_GENES    = os.path.join(HUMAN_DIR, f'{SECTION}.genes.csv')
OUT_GENES  = os.path.join(RES_DIR, '11.l23_selected_genes.tsv')
OUT_CELLS  = os.path.join(RES_DIR, '11.l23_selected_genes_ac_scores_cells.tsv')
OUT_FIG    = os.path.join(FIG_DIR, '11.l23_selected_genes_ac_scores_mouse_human.pdf')
OUT_PROFILES = os.path.join(RES_DIR, '11.l23_selected_genes_rel_depth_profiles.tsv')
OUT_FIG_PROFILE = os.path.join(FIG_DIR, '11.l23_selected_genes_rel_depth_profiles.pdf')
OUT_FIG_ALLSEC = os.path.join(FIG_DIR, '11.l23_selected_genes_in_situ_all_sections.pdf')
OUT_FIG_CURVE = os.path.join(FIG_DIR, '11.l23_top_fit_all_sections.pdf')
OUT_FIG_GENES = os.path.join(FIG_DIR, '11.l23_selected_genes_in_situ_mouse_human.pdf')

HUMAN_VX = ['VX3', 'VX6', 'VX8', 'VX9', 'VX10']           # Gate A at N_HVG = 4000, as 10
MOUSE_VX = ['VX1', 'VX3', 'VX4', 'VX6', 'VX9', 'VX10']
AXIS     = 'CCA1'
SELECTED = {'Cdh13': 'A', 'Robo1': 'A', 'Kcnh5': 'C', 'Tox': 'C'}   # summary score: mouse symbol -> side
INSITU_GENES = ['Cdh13', 'Kcnh5']   # in situ figures: one row (or column, all-sections figure) per gene
ALLSEC_YLIM  = (-0.35, 1.35)        # all-sections figure: same relative-depth window in every panel
H_MIN_DETECT = 0.10   # profile + all-sections figures: keep human sections where every INSITU_GENES gene is
                      # detected in > this fraction of the section's L2/3 cells
PROFILE_GENES = {   # depth-profile figure (row A, row C); all on both MERFISH panels
    'A': ['Cdh13', 'Robo1', 'Epha6'],
    'C': ['Kcnh5', 'Tox', 'Lrp1b'],
}
MOUSE_A_ARCH, MOUSE_C_ARCH = 'archetype_3', 'archetype_1'   # mouse L2/3 A', C' (C -> A', A -> C')
MOUSE_SUBCLASS = 'L2/3'
MOUSE_SUBCLASS_COL = 'gated_pred_subclass'
MOUSE_CONDITION = 'P28NR'           # sample-name prefix of the NR sections
MOUSE_DEPTH_RANGE = (50, 450)      # um from pia; ~ the 5-95% depth range of L2/3 MERFISH cells
MOUSE_LAYER    = 'jnorm'
H_SUBCLASS, H_SUBCLASS_COL = 'eL2/3.IT', 'cluster_L2'
H_DEEP = ['eL6.IT', 'eL6b', 'eL6.CT', 'eL6.IT.CAR3', 'eL5/6.NP', 'eL6.CAR3', 'eNP']   # orient the section
H_SRC_X, H_SRC_Y = 'adjusted.x', 'adjusted.y'   # pre-expansion um
H_MIN_ANCHOR = 5
PROFILE_MOUSE_SECTIONS = ['P28NRa_ant', 'P28NRa_pos', 'P28NRb_ant', 'P28NRb_pos']   # profile figure: all NR
PROFILE_HUMAN_SECTIONS = [   # profile figure: Fang22 human sections with >= 300 eL2/3.IT cells (all have the genes)
    'H18.06.006.MTG.250.expand.rep1', 'H18.06.006.MTG.4000.expand.rep1', 'H18.06.006.MTG.4000.expand.rep3',
    'H22.26.401.MTG.250.unexpand.rep1', 'H22.26.401.MTG.4000.expand.rep1',
    'H20.30.001.STG.250.unexpand.rep1', 'H20.30.001.STG.250.unexpand.rep2',
    'H20.30.001.STG.4000.expand.rep1', 'H20.30.001.STG.4000.expand.rep2',
]
CURVE_N_BINS  = 10    # top-of-L2/3 fit (every section): equal-count bins along the section ...
CURVE_PCTILE  = 5     # ... each giving (median x, this depth percentile); a 2nd-order polynomial fit to these
SCORE_MAX_PCTILE = 98
TUKEY_K     = 1.5
COLS = [('z_A', 'A side score (z)', 'viridis'), ('z_C', 'C side score (z)', 'viridis'),
        ('ac', 'A vs C: C z - A z', 'ac')]
CLIP_PCTILE = (5, 95)
CMAPS = {'ac': LinearSegmentedColormap.from_list('ac', ['C0', '0.8', 'C2'])}
DOT = 6
REL_PCTILE  = (2, 98)  # relative L2/3 depth: these depth percentiles of the section's kept L2/3 cells -> 0, 1
N_REL_BINS  = 5        # profile figure: equal-width bins over relative depth [0, 1]
SPECIES_COLORS = {'mouse': 'C3', 'human': 'k'}
GENE_PCTILE = 98       # gene figure: per gene and species, log expression / this percentile, capped at 1
GENE_CMAP   = sns.color_palette('rocket_r', as_cmap=True)   # gene figure

os.makedirs(os.path.dirname(OUT_GENES), exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

# --- per-gene CCA1 loadings ---
ortho = pd.read_csv(IN_ORTHO, sep='\t').drop_duplicates('human_symbol').drop_duplicates('mouse_symbol')
hvg_h = pd.read_csv(IN_H_HVG, sep='\t', index_col=0)
hvg_m = pd.read_csv(IN_M_HVG, sep='\t', index_col=0)
H = pd.read_csv(IN_H_FULL, sep='\t', index_col=0)
M = pd.read_csv(IN_M_FULL, sep='\t', index_col=0)
sh = ortho[ortho['human_symbol'].isin(H.index) & ortho['mouse_symbol'].isin(M.index)]
sh = sh[sh['human_symbol'].isin(hvg_h.index) | sh['mouse_symbol'].isin(hvg_m.index)].reset_index(drop=True)
X = H.loc[sh['human_symbol'], HUMAN_VX].to_numpy()
Y = M.loc[sh['mouse_symbol'], MOUSE_VX].to_numpy()
wh = pd.read_csv(IN_W_HUMAN, sep='\t', index_col=0).loc[AXIS, HUMAN_VX].to_numpy(dtype=float)
wm = pd.read_csv(IN_W_MOUSE, sep='\t', index_col=0).loc[AXIS, MOUSE_VX].to_numpy(dtype=float)
sh['human_load'] = (X - X.mean(0)) @ wh
sh['mouse_load'] = (Y - Y.mean(0)) @ wm
sh['human_z'] = (sh['human_load'] - sh['human_load'].mean()) / sh['human_load'].std()
sh['mouse_z'] = (sh['mouse_load'] - sh['mouse_load'].mean()) / sh['mouse_load'].std()

mk = pd.read_csv(IN_MOUSE_MARKERS, sep='\t')
med = {a: sh.loc[sh['mouse_symbol'].isin(mk.loc[mk['archetype'] == a, 'gene']), 'mouse_load'].median()
       for a in (MOUSE_A_ARCH, MOUSE_C_ARCH)}
if not med[MOUSE_A_ARCH] > 0 > med[MOUSE_C_ARCH]:
    raise ValueError(f"{AXIS} orientation: mouse A' / C' marker median loadings {med}; expected A' > 0 > C'")
print(f"{AXIS}: mouse A' markers median loading {med[MOUSE_A_ARCH]:.3f} (> 0), C' {med[MOUSE_C_ARCH]:.3f} (< 0)")

# --- panels; pick genes ---
msd = ad.read_h5ad(IN_MOUSE_SECTION, backed='r')
h_genes = pd.read_csv(IN_H_GENES)['name'].tolist()
SIDE = {**{g: s for s, gs in PROFILE_GENES.items() for g in gs}, **SELECTED}
for g in INSITU_GENES:
    if g not in SIDE:
        raise ValueError(f'INSITU_GENES {g} has no side (add it to SELECTED or PROFILE_GENES)')
conflict = {g for s, gs in PROFILE_GENES.items() for g in gs if g in SELECTED and SELECTED[g] != s}
if conflict:
    raise ValueError(f'side differs between SELECTED and PROFILE_GENES: {sorted(conflict)}')
sel = ortho[ortho['mouse_symbol'].isin(SIDE)].copy()
missing = set(SIDE) - set(sel['mouse_symbol'])
if missing:
    raise ValueError(f'no 1-to-1 human ortholog for {sorted(missing)} in {IN_ORTHO}')
off_panel = sel[~sel['mouse_symbol'].isin(msd.var_names) | ~sel['human_symbol'].isin(h_genes)]
if len(off_panel):
    raise ValueError(f'not on both MERFISH panels: {off_panel.to_dict("records")}')
sel['side'] = sel['mouse_symbol'].map(SIDE)
sel['in_score'] = sel['mouse_symbol'].isin(SELECTED)
genes = sel.merge(sh[['mouse_symbol', 'human_load', 'mouse_load', 'human_z', 'mouse_z']], on='mouse_symbol', how='left')
genes = genes.set_index('mouse_symbol').loc[list(SIDE)].reset_index()[
    ['side', 'in_score', 'human_symbol', 'mouse_symbol', 'human_load', 'mouse_load', 'human_z', 'mouse_z']]
genes.to_csv(OUT_GENES, sep='\t', index=False)
score_genes = genes[genes['in_score']]
for side in ['A', 'C']:
    g = score_genes[score_genes['side'] == side]
    print(f'score, {side} side ({len(g)}): ' + ', '.join(f'{h}/{m}' for h, m in zip(g['human_symbol'], g['mouse_symbol'])))
print(f'wrote {OUT_GENES}')


def tukey_keep(depth):
    q1, q3 = np.percentile(depth, [25, 75])
    return (depth >= q1 - TUKEY_K * (q3 - q1)) & (depth <= q3 + TUKEY_K * (q3 - q1))


curve_fits = {}   # section -> (raw positions, bin points, polynomial), for the fit-check figure


def flatten_top(pos, section):
    """Depth measured from a 2nd-order polynomial fit to the top of L2/3: per CURVE_N_BINS equal-count bins along
    x, the CURVE_PCTILE depth percentile; depth -> depth - top(x) (vertical offset; the curvature is gentle)."""
    b = pd.qcut(pos['x'], CURVE_N_BINS, labels=False)
    pts = pd.DataFrame({'x': pos['x'].groupby(b).median(),
                        'depth': pos['depth'].groupby(b).quantile(CURVE_PCTILE / 100)})
    top = np.poly1d(np.polyfit(pts['x'], pts['depth'], 2))
    curve_fits[section] = (pos, pts, top)
    print(f'{section}: top of L2/3 fit {top.coeffs[0]:.2e} x^2 + {top.coeffs[1]:.2e} x + {top.coeffs[2]:.1f}; '
          f'range of the fit over the section {np.ptp(top(pos["x"])):.0f} um')
    return pos.assign(depth=pos['depth'] - top(pos['x']))


def add_rel(pos):
    """Relative L2/3 depth (0 = REL_PCTILE[0] pct, 1 = REL_PCTILE[1] pct of the cells' depth) and position along
    the section in the same unit (L2/3 thickness), so mouse and human share one relative scale."""
    lo, hi = np.percentile(pos['depth'], REL_PCTILE)
    return pos.assign(rel_depth=(pos['depth'] - lo) / (hi - lo), rel_x=(pos['x'] - pos['x'].min()) / (hi - lo),
                      l23_thickness_um=hi - lo)


def scores(expr, side_genes):
    """Per gene x / SCORE_MAX_PCTILE percentile, capped at 1, mean over genes; returns z-scored A, C and C - A."""
    out = {}
    for side in ['A', 'C']:
        x = expr[side_genes[side]].to_numpy()
        hi = np.percentile(x, SCORE_MAX_PCTILE, axis=0)
        hi = np.where(hi > 0, hi, 1.0)
        s = np.clip(x / hi, 0.0, 1.0).mean(axis=1)
        out[f'z_{side}'] = (s - s.mean()) / s.std()
    out['ac'] = out['z_C'] - out['z_A']
    return pd.DataFrame(out, index=expr.index)


# --- mouse: L2/3 NR cells and their depth from pia ---
mobs = ad.read_h5ad(IN_MERFISH, backed='r').obs
depth = mobs['depth'].to_numpy()
keep_m = ((mobs['modality'] == 'merfish') & (mobs[MOUSE_SUBCLASS_COL] == MOUSE_SUBCLASS)
          & mobs['sample'].astype(str).str.startswith(MOUSE_CONDITION)).to_numpy() \
    & (depth >= MOUSE_DEPTH_RANGE[0]) & (depth <= MOUSE_DEPTH_RANGE[1])
mc_all = pd.DataFrame({'subclass': MOUSE_SUBCLASS, 'sample': mobs['sample'].astype(str).to_numpy()[keep_m],
                       'depth': depth[keep_m]}, index=pd.Index(mobs.index[keep_m], name='cell'))
print(f'mouse {MOUSE_CONDITION} {MOUSE_SUBCLASS} MERFISH cells in {MOUSE_DEPTH_RANGE} um: {len(mc_all)}')


def load_mouse(section):
    """L2/3 cells of one NR section: (positions, log2(jnorm+1) of the genes), depth outliers removed."""
    path = MOUSE_SECTION_FILE.format(section)
    sd = ad.read_h5ad(path, backed='r')
    absent = [g for g in genes['mouse_symbol'] if g not in sd.var_names]
    if absent:
        raise ValueError(f'{absent} not on the panel of {path}')
    mc = mc_all[(mc_all['subclass'] == MOUSE_SUBCLASS) & (mc_all['sample'] == section)]
    rows = sd.obs_names.get_indexer(mc.index.str.slice(len(section)))
    if (rows < 0).any():
        raise ValueError(f'{(rows < 0).sum()} {section} L2/3 cells not in {path}')
    order = np.argsort(rows)
    rows, mc = rows[order], mc.iloc[order]
    lay = sd.layers[MOUSE_LAYER][rows]
    lay = lay.toarray() if hasattr(lay, 'toarray') else np.asarray(lay)
    ex = pd.DataFrame(np.log2(lay[:, [sd.var_names.get_loc(g) for g in genes['mouse_symbol']]] + 1),
                      index=mc.index, columns=genes['mouse_symbol'])
    pos = flatten_top(pd.DataFrame({'x': sd.obs['width'].to_numpy()[rows], 'depth': mc['depth'].to_numpy()},
                                   index=mc.index), section)
    keep = tukey_keep(pos['depth'].to_numpy())
    print(f'mouse {section} {MOUSE_SUBCLASS}: {len(pos)} cells; {(~keep).sum()} depth outliers removed')
    return add_rel(pos[keep]), ex[keep]


mouse_sections = {sec: load_mouse(sec) for sec in MOUSE_SECTIONS}
mpos, mexpr = mouse_sections[MOUSE_SECTION]
mouse = mpos.join(scores(mexpr, {s: list(score_genes.loc[score_genes['side'] == s, 'mouse_symbol']) for s in 'AC'}))

# --- human ---
def load_human(section):
    """eL2/3.IT cells of one Fang22 section: (positions with relative depth, log1p median-normalized expression
    of the genes), depth outliers removed; rotated so pia is up (superficial eL2/3.IT above the deep H_DEEP cells)."""
    feat = pd.read_csv(os.path.join(HUMAN_DIR, f'{section}.features.csv'))
    hg = pd.read_csv(os.path.join(HUMAN_DIR, f'{section}.genes.csv'))['name'].tolist()
    absent = [g for g in genes['human_symbol'] if g not in hg]
    if absent:
        raise ValueError(f'{section}: {absent} not on the panel')
    trip = pd.read_csv(os.path.join(HUMAN_DIR, f'{section}.matrix.csv'))
    counts = coo_matrix((trip['val'].values, (trip['col'].values - 1, trip['row'].values - 1)),
                        shape=(len(feat), len(hg))).tocsr()
    total = np.asarray(counts.sum(axis=1)).ravel()
    if (total == 0).any():
        raise ValueError(f'{section}: {(total == 0).sum()} cells with zero total counts')
    norm = counts.multiply((np.median(total) / total)[:, None]).tocsr()
    norm.data = np.log1p(norm.data)
    sup, deep = feat[feat[H_SUBCLASS_COL] == H_SUBCLASS], feat[feat[H_SUBCLASS_COL].isin(H_DEEP)]
    if len(sup) < H_MIN_ANCHOR or len(deep) < H_MIN_ANCHOR:
        raise ValueError(f'{section}: too few anchor cells to orient ({len(sup)} {H_SUBCLASS}, {len(deep)} deep)')
    deg = (np.degrees(np.arctan2(sup[H_SRC_Y].mean() - deep[H_SRC_Y].mean(),
                                 sup[H_SRC_X].mean() - deep[H_SRC_X].mean())) - 90) % 360
    rad = np.deg2rad(deg)
    hx = feat[H_SRC_X] * np.cos(rad) + feat[H_SRC_Y] * np.sin(rad)
    hy = -feat[H_SRC_X] * np.sin(rad) + feat[H_SRC_Y] * np.cos(rad)
    m = (feat[H_SUBCLASS_COL] == H_SUBCLASS).to_numpy()
    ex = pd.DataFrame(norm[m][:, [hg.index(g) for g in genes['human_symbol']]].toarray(),
                      index=feat.loc[m, 'name'], columns=genes['human_symbol'])
    pos = pd.DataFrame({'x': hx[m].to_numpy(), 'depth': -hy[m].to_numpy()}, index=ex.index)
    pos = flatten_top(pos, section)
    keep = tukey_keep(pos['depth'].to_numpy())
    print(f'human {section} {H_SUBCLASS}: {len(pos)} cells (rotated {deg:.1f} deg); {(~keep).sum()} depth outliers removed')
    return add_rel(pos[keep]), ex[keep]



hpos, hexpr = load_human(SECTION)
human = hpos.join(scores(hexpr, {s: list(score_genes.loc[score_genes['side'] == s, 'human_symbol']) for s in 'AC'}))

tabs = {f'mouse {MOUSE_SECTION}': mouse, f'human {SECTION}': human}
pd.concat({k.split()[0]: t for k, t in tabs.items()}, names=['species', 'cell']).to_csv(OUT_CELLS, sep='\t')
print(f'wrote {OUT_CELLS}')
for name, t in tabs.items():
    print(f'{name}: Spearman rho with depth: ' + ', '.join(f'{c} {spearmanr(t[c], t["depth"])[0]:.2f}'
                                                         for c, _, _ in COLS))

# --- figure: rows = species, columns = A, C, A vs C ---
fig, axes = plt.subplots(len(tabs), len(COLS), figsize=(5.0 * len(COLS) + 1.0, 3.2 * len(tabs) + 1.0))
for i, (name, t) in enumerate(tabs.items()):
    x0, d0 = t['x'].min(), t['depth'].min()
    for j, (col, title, cmap) in enumerate(COLS):
        ax = axes[i, j]
        lo, hi = np.percentile(t[col], CLIP_PCTILE)
        if col == 'ac':
            if not lo < 0 < hi:
                raise ValueError(f'{name} {col}: {CLIP_PCTILE} pct range ({lo:.2f}, {hi:.2f}) does not contain 0')
            norm_c = TwoSlopeNorm(vmin=lo, vcenter=0, vmax=hi)
        else:
            norm_c = Normalize(vmin=lo, vmax=hi)
        c = t.loc[t[col].abs().sort_values().index]   # strongest values drawn on top
        sc_ = ax.scatter(c['x'] - x0, c['depth'] - d0, s=DOT, c=c[col], cmap=CMAPS.get(cmap, cmap), norm=norm_c,
                         lw=0, rasterized=True)
        ax.set_aspect('equal')
        ax.invert_yaxis()   # pia on top
        ax.set_title(f'{name.split()[0]}: {title}', fontsize=9, loc='left')
        ax.set_xlabel('position along the section (µm)', fontsize=8)
        if j == 0:
            ax.set_ylabel(f'{name}\ndepth from top L2/3 cell (µm)', fontsize=8)
        ax.spines[['top', 'right']].set_visible(False)
        cb = fig.colorbar(sc_, ax=ax, shrink=0.5 if i == 0 else 0.8, pad=0.02, extend='both')   # mouse panels are thin
        if col == 'ac':
            cb.ax.invert_yaxis()   # A side (C0, negative) on top
a_genes = ', '.join(score_genes.loc[score_genes['side'] == 'A', 'human_symbol'])
c_genes = ', '.join(score_genes.loc[score_genes['side'] == 'C', 'human_symbol'])
fig.suptitle(f'L2/3 A vs C score from a hand-picked gene set (same orthologs in both MERFISH panels); '
             f'L2/3 cells, depth outliers removed\nA side: {a_genes}   |   C side: {c_genes}', fontsize=9)
fig.tight_layout()
fig.savefig(OUT_FIG)
print(f'wrote {OUT_FIG}')


# --- figure 2: INSITU_GENES in situ; rows = genes, columns = mouse | human ---
g = genes.set_index('mouse_symbol', drop=False).loc[INSITU_GENES]
panels = [(f'mouse {sec}', pos, ex, 'mouse_symbol') for sec, (pos, ex) in mouse_sections.items()]
panels.append((f'human {SECTION}', human, hexpr.loc[human.index], 'human_symbol'))
with PdfPages(OUT_FIG_GENES) as pdf:
    for gs in [g]:
        fig_h = 2.6 * len(gs) + 1.2
        # all panels share the human panel's box shape; human keeps equal aspect, mouse is stretched to fit
        ht = panels[-1][1]
        box = (ht['rel_depth'].max() - ht['rel_depth'].min()) / (ht['rel_x'].max() - ht['rel_x'].min())
        fig, axes = plt.subplots(len(gs), len(panels), figsize=(4.6 * len(panels) + 1.0, fig_h), squeeze=False)
        fig.subplots_adjust(left=0.04, right=0.94, top=1 - 0.6 / fig_h, bottom=0.6 / fig_h,
                            hspace=0.5, wspace=0.12)
        for i, (_, r) in enumerate(gs.iterrows()):
            for j, (name, t, ex, sym_col) in enumerate(panels):
                ax = axes[i, j]
                v = ex[r[sym_col]]
                scale = np.percentile(v, GENE_PCTILE)
                if scale <= 0:
                    scale = v.max() if v.max() > 0 else 1.0   # gene detected in < 2% of cells
                x, d = t['rel_x'], t['rel_depth']   # one relative scale (L2/3 thickness) in both species
                o = v.sort_values().index   # zero counts = lowest expression (0); high expressers drawn on top
                sc_ = ax.scatter(x[o], d[o], s=DOT, c=np.clip(v[o] / scale, 0, 1), cmap=GENE_CMAP,
                                 vmin=0, vmax=1, lw=0, rasterized=True)
                if name.startswith('human'):
                    ax.set_aspect('equal')
                else:
                    ax.set_aspect('auto')
                    ax.set_box_aspect(box)   # mouse: same panel size as human, aspect not fixed
                ax.invert_yaxis()   # pia on top
                for yl in (0, 1):   # top and bottom of L2/3 (relative depth 0, 1)
                    ax.axhline(yl, color='k', lw=0.7, ls='--')
                ax.axis('off')
                ax.set_title(f"{r[sym_col]} ({name.split()[1]}): detected {(v > 0).mean():.0%}; "
                             f"{AXIS} z {r[sym_col.split('_')[0] + '_z']:.1f}", fontsize=8, loc='left')
        fig.colorbar(sc_, cax=fig.add_axes([0.955, 0.4, 0.008, 0.2]),
                     label=f'log expression / gene {GENE_PCTILE}th pct (per panel)')
        fig.suptitle(f"{' vs '.join(f'{m} ({s})' for m, s in zip(gs['mouse_symbol'], gs['side']))}, L2/3 cells, "
                     'relative L2/3 scale, pia up; '
                     'zero counts shown as 0 on the same scale', fontsize=9)
        pdf.savefig(fig)
        plt.close(fig)
print(f'wrote {OUT_FIG_GENES}')

# --- figure 3: depth profiles on the relative L2/3 scale; all mouse NR and human sections ---
edges = np.linspace(0, 1, N_REL_BINS + 1)
mids = (edges[:-1] + edges[1:]) / 2
prof_secs = {'mouse': {sec: mouse_sections[sec] if sec in mouse_sections else load_mouse(sec)
                       for sec in PROFILE_MOUSE_SECTIONS},
             'human': {sec: (hpos, hexpr) if sec == SECTION else load_human(sec) for sec in PROFILE_HUMAN_SECTIONS}}
sym_of = {'mouse': 'mouse_symbol', 'human': 'human_symbol'}
gi = genes.set_index('mouse_symbol', drop=False)
h_detect = pd.DataFrame({sec: {gi.loc[g, 'human_symbol']: (ex[gi.loc[g, 'human_symbol']] > 0).mean() for g in INSITU_GENES}
                         for sec, (_, ex) in prof_secs['human'].items()}).T
h_keep = (h_detect > H_MIN_DETECT).all(axis=1)
print(f'human sections, fraction of L2/3 cells detecting {list(h_detect.columns)} (keep if all > {H_MIN_DETECT}):')
print(h_detect.assign(keep=h_keep).round(3).to_string())
if not h_keep.any():
    raise ValueError(f'no human section detects all of {list(h_detect.columns)} in > {H_MIN_DETECT} of L2/3 cells')
prof_secs['human'] = {sec: v for sec, v in prof_secs['human'].items() if h_keep[sec]}


def profile(pos, ex, sym):
    """Mean expression per relative-depth bin (cells within [0, 1]), divided by the peak bin."""
    inb = pos['rel_depth'].between(0, 1)
    b = np.clip(np.searchsorted(edges, pos.loc[inb, 'rel_depth'], side='right') - 1, 0, N_REL_BINS - 1)
    prof = ex.loc[pos.index[inb], sym].groupby(b).mean().reindex(range(N_REL_BINS))
    return prof / prof.max()


n_col = max(len(v) for v in PROFILE_GENES.values())
fig, axes = plt.subplots(len(PROFILE_GENES), n_col, figsize=(2.6 * n_col, 3.4 * len(PROFILE_GENES)), squeeze=False)
gord, gax = [], []
for i, (side, gl) in enumerate(PROFILE_GENES.items()):
    gord += gl
    gax += list(axes[i, :len(gl)])
    for ax in axes[i, len(gl):]:
        ax.set_visible(False)
prof_rows = []
for ax, mg in zip(gax, gord):
    for sp, secs in prof_secs.items():
        sym = gi.loc[mg, sym_of[sp]]
        per = pd.DataFrame({sec: profile(pos, ex, sym) for sec, (pos, ex) in secs.items()})
        for sec in per:
            prof_rows.append(pd.DataFrame({'gene': mg, 'species': sp, 'section': sec, 'bin_mid': mids,
                                           'rel_expr': per[sec].to_numpy()}))
        ax.errorbar(per.mean(axis=1), mids, xerr=per.sem(axis=1), color=SPECIES_COLORS[sp], lw=1.8, marker='o',
                    ms=3.5, capsize=2, elinewidth=1, label=f'{sp} (n = {per.shape[1]} sections)')
    side = gi.loc[mg, 'side']
    ax.set_title(f"{gi.loc[mg, 'human_symbol']} / {mg} ({side} side)", fontsize=9, color='C0' if side == 'A' else 'C2')
    ax.set_xlabel('mean expression / peak bin', fontsize=8)
    ax.set_xlim(0, 1.1)
    ax.set_ylim(1, 0)   # top of L2/3 on top
    ax.spines[['top', 'right']].set_visible(False)
gax[0].legend(fontsize=7, frameon=False, loc='lower left')
for ax in axes[:, 0]:
    ax.set_ylabel('relative L2/3 depth (0 = top, 1 = bottom)', fontsize=8)
pd.concat(prof_rows).to_csv(OUT_PROFILES, sep='\t', index=False)
print(f'wrote {OUT_PROFILES}')
fig.suptitle(f"L2/3 depth profiles on one relative scale: {len(prof_secs['mouse'])} mouse NR vs "
             f"{len(prof_secs['human'])} human sections\n(human: {' and '.join(h_detect.columns)} each detected in > "
             f"{H_MIN_DETECT:.0%} of L2/3 cells)\nper section: bin means divided by the section's peak bin; "
             f'mean ± SEM over sections; {N_REL_BINS} bins', fontsize=9)
fig.tight_layout()
fig.savefig(OUT_FIG_PROFILE)
print(f'wrote {OUT_FIG_PROFILE}')

# --- figure 4: INSITU_GENES in every profile section, to inspect depth alignment ---
rows = [(sp, sec, pos, ex) for sp, secs in prof_secs.items() for sec, (pos, ex) in secs.items()]
fig, axes = plt.subplots(len(rows), len(INSITU_GENES), figsize=(4.2 * len(INSITU_GENES) + 1.6, 1.25 * len(rows) + 0.8),
                         squeeze=False)
fig.subplots_adjust(left=0.24, right=0.86, top=1 - 0.9 / (1.25 * len(rows) + 0.8), bottom=0.02, hspace=0.35, wspace=0.16)
for i, (sp, sec, pos, ex) in enumerate(rows):
    for j, mg in enumerate(INSITU_GENES):
        ax = axes[i, j]
        v = ex[gi.loc[mg, sym_of[sp]]]
        scale = np.percentile(v, GENE_PCTILE)
        if scale <= 0:
            scale = v.max() if v.max() > 0 else 1.0   # gene detected in < 2% of cells
        o = v.sort_values().index   # zero counts = lowest expression (0); high expressers drawn on top
        sc_ = ax.scatter(pos.loc[o, 'rel_x'], pos.loc[o, 'rel_depth'], s=2.5, c=np.clip(v[o] / scale, 0, 1),
                         cmap=GENE_CMAP, vmin=0, vmax=1, lw=0, rasterized=True)
        for yl in (0, 1):
            ax.axhline(yl, color='k', lw=0.6, ls='--')
        ax.set_ylim(ALLSEC_YLIM[1], ALLSEC_YLIM[0])   # top of L2/3 on top
        ax.axis('off')
        if i == 0:
            ax.set_title(f"{gi.loc[mg, 'human_symbol']} / {mg} ({gi.loc[mg, 'side']} side)", fontsize=9)
        ax.text(1.01, 0.5, f'{(v > 0).mean():.0%}', transform=ax.transAxes, ha='left', va='center', fontsize=6)
    axes[i, 0].text(-0.03, 0.5, f'{sp}\n{sec}\n(n = {len(pos)})', transform=axes[i, 0].transAxes,
                    ha='right', va='center', fontsize=7)
fig.colorbar(sc_, cax=fig.add_axes([0.9, 0.45, 0.012, 0.1]), label=f'log expr / {GENE_PCTILE}th pct')
fig.suptitle('Depth alignment check: L2/3 cells of every profile section on the relative scale\n'
             '(dashed = top / bottom of L2/3; % right of each panel = cells with the gene detected)', fontsize=8)
fig.savefig(OUT_FIG_ALLSEC)
print(f'wrote {OUT_FIG_ALLSEC}')

# --- figure 5: every profile section, fitted top of L2/3 (left: raw, right: flattened, relative scale) ---
fig, axes = plt.subplots(len(rows), 2, figsize=(11, 2.6 * len(rows)), squeeze=False)
for i, (sp, sec, pos, _) in enumerate(rows):
    raw, pts, top = curve_fits[sec]
    xs = np.linspace(raw['x'].min(), raw['x'].max(), 200)
    ax = axes[i, 0]
    ax.scatter(raw['x'], raw['depth'], s=2, c='0.6', lw=0, rasterized=True)
    ax.scatter(pts['x'], pts['depth'], s=14, c='C3', zorder=3, label=f'{CURVE_PCTILE}th pct depth per bin')
    ax.plot(xs, top(xs), color='C3', lw=1.2, label='2nd-order fit (top of L2/3)')
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.set_title(f'{sp} {sec}: raw L2/3 cells (µm)', fontsize=8, loc='left')
    ax = axes[i, 1]
    ax.scatter(pos['rel_x'], pos['rel_depth'], s=2, c='0.6', lw=0, rasterized=True)
    for yl in (0, 1):
        ax.axhline(yl, color='k', lw=0.6, ls='--')
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.set_title('flattened, outliers removed, relative scale', fontsize=8, loc='left')
axes[0, 0].legend(fontsize=7, frameon=False, loc='lower left')
fig.tight_layout()
fig.savefig(OUT_FIG_CURVE)
print(f'wrote {OUT_FIG_CURVE}')
