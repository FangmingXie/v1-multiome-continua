"""Mouse L2/3 archetype scores (A', B', C') along the conserved CCA1 and CCA2 axes, human vs mouse.

A cell's CCA coordinate is `(VX coords - mean) . w` with the Gate-A canonical weights of script 05;
human cells carry the mouse-derived scores (script 04), mouse cells the native scores (script 03).
Letters: 03's A/B/C are relabelled to the published primed letters (RELABEL) and keep their colors
(A'->C0, B'->C1, C'->C2). CCA axes are in their raw orientation.

Profiles normalized for the cross-species comparison, human and mouse in separate rows:
  - x = CCA percentile within species (bin k is the same slice of each species' distribution);
  - y = archetype score z-scored over all cells of that species (removes the level / spread gap of
    human scores imputed from mouse), mean +/- SEM in N_BINS equal-count bins;
  - per line, Pearson r between the N_BINS paired human and mouse bin means, with a one-sided
    permutation p (scores shuffled across cells within each species, N_PERM times).
Two compact figures show the same profiles and statistics in one row, mouse and human on twin y axes
(mouse left, darker shade of the archetype color; human right, lighter shade): PAIRS (A' and C' along
CCA1, B' along CCA2) and ABC_PANELS (one archetype per panel).

Reads:
  results/1_mouse_human/01.mouse_l23_varimax_coords.tsv
  results/1_mouse_human/02.human_l23_varimax_coords.tsv
  results/1_mouse_human/03.mouse_l23_archetype_scores.tsv
  results/1_mouse_human/04.mouse_l23_archetype_scores_on_human_cells.tsv
  results/1_mouse_human/05.l23_axis_cca_weights_{human,mouse}.tsv
Outputs:
  figures/1_mouse_human/06.l23_ac_cca1_b_cca2_normalized.pdf
  figures/1_mouse_human/06.l23_ac_cca1_b_cca2_normalized_compact.pdf
  figures/1_mouse_human/06.l23_abc_normalized_compact.pdf
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import to_rgb

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

RES_DIR      = os.path.join(REPO_ROOT, 'results', '1_mouse_human')
FIG_DIR      = os.path.join(REPO_ROOT, 'figures', '1_mouse_human')
IN_M_COORDS  = os.path.join(RES_DIR, '01.mouse_l23_varimax_coords.tsv')
IN_H_COORDS  = os.path.join(RES_DIR, '02.human_l23_varimax_coords.tsv')
IN_M_SCORES  = os.path.join(RES_DIR, '03.mouse_l23_archetype_scores.tsv')
IN_H_SCORES  = os.path.join(RES_DIR, '04.mouse_l23_archetype_scores_on_human_cells.tsv')
IN_H_WEIGHTS = os.path.join(RES_DIR, '05.l23_axis_cca_weights_human.tsv')
IN_M_WEIGHTS = os.path.join(RES_DIR, '05.l23_axis_cca_weights_mouse.tsv')
OUT_PDF_NORM        = os.path.join(FIG_DIR, '06.l23_ac_cca1_b_cca2_normalized.pdf')
OUT_PDF_COMPACT     = os.path.join(FIG_DIR, '06.l23_ac_cca1_b_cca2_normalized_compact.pdf')
OUT_PDF_COMPACT_ABC = os.path.join(FIG_DIR, '06.l23_abc_normalized_compact.pdf')

HUMAN_VX   = ['VX2', 'VX6', 'VX7', 'VX8', 'VX9', 'VX10']   # Gate A, as script 05
MOUSE_VX   = ['VX1', 'VX2', 'VX5', 'VX7', 'VX8', 'VX9']
RELABEL    = {'C': "A'", 'B': "B'", 'A': "C'"}             # 03's letters -> published letters
CCA_AXES   = ['CCA1', 'CCA2']
ARCH_COLOR = {"A'": 'C0', "B'": 'C1', "C'": 'C2'}
SPECIES_STYLE = {'Human Jorstad23': dict(ls='-', marker='o'),
                 'Mouse Cheng22':   dict(ls='--', marker='s')}
N_BINS     = 20                                            # equal-count bins along each CCA axis
PAIRS      = [('CCA1', ["A'", "C'"]), ('CCA2', ["B'"])]     # panel: axis -> archetypes
ABC_PANELS = [('CCA1', ["A'"]), ('CCA2', ["B'"]), ('CCA1', ["C'"])]
N_PERM     = 1000                                          # permutations for the human-mouse r null
SEED       = 0
SHADE      = {'Mouse Cheng22': ('black', 0.25), 'Human Jorstad23': ('white', 0.45)}  # compact: blend

os.makedirs(FIG_DIR, exist_ok=True)


def cca_coords(coords_path, vx_cols, weights_path):
    """(cell index, cells x 2 CCA coordinates, [CCA1 r, CCA2 r]) for one species."""
    vx_df   = pd.read_csv(coords_path, sep='\t', index_col=0)
    weights = pd.read_csv(weights_path, sep='\t', index_col=0)
    W = weights.loc[CCA_AXES, vx_cols].values.T
    C = vx_df[vx_cols].values
    return vx_df.index, (C - C.mean(axis=0)) @ W, weights.loc[CCA_AXES, 'canonical_r'].values


def species_shade(color, name):
    """Archetype color blended toward SHADE[name] (dark for mouse, light for human)."""
    target, frac = SHADE[name]
    return tuple((1 - frac) * np.array(to_rgb(color)) + frac * np.array(to_rgb(target)))


m_scores = pd.read_csv(IN_M_SCORES, sep='\t', index_col=0)
h_scores = pd.read_csv(IN_H_SCORES, sep='\t', index_col=0)
if list(h_scores.columns) != list(m_scores.columns):
    raise ValueError(f'Score columns differ between human ({list(h_scores.columns)}) '
                     f'and mouse ({list(m_scores.columns)}).')
letters = sorted((c.split('_')[1] for c in m_scores.columns), key=lambda L: RELABEL[L])

h_index, h_cca, r = cca_coords(IN_H_COORDS, HUMAN_VX, IN_H_WEIGHTS)
m_index, m_cca, _ = cca_coords(IN_M_COORDS, MOUSE_VX, IN_M_WEIGHTS)
if not h_scores.index.equals(h_index):
    raise ValueError(f'Human scores ({IN_H_SCORES}) and coords ({IN_H_COORDS}) cover different cells/order.')
if not m_scores.index.equals(m_index):
    raise ValueError(f'Mouse scores ({IN_M_SCORES}) and coords ({IN_M_COORDS}) cover different cells/order.')

species = [('Human Jorstad23', h_cca, h_scores), ('Mouse Cheng22', m_cca, m_scores)]
print(f'--- L2/3 {[RELABEL[L] for L in letters]}; human cells {len(h_cca)}, mouse cells {len(m_cca)}; '
      f'CCA r = {r[0]:.3f}, {r[1]:.3f} ---')
old_letter = {RELABEL[L]: L for L in letters}
plt.rcParams['pdf.fonttype'] = 42

# --- normalized figure: x = CCA percentile, y = within-species z-score, human-mouse r ---
rng = np.random.default_rng(SEED)
pct = (np.arange(N_BINS) + 0.5) * 100 / N_BINS                # bin centers on the percentile scale


def bin_means(bins, y):
    return np.bincount(bins, weights=y, minlength=N_BINS) / np.bincount(bins, minlength=N_BINS)


# rows = species (each with its own y scale), columns = PAIRS
profiles, stat_txt = {}, {}                                  # reused by the compact figures
fig, axes = plt.subplots(len(species), len(PAIRS), figsize=(4.0 * len(PAIRS), 2.8 * len(species)),
                         sharex=True, squeeze=False)
for j, (axis, disps) in enumerate(PAIRS):
    i = CCA_AXES.index(axis)
    stats = []
    for disp in disps:
        L = old_letter[disp]
        prof, z_bins = {}, {}
        for row, (name, cca, scores) in enumerate(species):
            y = scores[f'score_{L}'].values
            z = (y - y.mean()) / y.std()
            bins = pd.qcut(cca[:, i], N_BINS, labels=False)
            g = pd.Series(z).groupby(bins)
            prof[name], z_bins[name] = g.mean().values, (bins, z)
            profiles[(axis, disp, name)] = (prof[name], g.sem().values)
            ax = axes[row][j]
            ax.plot(pct, prof[name], color=ARCH_COLOR[disp], ms=3, lw=1.5,
                    label=f'Score {disp}', **SPECIES_STYLE[name])
            ax.fill_between(pct, prof[name] - g.sem().values, prof[name] + g.sem().values,
                            color=ARCH_COLOR[disp], alpha=0.25, lw=0)
        (hb, hz), (mb, mz) = z_bins['Human Jorstad23'], z_bins['Mouse Cheng22']
        r_obs = np.corrcoef(prof['Human Jorstad23'], prof['Mouse Cheng22'])[0, 1]
        r_null = np.array([np.corrcoef(bin_means(hb, rng.permutation(hz)),
                                       bin_means(mb, rng.permutation(mz)))[0, 1]
                           for _ in range(N_PERM)])
        p_perm = (1 + np.sum(r_null >= r_obs)) / (1 + N_PERM)
        print(f'  {disp} along {axis}: human-mouse r = {r_obs:.3f}, permutation p = {p_perm:.4f}')
        p_txt = f'p < {1 / N_PERM:g}' if p_perm == 1 / (1 + N_PERM) else f'p = {p_perm:.3g}'
        stat_txt[(axis, disp)] = f'r = {r_obs:.2f}, {p_txt}'
        stats.append(f'{disp}: {stat_txt[(axis, disp)]}')
    for row, (name, _, _) in enumerate(species):
        ax = axes[row][j]
        ax.axhline(0, color='lightgray', lw=0.8, zorder=0)
        ax.set_ylabel(f'{name.split()[0]}: score (z)')
        ax.legend(frameon=False, fontsize=7)
        sns.despine(ax=ax)
    axes[0][j].set_title(' & '.join(f'Score {d}' for d in disps) + f' along {axis}\n'
                         + 'human vs mouse ' + '; '.join(stats), fontsize=9)
    axes[-1][j].set_xlabel(f'{axis} percentile within species (r={r[i]:.3f})')
fig.suptitle(f'L2/3: human (top) vs mouse (bottom), z within species; mean +/- SEM in {N_BINS} '
             f'equal-count bins')
fig.tight_layout()
fig.savefig(OUT_PDF_NORM, bbox_inches='tight', dpi=300)
plt.close(fig)
print(f'  Saved {OUT_PDF_NORM}')

# --- compact figures: one row; mouse (left y, dark shade) and human (right y, light shade) ---
def compact_figure(panels, out_pdf):
    """One panel per (CCA axis, [archetypes]); mouse and human on twin y axes."""
    fig, axes = plt.subplots(1, len(panels), figsize=(3.4 * len(panels), 2.8), squeeze=False)
    for ax, (axis, disps) in zip(axes[0], panels):
        twin = {'Mouse Cheng22': ax, 'Human Jorstad23': ax.twinx()}
        handles = []
        for disp in disps:
            for name, a in twin.items():
                m, se = profiles[(axis, disp, name)]
                c = species_shade(ARCH_COLOR[disp], name)
                h, = a.plot(pct, m, color=c, marker=SPECIES_STYLE[name]['marker'], ms=2.5, lw=1.3,
                            label=f"{disp} {name.split()[0].lower()}"
                                  + (f' ({stat_txt[(axis, disp)]})' if name.startswith('Human') else ''))
                a.fill_between(pct, m - se, m + se, color=c, alpha=0.25, lw=0)
                handles.append(h)
        for a in twin.values():                               # center 0 on both y axes
            lim = np.abs(a.get_ylim()).max()
            a.set_ylim(-lim, lim)
            a.tick_params(labelsize=7)
            a.spines['top'].set_visible(False)
        ax.axhline(0, color='lightgray', lw=0.8, zorder=0)
        ax.set_ylabel('mouse score (z, dark)', fontsize=8)
        twin['Human Jorstad23'].set_ylabel('human score (z, light)', fontsize=8)
        ax.set_xlabel(f'{axis} percentile within species', fontsize=8)
        ax.legend(handles=handles, frameon=False, fontsize=7, loc='lower center',
                  bbox_to_anchor=(0.5, 1.0))
    fig.suptitle('L2/3 mouse vs human along conserved CCA axes', fontsize=9)
    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f'  Saved {out_pdf}')


compact_figure(PAIRS, OUT_PDF_COMPACT)
compact_figure(ABC_PANELS, OUT_PDF_COMPACT_ABC)
print('\nDone.')
