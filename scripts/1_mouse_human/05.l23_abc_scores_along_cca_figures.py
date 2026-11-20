"""Figures of 03: mouse L2/3 archetype scores A, B, C along the conserved CCA axes, mouse vs human, and
on the CCA1 x CCA2 embedding.

05.l23_abc_normalized_compact.pdf: one panel per archetype (A along CCA1, B along CCA2, C along CCA1):
x = CCA percentile within species; y = score z-scored within species, mean +/- SEM per bin; mouse and human on twin y axes (mouse
left, darker shade of the archetype color; human right, lighter shade), 0 centered on both. The human
legend entry carries 03's human-mouse r and permutation p.

05.l23_abc_scores_cca_embedding.pdf: ported from scomics it_evo/18b; the mouse archetype triangle is drawn on
the mouse row, the human archetypes (18b's human PCHA vertices) are not shown.

Cells at their CCA1, CCA2 coordinates (script 03, from its canonical weights): rows = species
(human Jorstad23 cells with the mouse-derived scores of script 02; mouse Cheng22 cells with the native
scores of script 02), columns = score B and score C - A. CCA1 is reflected in both rows (CCA1') so the
two embeddings share one orientation; the canonical r of each axis is in the axis label. Color range per
panel: B = its 5th-95th percentile; C - A = +/- the 95th percentile of |C - A| (gray = 0).
Mouse row: the mouse PCHA archetypes (script 02) projected like the mouse cells (script 03), as a triangle.

Reads:
  results/1_mouse_human/03.l23_abc_scores_along_cca_profiles.tsv
  results/1_mouse_human/03.l23_abc_scores_along_cca_stats.tsv
  results/1_mouse_human/03.l23_cca_cells.tsv
  results/1_mouse_human/03.l23_cca_mouse_archetypes.tsv
  results/1_mouse_human/03.l23_axis_cca_weights_human.tsv
Outputs:
  figures/1_mouse_human/05.l23_abc_normalized_compact.pdf
  figures/1_mouse_human/05.l23_abc_scores_cca_embedding.pdf
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgb

from l23.constants import ARCH_COLORS

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

RES_DIR     = os.path.join(REPO_ROOT, 'results', '1_mouse_human')
IN_PROFILES = os.path.join(RES_DIR, '03.l23_abc_scores_along_cca_profiles.tsv')
IN_STATS    = os.path.join(RES_DIR, '03.l23_abc_scores_along_cca_stats.tsv')
IN_CELLS    = os.path.join(RES_DIR, '03.l23_cca_cells.tsv')
IN_M_AA     = os.path.join(RES_DIR, '03.l23_cca_mouse_archetypes.tsv')
IN_WEIGHTS  = os.path.join(RES_DIR, '03.l23_axis_cca_weights_human.tsv')
OUT_PDF     = os.path.join(REPO_ROOT, 'figures', '1_mouse_human', '05.l23_abc_normalized_compact.pdf')
OUT_EMBED   = os.path.join(REPO_ROOT, 'figures', '1_mouse_human', '05.l23_abc_scores_cca_embedding.pdf')

PANELS  = [('CCA1', 'A'), ('CCA2', 'B'), ('CCA1', 'C')]   # (axis, archetype), left to right
MARKERS = {'mouse': 's', 'human': 'o'}
SHADE   = {'mouse': ('black', 0.25), 'human': ('white', 0.45)}   # blend of the archetype color

ROWS         = [('human', 'Human (Jorstad23), mouse scores'), ('mouse', 'Mouse (Cheng22)')]
SCORE_PCTILE = (5, 95)   # color range of score B
DIFF_PCTILE  = 95        # color range of C - A: +/- this percentile of |C - A|
POINT_S      = 3
FS_TITLE     = 10
FS_LABEL     = 9
FS_TICK      = 8
CMAP_B    = LinearSegmentedColormap.from_list('b', ['gainsboro', 'lightgray', ARCH_COLORS['B']])
CMAP_DIFF = LinearSegmentedColormap.from_list('ca', [ARCH_COLORS['A'], 'lightgray', ARCH_COLORS['C']])

os.makedirs(os.path.dirname(OUT_PDF), exist_ok=True)
plt.rcParams['pdf.fonttype'] = 42


def species_shade(color, name):
    """Archetype color blended toward SHADE[name] (dark for mouse, light for human)."""
    target, frac = SHADE[name]
    return tuple((1 - frac) * np.array(to_rgb(color)) + frac * np.array(to_rgb(target)))


def embedding_panels(c):
    """(title, values, cmap, (vmin, vmax), colorbar label) for one species."""
    b = c["score_B"].to_numpy()
    ca = (c["score_C"] - c["score_A"]).to_numpy()
    lim = np.percentile(np.abs(ca), DIFF_PCTILE)
    return [("score B", b, CMAP_B, tuple(np.percentile(b, SCORE_PCTILE)), 'archetype score'),
            ("score C − A", ca, CMAP_DIFF, (-lim, lim), 'score difference')]


def draw_triangle(ax, aa):
    """Closed triangle through the archetypes (CCA1 reflected like the cells), labeled."""
    x, y = -aa['CCA1'].to_numpy(), aa['CCA2'].to_numpy()
    ax.plot(np.r_[x, x[0]], np.r_[y, y[0]], '-', color='black', lw=1)
    ax.scatter(x, y, marker='D', color='black', s=30, zorder=3)
    for xi, yi, label in zip(x, y, aa['archetype']):
        ax.annotate(label, (xi, yi), textcoords='offset points', xytext=(5, 5),
                    fontsize=FS_LABEL, fontweight='bold', zorder=4)


def draw_embedding(fig, ax, c, vals, cmap, vlim, cbar_label, r):
    sc = ax.scatter(-c['CCA1'], c['CCA2'], c=vals, cmap=cmap, vmin=vlim[0], vmax=vlim[1],
                    s=POINT_S, linewidths=0, rasterized=True)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel(f"CCA1' (r = {r['CCA1']:.3f})", fontsize=FS_LABEL)
    ax.set_ylabel(f"CCA2 (r = {r['CCA2']:.3f})", fontsize=FS_LABEL)
    ax.tick_params(labelsize=FS_TICK)
    ax.spines[['top', 'right']].set_visible(False)
    cb = fig.colorbar(sc, ax=ax, shrink=0.8)
    cb.set_label(cbar_label, fontsize=FS_LABEL)
    cb.ax.tick_params(labelsize=FS_TICK)


# --- scores along the CCA axes ---
prof = pd.read_csv(IN_PROFILES, sep='\t').set_index(['axis', 'archetype', 'species'])
st = pd.read_csv(IN_STATS, sep='\t').set_index(['axis', 'archetype'])

fig, axes = plt.subplots(1, len(PANELS), figsize=(3.4 * len(PANELS), 2.8), squeeze=False)
for ax, (axis, disp) in zip(axes[0], PANELS):
    s = st.loc[(axis, disp)]
    n_perm = int(s['n_perm'])
    p_txt = f'p < {1 / n_perm:g}' if np.isclose(s['p_perm'], 1 / (1 + n_perm)) else f'p = {s["p_perm"]:.3g}'
    stat_txt = f'r = {s["r"]:.2f}, {p_txt}'
    twin = {'mouse': ax, 'human': ax.twinx()}
    handles = []
    for name, a in twin.items():
        t = prof.loc[(axis, disp, name)]
        c = species_shade(ARCH_COLORS[disp], name)
        h, = a.plot(t['pct'], t['mean'], color=c, marker=MARKERS[name], ms=2.5, lw=1.3,
                    label=f'{disp} {name}' + (f' ({stat_txt})' if name == 'human' else ''))
        a.fill_between(t['pct'], t['mean'] - t['sem'], t['mean'] + t['sem'], color=c, alpha=0.25, lw=0)
        handles.append(h)
    for a in twin.values():                               # center 0 on both y axes
        lim = np.abs(a.get_ylim()).max()
        a.set_ylim(-lim, lim)
        a.tick_params(labelsize=7)
        a.spines['top'].set_visible(False)
    ax.axhline(0, color='lightgray', lw=0.8, zorder=0)
    ax.set_ylabel('mouse score (z, dark)', fontsize=8)
    twin['human'].set_ylabel('human score (z, light)', fontsize=8)
    ax.set_xlabel(f'{axis} percentile within species', fontsize=8)
    ax.legend(handles=handles, frameon=False, fontsize=7, loc='lower center', bbox_to_anchor=(0.5, 1.0))
fig.suptitle('L2/3 mouse vs human along conserved CCA axes', fontsize=9)
fig.tight_layout()
fig.savefig(OUT_PDF, bbox_inches='tight', dpi=300)
plt.close(fig)
print(f'Saved {OUT_PDF}')

# --- cells on the CCA1 x CCA2 embedding ---
cells = pd.read_csv(IN_CELLS, sep='\t')
missing = {s for s, _ in ROWS} - set(cells['species'])
if missing:
    raise ValueError(f'species {sorted(missing)} absent from {IN_CELLS}')
m_aa = pd.read_csv(IN_M_AA, sep='\t')
r = pd.read_csv(IN_WEIGHTS, sep='\t', index_col=0)['canonical_r']

fig, axs = plt.subplots(len(ROWS), 2, figsize=(8.4, 7), layout='constrained')
for i, (species, row_title) in enumerate(ROWS):
    c = cells[cells['species'] == species]
    for j, (title, vals, cmap, vlim, cbar_label) in enumerate(embedding_panels(c)):
        draw_embedding(fig, axs[i, j], c, vals, cmap, vlim, cbar_label, r)
        axs[i, j].set_title(f'{row_title}\n{title}', fontsize=FS_TITLE)
        if species == 'mouse':
            draw_triangle(axs[i, j], m_aa)
    print(f'{species}: {len(c)} cells')
fig.savefig(OUT_EMBED, dpi=300)
print(f'Saved {OUT_EMBED}')
