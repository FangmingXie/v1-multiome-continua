"""Figure of 05: mouse L2/3 archetype scores A', B', C' along the conserved CCA axes, mouse vs human.

One panel per archetype (A' along CCA1, B' along CCA2, C' along CCA1): x = CCA percentile within
species; y = score z-scored within species, mean +/- SEM per bin; mouse and human on twin y axes (mouse
left, darker shade of the archetype color; human right, lighter shade), 0 centered on both. The human
legend entry carries 05's human-mouse r and permutation p.

Reads:
  results/1_mouse_human/05.l23_abc_scores_along_cca_profiles.tsv
  results/1_mouse_human/05.l23_abc_scores_along_cca_stats.tsv
Outputs:
  figures/1_mouse_human/07.l23_abc_normalized_compact.pdf
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb

from l23.constants import ARCH_COLORS

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

RES_DIR     = os.path.join(REPO_ROOT, 'results', '1_mouse_human')
IN_PROFILES = os.path.join(RES_DIR, '05.l23_abc_scores_along_cca_profiles.tsv')
IN_STATS    = os.path.join(RES_DIR, '05.l23_abc_scores_along_cca_stats.tsv')
OUT_PDF     = os.path.join(REPO_ROOT, 'figures', '1_mouse_human', '07.l23_abc_normalized_compact.pdf')

PANELS  = [('CCA1', "A'"), ('CCA2', "B'"), ('CCA1', "C'")]   # (axis, archetype), left to right
MARKERS = {'mouse': 's', 'human': 'o'}
SHADE   = {'mouse': ('black', 0.25), 'human': ('white', 0.45)}   # blend of the archetype color

os.makedirs(os.path.dirname(OUT_PDF), exist_ok=True)
plt.rcParams['pdf.fonttype'] = 42


def species_shade(color, name):
    """Archetype color blended toward SHADE[name] (dark for mouse, light for human)."""
    target, frac = SHADE[name]
    return tuple((1 - frac) * np.array(to_rgb(color)) + frac * np.array(to_rgb(target)))


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
