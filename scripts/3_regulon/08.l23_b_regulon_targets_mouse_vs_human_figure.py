"""Figure of 04: B'-enriched AP-1 and EGR regulons, target-gene activity over normal development, mouse vs human.

Left: mouse Yoo25 L2/3 NR, per age; right: human Wang25 V1 L2/3, per developmental stage. Mean +/- SEM
over cells, divided by the peak; one colour per TF in both panels (tab20 colour i for TF i, AP-1 then
EGR), marker = member order within its family.

Reads:
  results/3_regulon/04.l23_b_regulon_targets_nr_mouse.tsv
  results/3_regulon/04.wang25_l23_ap1_egr_regulon_targets_over_stage.tsv
Outputs:
  figures/3_regulon/08.l23_b_regulon_targets_mouse_vs_human.pdf
"""

import os

import pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES_DIR   = os.path.join(REPO_ROOT, 'results', '3_regulon')
IN_MOUSE  = os.path.join(RES_DIR, '04.l23_b_regulon_targets_nr_mouse.tsv')
IN_HUMAN  = os.path.join(RES_DIR, '04.wang25_l23_ap1_egr_regulon_targets_over_stage.tsv')
OUT_PDF   = os.path.join(REPO_ROOT, 'figures', '3_regulon', '08.l23_b_regulon_targets_mouse_vs_human.pdf')

SIGN = '+/+'
HUMAN_STAGES = {'Second_trimester': '2nd trim.', 'Third_trimester': '3rd trim.',
                'Infancy': 'Infancy', 'Adolescence': 'Adolescence'}
CMAP = 'tab20'
LINE_MARKERS = ['o', 's', '^', 'D']

os.makedirs(os.path.dirname(OUT_PDF), exist_ok=True)


def draw(ax, s, color, marker, label):
    ax.errorbar(s['x'], s['mean_norm'], yerr=s['sem_norm'], color=color, marker=marker,
                markersize=4.5, capsize=2.0, linewidth=1.3, elinewidth=0.9, label=label)


mouse = pd.read_csv(IN_MOUSE, sep='\t')
hum = pd.read_csv(IN_HUMAN, sep='\t')
if set(hum['stage']) != set(HUMAN_STAGES):
    raise ValueError(f'unexpected human stages: {sorted(set(hum["stage"]))}')
tfs = list(dict.fromkeys(mouse['TF']))
family = mouse.drop_duplicates('TF').set_index('TF')['family']
style = {}
for fam in dict.fromkeys(family):
    for j, tf in enumerate(family.index[family == fam]):
        style[tf] = (plt.get_cmap(CMAP)(len(style)), LINE_MARKERS[j])

days = sorted(mouse['day'].unique())
stage_x = {s: i for i, s in enumerate(HUMAN_STAGES)}
fig, (axm, axh) = plt.subplots(1, 2, figsize=(6.0, 3.6), sharey=True)
for tf in tfs:
    color, marker = style[tf]
    m = mouse[mouse['TF'] == tf].assign(x=lambda d: d['day'])
    draw(axm, m, color, marker, None)
    h = hum[hum['TF'] == tf.upper()].assign(x=lambda d: d['stage'].map(stage_x)).sort_values('x')
    draw(axh, h, color, marker, f'{tf} ({m["n_targets_used"].iloc[0]} / {h["n_targets_used"].iloc[0]})')
axm.set_xticks(days)
axm.set_xticklabels([f'P{d}' for d in days], fontsize=7)
axm.set_xlim(days[0] - 1, days[-1] + 1)
axm.set_xlabel('postnatal age')
axm.set_title('mouse V1 L2/3 (Yoo25, NR)', fontsize=9)
axm.set_ylabel('regulon target activity\n(per-gene max-scaled mean, /peak)', fontsize=8)
axh.set_xticks(range(len(HUMAN_STAGES)))
axh.set_xticklabels(list(HUMAN_STAGES.values()), rotation=45, ha='right', fontsize=7)
axh.set_xlim(-0.4, len(HUMAN_STAGES) - 0.6)
axh.set_xlabel('developmental stage')
axh.set_title('human V1 L2/3 (Wang25)', fontsize=9)
for ax in (axm, axh):
    ax.set_ylim(bottom=0)
    ax.axhline(1.0, color='0.8', lw=0.6, zorder=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
fig.suptitle(f'AP-1 & EGR {SIGN} regulon target activity over development, mouse vs human\n'
             'mean ± SEM over cells', fontsize=10)
axh.legend(fontsize=7, frameon=False, title='TF (n targets:\nmouse / human)', title_fontsize=7,
           loc='center left', bbox_to_anchor=(1.02, 0.5))
fig.tight_layout()
fig.savefig(OUT_PDF, bbox_inches='tight')
plt.close(fig)
print(f'wrote {OUT_PDF}')
