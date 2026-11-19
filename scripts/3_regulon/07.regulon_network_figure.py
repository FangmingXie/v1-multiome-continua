"""Figure of 04: TF -> target-gene network of the L2/3-archetype-enriched Yoo25 regulons (activating, +/+).

TFs are squares coloured by their enriched archetype; marker genes are circles coloured by the archetype
they mark; a target that is a TF but not enriched is a small square (always labelled). All edges are
drawn alike (thin, arrowed, in the TF's colour); self-regulation is a small loop on the square. Labels:
TF names only. PDF: the network (edges, nodes) is one rasterized image; all text stays vector.

Reads:
  results/3_regulon/04.regulon_network_nodes.tsv
  results/3_regulon/04.regulon_network_edges.tsv
Outputs:
  figures/3_regulon/07.regulon_network_free.pdf
"""

import os

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import offset_copy

from l23.constants import ARCH_COLORS

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_NODES = os.path.join(REPO_ROOT, 'results', '3_regulon', '04.regulon_network_nodes.tsv')
INPUT_EDGES = os.path.join(REPO_ROOT, 'results', '3_regulon', '04.regulon_network_edges.tsv')
OUT_PDF = os.path.join(REPO_ROOT, 'figures', '3_regulon', '07.regulon_network_free.pdf')

SIGN = '+/+'
TF_NODE_SIZE = 750
GENE_NODE_SIZE = 160
EDGE_WIDTH = 2.4
EDGE_ALPHA = 0.35
LOOP_SPAN_PT = 8      # self-loop: arc between two points this far apart on the square's right edge
LOOP_RAD = -1.6       # arc3 bulge; |rad| * span / 2 ~ how far the loop sticks out (points)
# PDF: edges and nodes (networkx zorder 1-2, loops 1.5) are rasterized as ONE image below
# RASTER_ZORDER; labels are lifted to LABEL_ZORDER so text stays vector (legend/title already are)
RASTER_ZORDER = 3
LABEL_ZORDER = 4
PDF_RASTER_DPI = 300

os.makedirs(os.path.dirname(OUT_PDF), exist_ok=True)


def load_graph():
    """The 04 graph, nodes and edges in their original order (draw order follows it)."""
    nodes = pd.read_csv(INPUT_NODES, sep='\t')
    edges = pd.read_csv(INPUT_EDGES, sep='\t')
    g = nx.DiGraph()
    for r in nodes.itertuples():
        g.add_node(r.node, kind=r.kind, arch=r.arch, is_tf=r.is_tf)
    g.add_edges_from(zip(edges['TF'], edges['Gene']))
    if set(g) != set(nodes['node']):
        raise ValueError(f'edges of {INPUT_EDGES} reach nodes missing from {INPUT_NODES}')
    return g, {r.node: np.array([r.x, r.y]) for r in nodes.itertuples()}


def plot(g, pos):
    tf_nodes = [n for n, d in g.nodes(data=True) if d['kind'] == 'tf']
    gene_nodes = [n for n, d in g.nodes(data=True) if d['kind'] == 'gene']
    tgt_tf_nodes = [n for n in gene_nodes if g.nodes[n]['is_tf']]
    dot_nodes = [n for n in gene_nodes if not g.nodes[n]['is_tf']]
    color = {n: ARCH_COLORS[d['arch']] for n, d in g.nodes(data=True)}

    n_tf_edges = sum(g.nodes[v]['kind'] == 'tf' for _, v in g.edges())
    node_size = {n: TF_NODE_SIZE if g.nodes[n]['kind'] == 'tf' else GENE_NODE_SIZE for n in g}

    fig, ax = plt.subplots(figsize=(16, 16))
    # all edges alike: thin, arrowed; slight curve keeps reciprocal TF <-> TF pairs apart.
    # nodelist/node_size make the arrowheads stop at each target's border
    nodelist = list(g)
    loops = list(nx.selfloop_edges(g))
    edgelist = [e for e in g.edges() if e[0] != e[1]]
    nx.draw_networkx_edges(g, pos, ax=ax, edgelist=edgelist,
                           edge_color=[color[u] for u, _ in edgelist],
                           alpha=EDGE_ALPHA, width=EDGE_WIDTH, arrows=True, arrowstyle='-|>',
                           arrowsize=12, connectionstyle='arc3,rad=0.1', nodelist=nodelist,
                           node_size=[node_size[n] for n in nodelist])
    # networkx sizes self-loops at 10% of the plot height; draw them instead as a small arc on
    # the square's right side, fixed in points so the size is independent of the layout's scale
    half = np.sqrt(TF_NODE_SIZE) / 2
    for n, _ in loops:
        ax.annotate('', xy=pos[n], xytext=pos[n],
                    xycoords=offset_copy(ax.transData, fig=fig, x=half, y=-LOOP_SPAN_PT / 2,
                                         units='points'),
                    textcoords=offset_copy(ax.transData, fig=fig, x=half, y=LOOP_SPAN_PT / 2,
                                           units='points'),
                    zorder=1.5,
                    arrowprops=dict(arrowstyle='-|>', mutation_scale=6, color=color[n],
                                    alpha=EDGE_ALPHA, lw=EDGE_WIDTH, shrinkA=0, shrinkB=0,
                                    connectionstyle=f'arc3,rad={LOOP_RAD}'))
    nx.draw_networkx_nodes(g, pos, nodelist=dot_nodes, ax=ax, node_shape='o', node_size=GENE_NODE_SIZE,
                           node_color=[color[n] for n in dot_nodes], edgecolors='white',
                           linewidths=0.4)
    nx.draw_networkx_nodes(g, pos, nodelist=tgt_tf_nodes, ax=ax, node_shape='s',
                           node_size=GENE_NODE_SIZE, node_color=[color[n] for n in tgt_tf_nodes],
                           edgecolors='black', linewidths=0.5)
    nx.draw_networkx_nodes(g, pos, nodelist=tf_nodes, ax=ax, node_shape='s', node_size=TF_NODE_SIZE,
                           node_color=[color[n] for n in tf_nodes], edgecolors='black',
                           linewidths=0.8)
    # only TF names are labelled (enriched TFs and target TFs); plain genes stay unlabelled
    texts = {}
    texts.update(nx.draw_networkx_labels(g, pos, labels={n: n for n in tf_nodes}, ax=ax,
                                         font_size=9, font_weight='bold',
                                         verticalalignment='bottom'))
    texts.update(nx.draw_networkx_labels(g, pos, labels={n: n for n in tgt_tf_nodes}, ax=ax,
                                         font_size=7, font_weight='bold', verticalalignment='top'))
    for t in texts.values():
        t.set_zorder(LABEL_ZORDER)
    ax.set_rasterization_zorder(RASTER_ZORDER)

    handles = ([Line2D([], [], marker='s', ls='', mfc='gray', mec='black', ms=15, label='TF (+/+ regulon)'),
                Line2D([], [], marker='o', ls='', mfc='gray', mec='white', ms=10, label='target gene'),
                Line2D([], [], marker='s', ls='', mfc='gray', mec='black', ms=10,
                       label='target gene that is a TF (not enriched)'),
                Line2D([], [], color='gray', lw=EDGE_WIDTH, marker='>', ms=6, label='TF -> target (loop: self-regulation)')]
               + [Line2D([], [], marker='o', ls='', mfc=c, mec=c, ms=8, label=f'L2/3 {a}')
                  for a, c in ARCH_COLORS.items()])
    ax.legend(handles=handles, loc='upper left', frameon=False, fontsize=9)
    ax.set_title(f'Yoo25 L2/3 enriched {SIGN} regulons -> archetype marker targets\n'
                 f'{len(tf_nodes)} TFs, {len(gene_nodes)} genes ({len(tgt_tf_nodes)} TFs), '
                 f'{g.number_of_edges()} edges ({n_tf_edges} TF -> TF); '
                 'layout: free force-directed (all nodes move)',
                 fontsize=11)
    ax.axis('off')
    fig.tight_layout()
    fig.savefig(OUT_PDF, dpi=PDF_RASTER_DPI)
    plt.close(fig)
    print(f'written {OUT_PDF}')


if __name__ == '__main__':
    plot(*load_graph())
