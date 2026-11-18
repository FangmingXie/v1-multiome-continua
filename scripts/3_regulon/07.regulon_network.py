"""TF -> target-gene network of L2/3-archetype-enriched Yoo25 regulons (activating, +/+).

Nodes
  - TF (square): every Yoo25 L2/3 `+/+` regulon that is significantly enriched in at least
    one L2/3 archetype in 04's stratified test, under the star rule (fdr_strat < 0.05, log2_enr > 1, overlap >= 5). A TF starred in several
    archetypes is coloured by its strongest one (max log2_enr).
  - gene (circle): an L2/3 archetype marker (NR union, 04) that is a target of at least one selected TF,
    coloured by the archetype it marks. Markers are unique to one archetype. A target that is
    itself a selected TF is NOT drawn as a circle: it is merged into that TF's square (which
    keeps its enrichment colour), so TF -> TF regulation lands on the square.
Edges
  - TF -> gene for every SCENIC+ target (01) that is an L2/3 marker of ANY archetype, so
    cross-archetype regulation is visible. All edges are drawn alike (thin, arrowed, in the
    TF's colour); self-regulation is a loop on the square.
  - target TF (small square, gene-dot size, always labelled): a target gene that is itself a TF
    -- i.e. has a regulon (any sign) in the Yoo25 L2/3 SCENIC+ output (01) -- but is not one of
    the selected enriched TFs. Coloured like any gene (archetype, or grey if not a marker).
Labels: TF names only (enriched TFs and target TFs); genes are unlabelled.
PDF: the network (edges, nodes) is one rasterized image; all text stays vector.

All-targets version (OUT_ALL_*): same TFs, but every SCENIC+ target of their regulons is a
node; targets that are not L2/3 archetype markers are grey (unlabelled). Free layout only.

Reads:
  results/3_regulon/04.l23_stratified_enrichment.tsv
  results/3_regulon/01.yoo25_l23_regulon_targets.tsv
  results/3_regulon/04.l23_nr_union_markers.tsv
Outputs:
  results/3_regulon/07.regulon_network_edges.tsv
  figures/3_regulon/07.regulon_network_free.{pdf,png}
  results/3_regulon/07.regulon_network_alltargets_edges.tsv
  figures/3_regulon/07.regulon_network_alltargets.{pdf,png}
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

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES_DIR = os.path.join(REPO_ROOT, 'results', '3_regulon')
FIG_DIR = os.path.join(REPO_ROOT, 'figures', '3_regulon')

INPUT_ENRICHMENT = os.path.join(RES_DIR, '04.l23_stratified_enrichment.tsv')
INPUT_REGULONS = os.path.join(RES_DIR, '01.yoo25_l23_regulon_targets.tsv')
INPUT_MARKERS = os.path.join(RES_DIR, '04.l23_nr_union_markers.tsv')
OUT_EDGES = os.path.join(RES_DIR, '07.regulon_network_edges.tsv')
OUT_FREE_PDF = os.path.join(FIG_DIR, '07.regulon_network_free.pdf')
OUT_FREE_PNG = os.path.join(FIG_DIR, '07.regulon_network_free.png')
OUT_ALL_EDGES = os.path.join(RES_DIR, '07.regulon_network_alltargets_edges.tsv')
OUT_ALL_PDF = os.path.join(FIG_DIR, '07.regulon_network_alltargets.pdf')
OUT_ALL_PNG = os.path.join(FIG_DIR, '07.regulon_network_alltargets.png')

SIGN = '+/+'
# star rule (as 04)
STAR_FDR = 0.05
STAR_LOG2ENR = 1.0
MIN_OVERLAP = 5
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

ARCH_COLORS = {"A'": 'C0', "B'": 'C1', "C'": 'C2'}
OTHER_COLOR = 'lightgray'   # all-targets version: targets that are not archetype markers
OTHER_INIT_SCALE = 0.4      # ... start near the centre with this jitter
# Free force-directed layout (every node moves). It is seeded from a fixed start -- TFs on a
# unit ring, each archetype on its own sector centred on ARCH_ANGLE, spread over TF_ARC_DEG in
# descending log2_enr; genes inside -- so the result is deterministic.
ARCH_ANGLE = {"A'": 90, "B'": 210, "C'": 330}
TF_ARC_DEG = 100
GENE_INIT_R = 0.6     # genes start inside the ring, near their own archetype's sector
LAYOUT_SEED = 0
LAYOUT_K = None       # networkx default 1/sqrt(n_nodes); output rescaled to [-1, 1]
LAYOUT_ITER = 500
GRAVITY_NODE = '__gravity__'
GRAVITY_W = 1.0       # pull of the invisible hub on each TF; 0 would reproduce plain FR
# post-hoc orientation: A' -> C' TF centroids horizontal (A' left), B' centroid below that line
ORIENT_LEFT, ORIENT_RIGHT, ORIENT_BELOW = "A'", "C'", "B'"


os.makedirs(RES_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)


def select_tfs(enr):
    """One row per enriched +/+ TF, assigned to its max-log2_enr starred archetype."""
    star = enr[(enr['regulation_direction'] == SIGN)
               & (enr['fdr_strat'] < STAR_FDR)
               & (enr['log2_enr'] > STAR_LOG2ENR)
               & (enr['overlap'] >= MIN_OVERLAP)]
    if star.empty:
        raise ValueError('no enriched +/+ regulons under the star rule')
    tfs = (star.sort_values('log2_enr', ascending=False)
               .groupby('TF', as_index=False).first()
               [['TF', 'regulon', 'arch_letter', 'log2_enr', 'fdr_strat']])
    multi = star.groupby('TF')['arch_letter'].apply(sorted)
    multi = multi[multi.str.len() > 1]
    for tf, letters in multi.items():
        print(f'  {tf}: starred in {letters}, coloured by '
              f'{tfs.set_index("TF").loc[tf, "arch_letter"]}')
    return tfs


def main():
    enr = pd.read_csv(INPUT_ENRICHMENT, sep='\t')
    targets = pd.read_csv(INPUT_REGULONS, sep='\t')
    markers = pd.read_csv(INPUT_MARKERS, sep='\t')

    for name, labels in [('enrichment', enr['arch_letter']), ('markers', markers['arch_letter'])]:
        if set(labels) != set(ARCH_COLORS):
            raise ValueError(f'unexpected archetype labels in {name}: {sorted(set(labels))}')
    if markers['gene'].duplicated().any():
        raise ValueError('L2/3 marker genes are expected to be unique to one archetype')
    gene_arch = markers.set_index('gene')['arch_letter']

    tfs = select_tfs(enr)
    print(f'{len(tfs)} enriched {SIGN} TFs: '
          + ', '.join(f'{a}={n}' for a, n in tfs['arch_letter'].value_counts().sort_index().items()))

    tf_arch = tfs.set_index('TF')['arch_letter']

    edges = build_edges(tfs, targets, gene_arch, markers_only=True)
    edges.to_csv(OUT_EDGES, sep='\t', index=False)
    print(f'{len(edges)} edges -> {edges["Gene"].nunique()} marker genes; written {OUT_EDGES}')
    print(pd.crosstab(edges['tf_arch'], edges['gene_arch']).to_string())
    tt = edges[edges['target_is_tf']]
    print(f'{len(tt)} TF -> TF edges ({(tt["TF"] == tt["Gene"]).sum()} self-loops) onto '
          f'{tt["Gene"].nunique()} TFs: {sorted(tt["Gene"].unique())}')

    tf_universe = set(targets['TF'])   # every TF with an L2/3 regulon (01), any sign
    g = build_graph(edges, tf_arch, tf_universe)
    plot(g, free_layout(g, tfs, initial_positions(g, tfs)), OUT_FREE_PDF, OUT_FREE_PNG,
         'archetype marker targets', 'free force-directed (all nodes move)')

    edges_all = build_edges(tfs, targets, gene_arch, markers_only=False)
    edges_all.to_csv(OUT_ALL_EDGES, sep='\t', index=False)
    print(f'all targets: {len(edges_all)} edges -> {edges_all["Gene"].nunique()} genes '
          f'({edges_all.loc[edges_all["gene_arch"].notna(), "Gene"].nunique()} markers); '
          f'written {OUT_ALL_EDGES}')
    g_all = build_graph(edges_all, tf_arch, tf_universe)
    plot(g_all, free_layout(g_all, tfs, initial_positions(g_all, tfs)), OUT_ALL_PDF, OUT_ALL_PNG,
         'all targets (grey: not an archetype marker)', 'free force-directed (all nodes move)')


def build_edges(tfs, targets, gene_arch, markers_only):
    """Unique TF -> target pairs of the selected regulons; gene_arch is NaN for non-markers."""
    keep = targets['regulon'].isin(tfs['regulon'])
    if markers_only:
        keep &= targets['Gene'].isin(gene_arch.index)
    edges = targets.loc[keep, ['TF', 'Gene']].drop_duplicates().reset_index(drop=True)
    tf_arch = tfs.set_index('TF')['arch_letter']
    edges['tf_arch'] = edges['TF'].map(tf_arch)
    edges['gene_arch'] = edges['Gene'].map(gene_arch)
    edges['target_is_tf'] = edges['Gene'].isin(tf_arch.index)
    return edges


def build_graph(edges, tf_arch, tf_universe):
    """A target that is a selected TF shares the TF's node id, i.e. it is merged into the square.
    Non-marker genes get arch=None; genes in tf_universe are flagged is_tf (target TFs)."""
    g = nx.DiGraph()
    for tf, arch in tf_arch.items():
        g.add_node(tf, kind='tf', arch=arch)
    genes = edges.loc[~edges['target_is_tf'], ['Gene', 'gene_arch']].drop_duplicates('Gene')
    for gene, arch in zip(genes['Gene'], genes['gene_arch']):
        g.add_node(gene, kind='gene', arch=None if pd.isna(arch) else arch,
                   is_tf=gene in tf_universe)
    g.add_edges_from(zip(edges['TF'], edges['Gene']))
    isolated = [n for n in g if g.degree(n) == 0]
    if isolated:
        raise ValueError(f'selected TFs with no targets: {isolated}')
    return g


def free_layout(g, tfs, init):
    """FR shape is invariant to k and to a uniform edge weight (both only rescale), so cluster
    gaps are closed with a "gravity" node: an invisible hub tied to every TF with weight
    GRAVITY_W, dropped after layout. It pulls the TF clusters (and their targets) together."""
    ug = g.to_undirected()
    nx.set_edge_attributes(ug, 1.0, 'weight')
    ug.add_edges_from(((GRAVITY_NODE, tf) for tf in tfs['TF']), weight=GRAVITY_W)
    pos = nx.spring_layout(ug, pos={**init, GRAVITY_NODE: np.zeros(2)}, k=LAYOUT_K,
                           iterations=LAYOUT_ITER, seed=LAYOUT_SEED, weight='weight')
    del pos[GRAVITY_NODE]
    return orient(pos, tfs)


def orient(pos, tfs):
    """Rigidly rotate the layout so the ORIENT_LEFT -> ORIENT_RIGHT TF centroids lie on a
    horizontal line, then mirror top-bottom if needed so ORIENT_BELOW's centroid sits below it.
    Distances are preserved; FR layouts have no intrinsic orientation."""
    cen = {a: np.mean([pos[tf] for tf in grp['TF']], axis=0)
           for a, grp in tfs.groupby('arch_letter')}
    d = cen[ORIENT_RIGHT] - cen[ORIENT_LEFT]
    t = -np.arctan2(d[1], d[0])
    rot = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    origin = (cen[ORIENT_LEFT] + cen[ORIENT_RIGHT]) / 2
    flip = np.diag([1.0, -1.0 if (rot @ (cen[ORIENT_BELOW] - origin))[1] > 0 else 1.0])
    m = flip @ rot
    print(f'  orient: rotated {np.rad2deg(t):+.0f} deg'
          f'{", mirrored top-bottom" if flip[1, 1] < 0 else ""}')
    return {n: m @ (xy - origin) for n, xy in pos.items()}


def initial_positions(g, tfs):
    """TFs on their archetype's ring sector (descending log2_enr); genes inside, near their own
    archetype's sector, with seeded jitter."""
    init = {}
    for arch, grp in tfs.sort_values('log2_enr', ascending=False).groupby('arch_letter'):
        offs = np.linspace(-TF_ARC_DEG / 2, TF_ARC_DEG / 2, len(grp)) if len(grp) > 1 else [0.0]
        for tf, off in zip(grp['TF'], offs):
            t = np.deg2rad(ARCH_ANGLE[arch] + off)
            init[tf] = np.array([np.cos(t), np.sin(t)])
    rng = np.random.default_rng(LAYOUT_SEED)
    for n, d in g.nodes(data=True):
        if d['kind'] == 'gene' and d['arch'] is None:
            init[n] = rng.normal(scale=OTHER_INIT_SCALE, size=2)
        elif d['kind'] == 'gene':
            t = np.deg2rad(ARCH_ANGLE[d['arch']])
            init[n] = GENE_INIT_R * np.array([np.cos(t), np.sin(t)]) + rng.normal(scale=0.15, size=2)
    return init


def plot(g, pos, out_pdf, out_png, target_desc, layout_desc):
    tf_nodes = [n for n, d in g.nodes(data=True) if d['kind'] == 'tf']
    # non-markers first so the coloured markers are drawn on top of them
    gene_nodes = sorted((n for n, d in g.nodes(data=True) if d['kind'] == 'gene'),
                        key=lambda n: g.nodes[n]['arch'] is not None)
    other_nodes = [n for n in gene_nodes if g.nodes[n]['arch'] is None]
    tgt_tf_nodes = [n for n in gene_nodes if g.nodes[n]['is_tf']]
    dot_nodes = [n for n in gene_nodes if not g.nodes[n]['is_tf']]
    color = {n: OTHER_COLOR if d['arch'] is None else ARCH_COLORS[d['arch']]
             for n, d in g.nodes(data=True)}

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
                  for a, c in ARCH_COLORS.items()]
               + ([Line2D([], [], marker='o', ls='', mfc=OTHER_COLOR, mec=OTHER_COLOR, ms=8,
                          label='not an archetype marker')] if other_nodes else []))
    ax.legend(handles=handles, loc='upper left', frameon=False, fontsize=9)
    ax.set_title(f'Yoo25 L2/3 enriched {SIGN} regulons -> {target_desc}\n'
                 f'{len(tf_nodes)} TFs, {len(gene_nodes)} genes ({len(other_nodes)} non-marker, {len(tgt_tf_nodes)} TFs), {g.number_of_edges()} edges ({n_tf_edges} TF -> TF); '
                 f'layout: {layout_desc}',
                 fontsize=11)
    ax.axis('off')
    fig.tight_layout()
    fig.savefig(out_pdf, dpi=PDF_RASTER_DPI)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f'written {out_pdf}\nwritten {out_png}')


if __name__ == '__main__':
    main()
