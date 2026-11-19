"""TF -> target-gene network of the L2/3-archetype-enriched Yoo25 regulons (activating, +/+): nodes, edges, layout.

Nodes
  - TF (square): every TF of 03's selection (+/+ regulon starred in at least one L2/3 archetype:
    fdr_strat < 0.05, log2_enr > 1, overlap >= 5), coloured by its strongest archetype (max log2_enr).
  - gene (circle): an L2/3 archetype marker (NR union, 03) that is a target of at least one selected TF,
    coloured by the archetype it marks. Markers are unique to one archetype. A target that is
    itself a selected TF is NOT a separate node: it is merged into that TF's square, so TF -> TF
    regulation lands on the square.
  - target TF (is_tf): a target gene that is itself a TF -- has a regulon (any sign) in the Yoo25 L2/3
    SCENIC+ output (01) -- but is not one of the selected TFs.
Edges: TF -> gene for every SCENIC+ target (01) that is an L2/3 marker of ANY archetype, so
cross-archetype regulation is visible; self-regulation is a loop.
Layout: free force-directed (every node moves), seeded from a fixed start -- TFs on a unit ring, each
archetype on its own sector centred on ARCH_ANGLE, spread over TF_ARC_DEG in descending log2_enr; genes
inside -- so the result is deterministic; then rotated so the A' -> C' TF centroids are horizontal (A'
left) with B' below. Figure: 09.

Reads:
  results/3_regulon/03.enriched_regulon_selection.tsv
  results/3_regulon/01.yoo25_l23_regulon_targets.tsv
  results/3_regulon/03.l23_nr_union_markers.tsv
Outputs:
  results/3_regulon/04.regulon_network_edges.tsv   (TF, Gene, tf_arch, gene_arch, target_is_tf)
  results/3_regulon/04.regulon_network_nodes.tsv   (node, kind, arch, is_tf, x, y; in graph order)
"""

import os

import numpy as np
import pandas as pd
import networkx as nx

from l23.constants import ARCH_COLORS

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES_DIR = os.path.join(REPO_ROOT, 'results', '3_regulon')

INPUT_SELECTION = os.path.join(RES_DIR, '03.enriched_regulon_selection.tsv')
INPUT_REGULONS = os.path.join(RES_DIR, '01.yoo25_l23_regulon_targets.tsv')
INPUT_MARKERS = os.path.join(RES_DIR, '03.l23_nr_union_markers.tsv')
OUT_EDGES = os.path.join(RES_DIR, '04.regulon_network_edges.tsv')
OUT_NODES = os.path.join(RES_DIR, '04.regulon_network_nodes.tsv')

# initial positions and free force-directed layout
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


def build_edges(tfs, targets, gene_arch):
    """Unique TF -> marker-target pairs of the selected regulons."""
    keep = targets['regulon'].isin(tfs['regulon']) & targets['Gene'].isin(gene_arch.index)
    edges = targets.loc[keep, ['TF', 'Gene']].drop_duplicates().reset_index(drop=True)
    tf_arch = tfs.set_index('TF')['arch_letter']
    edges['tf_arch'] = edges['TF'].map(tf_arch)
    edges['gene_arch'] = edges['Gene'].map(gene_arch)
    edges['target_is_tf'] = edges['Gene'].isin(tf_arch.index)
    return edges


def build_graph(edges, tf_arch, tf_universe):
    """A target that is a selected TF shares the TF's node id, i.e. it is merged into the square;
    genes in tf_universe are flagged is_tf (target TFs)."""
    g = nx.DiGraph()
    for tf, arch in tf_arch.items():
        g.add_node(tf, kind='tf', arch=arch, is_tf=True)
    genes = edges.loc[~edges['target_is_tf'], ['Gene', 'gene_arch']].drop_duplicates('Gene')
    for gene, arch in zip(genes['Gene'], genes['gene_arch']):
        g.add_node(gene, kind='gene', arch=arch, is_tf=gene in tf_universe)
    g.add_edges_from(zip(edges['TF'], edges['Gene']))
    isolated = [n for n in g if g.degree(n) == 0]
    if isolated:
        raise ValueError(f'selected TFs with no targets: {isolated}')
    return g


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
        if d['kind'] == 'gene':
            t = np.deg2rad(ARCH_ANGLE[d['arch']])
            init[n] = GENE_INIT_R * np.array([np.cos(t), np.sin(t)]) + rng.normal(scale=0.15, size=2)
    return init


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


def main():
    # TFs in name order: the node order of the graph, which the layout depends on
    tfs = pd.read_csv(INPUT_SELECTION, sep='\t').sort_values('TF').reset_index(drop=True)
    targets = pd.read_csv(INPUT_REGULONS, sep='\t')
    markers = pd.read_csv(INPUT_MARKERS, sep='\t')

    for name, labels in [('selection', tfs['arch_letter']), ('markers', markers['arch_letter'])]:
        if set(labels) != set(ARCH_COLORS):
            raise ValueError(f'unexpected archetype labels in {name}: {sorted(set(labels))}')
    if markers['gene'].duplicated().any():
        raise ValueError('L2/3 marker genes are expected to be unique to one archetype')
    gene_arch = markers.set_index('gene')['arch_letter']
    for tf, starred in zip(tfs['TF'], tfs['starred_in']):
        if ',' in starred:
            print(f'  {tf}: starred in {starred}, coloured by {tfs.set_index("TF").loc[tf, "arch_letter"]}')
    print(f'{len(tfs)} enriched TFs: '
          + ', '.join(f'{a}={n}' for a, n in tfs['arch_letter'].value_counts().sort_index().items()))

    edges = build_edges(tfs, targets, gene_arch)
    edges.to_csv(OUT_EDGES, sep='\t', index=False)
    print(f'{len(edges)} edges -> {edges["Gene"].nunique()} marker genes; written {OUT_EDGES}')
    print(pd.crosstab(edges['tf_arch'], edges['gene_arch']).to_string())
    tt = edges[edges['target_is_tf']]
    print(f'{len(tt)} TF -> TF edges ({(tt["TF"] == tt["Gene"]).sum()} self-loops) onto '
          f'{tt["Gene"].nunique()} TFs: {sorted(tt["Gene"].unique())}')

    tf_universe = set(targets['TF'])   # every TF with an L2/3 regulon (01), any sign
    g = build_graph(edges, tfs.set_index('TF')['arch_letter'], tf_universe)
    pos = free_layout(g, tfs, initial_positions(g, tfs))
    nodes = pd.DataFrame([{'node': n, **d, 'x': pos[n][0], 'y': pos[n][1]} for n, d in g.nodes(data=True)])
    nodes.to_csv(OUT_NODES, sep='\t', index=False)
    print(f'{len(nodes)} nodes; written {OUT_NODES}')


if __name__ == '__main__':
    main()
