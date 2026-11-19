"""Expression-stratified regulon x archetype enrichment, L2/3, on the developmental archetype genes.

Test (src/l23/regulon.py): each Yoo25 L2/3 regulon's targets (0_prep/03) against each archetype's marker
set, over the adult L2/3 gene universe (01), with the closed-form null stratified by mean expression
(20 strata) and BH over mid-p.

Marker sets: the archetype genes of 2_development/01 (= Table S4), union per archetype over the
normally-reared timepoints P6..P21 (every DR condition excluded). A gene called for more than one
archetype at different timepoints is dropped, so each marker belongs to exactly one archetype.
Archetype letters are carried as the adult PCHA labels (A -> archetype_3, B -> archetype_2,
C -> archetype_1; column archetype), whose display labels are A'/B'/C' (arch_letter).

Selection: the star rule (fdr_strat < 0.05, log2_enr > 1, overlap >= 5) on the +/+ regulons;
each TF is assigned to its max-log2_enr starred archetype.

Reads:
  results/2_development/01.yoo25_l23_abc_genes.tsv
  results/0_prep/03.yoo25_l23_regulon_targets.tsv
  results/3_regulon/01.l23_gene_universe_stats.tsv
Outputs:
  results/3_regulon/02.l23_nr_union_markers.tsv
  results/3_regulon/02.l23_stratified_enrichment.tsv
  results/3_regulon/02.enriched_regulon_selection.tsv
"""

import os

import numpy as np
import pandas as pd

from l23.regulon import stratified_enrichment

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_GENES    = os.path.join(REPO_ROOT, 'results', '2_development', '01.yoo25_l23_abc_genes.tsv')
INPUT_REGULONS = os.path.join(REPO_ROOT, 'results', '0_prep', '03.yoo25_l23_regulon_targets.tsv')
INPUT_UNIVERSE = os.path.join(REPO_ROOT, 'results', '3_regulon', '01.l23_gene_universe_stats.tsv')
OUT_MARKERS    = os.path.join(REPO_ROOT, 'results', '3_regulon', '02.l23_nr_union_markers.tsv')
OUT_ENRICHMENT = os.path.join(REPO_ROOT, 'results', '3_regulon', '02.l23_stratified_enrichment.tsv')
OUT_SELECTION  = os.path.join(REPO_ROOT, 'results', '3_regulon', '02.enriched_regulon_selection.tsv')

LAYER = 'L2_3'
NR_CONDITIONS = {'P6', 'P8', 'P10', 'P12', 'P14', 'P17', 'P21'}
# archetype letter -> adult PCHA label; ARCH_RELABEL maps those to the display labels
S4_TO_PCHA = {'A': 'archetype_3', 'B': 'archetype_2', 'C': 'archetype_1'}
ARCH_LABELS = ['archetype_1', 'archetype_2', 'archetype_3']
ARCH_RELABEL = {'archetype_1': "C'", 'archetype_2': "B'", 'archetype_3': "A'"}
MIN_MARKERS_IN_UNIVERSE = 20
SIGN = '+/+'
# star rule
STAR_FDR = 0.05
STAR_LOG2ENR = 1.0
MIN_OVERLAP = 5

os.makedirs(os.path.dirname(OUT_MARKERS), exist_ok=True)


def build_markers():
    """NR-timepoint union per archetype, genes in >1 archetype dropped."""
    genes = pd.read_csv(INPUT_GENES, sep='\t', usecols=['cond', 'archetype', 'gene'])
    nr = genes[~genes['cond'].str.endswith('DR')]
    if set(nr['cond']) != NR_CONDITIONS:
        raise ValueError(f'unexpected NR conditions in {INPUT_GENES}: {sorted(set(nr["cond"]))}')
    if set(nr['archetype']) != set(S4_TO_PCHA):
        raise ValueError(f'unexpected archetypes in {INPUT_GENES}: {sorted(set(nr["archetype"]))}')

    per_gene = nr.groupby('gene')['archetype'].agg(lambda x: ''.join(sorted(set(x))))
    multi = per_gene[per_gene.str.len() > 1]
    print(f'NR union: {per_gene.size} genes; dropping {len(multi)} called for >1 '
          f'archetype ({multi.value_counts().to_dict()})')
    keep = nr[~nr['gene'].isin(multi.index)]
    mk = (keep.groupby(['gene', 'archetype']).size().rename('n_conditions').reset_index()
              .rename(columns={'archetype': 'letter'}))
    mk['archetype'] = mk['letter'].map(S4_TO_PCHA)
    mk['arch_letter'] = mk['letter'] + "'"
    mk = mk[['gene', 'archetype', 'arch_letter', 'n_conditions']].sort_values(['arch_letter', 'gene'])
    mk.to_csv(OUT_MARKERS, sep='\t', index=False)
    print(f'  markers per archetype: {mk["arch_letter"].value_counts().sort_index().to_dict()}; '
          f'written {OUT_MARKERS}')
    return mk


def select(enr):
    """One row per starred +/+ TF, at its max-log2_enr starred archetype."""
    star = enr[(enr['regulation_direction'] == SIGN)
               & (enr['fdr_strat'] < STAR_FDR)
               & (enr['log2_enr'] > STAR_LOG2ENR)
               & (enr['overlap'] >= MIN_OVERLAP)]
    sel = (star.sort_values('log2_enr', ascending=False)
               .groupby('TF', as_index=False).first()
               [['TF', 'regulon', 'arch_letter', 'log2_enr', 'fdr_strat', 'overlap',
                 'n_markers', 'n_targets']])
    starred_in = star.groupby('TF')['arch_letter'].apply(lambda x: ','.join(sorted(x)))
    sel['starred_in'] = sel['TF'].map(starred_in)
    return sel.sort_values(['arch_letter', 'log2_enr'], ascending=[True, False])


mk = build_markers()

stats = pd.read_csv(INPUT_UNIVERSE, sep='\t', index_col='gene')
universe = set(stats.index)
print(f'gene universe: {len(universe)} genes')
regulons = pd.read_csv(INPUT_REGULONS, sep='\t')

enr, M_raw = stratified_enrichment(LAYER, universe, stats, mk, regulons, ARCH_LABELS, ARCH_RELABEL,
                                   np.random.default_rng(0))
enr.to_csv(OUT_ENRICHMENT, sep='\t', index=False)
print(f'wrote {OUT_ENRICHMENT} ({len(enr)} pairs)')

n_all = mk['archetype'].value_counts().to_dict()
for a in sorted(M_raw, key=ARCH_RELABEL.get):
    print(f'  {ARCH_RELABEL[a]}: {len(M_raw[a])} of {n_all[a]} markers in the universe')
    if len(M_raw[a]) < MIN_MARKERS_IN_UNIVERSE:
        raise ValueError(f'{ARCH_RELABEL[a]} keeps only {len(M_raw[a])} in-universe markers')

sel = select(enr)
sel.to_csv(OUT_SELECTION, sep='\t', index=False)
print(f'\n{len(sel)} enriched {SIGN} TFs; written {OUT_SELECTION}')
print(sel.to_string(index=False, float_format=lambda v: f'{v:.3g}'))
