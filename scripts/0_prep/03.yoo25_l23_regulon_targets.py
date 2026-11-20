"""Clean regulon -> target-gene table from the Yoo25 SCENIC+ eRegulons of L2/3.

The SCENIC+ output is split into two metadata files, "direct" and "extended". Each row is a
TF -> enhancer-region -> target-gene triplet, and `eRegulon_name` encodes the regulon as
`<TF>_<direct|extended>_<sign>` (e.g. `Arnt_direct_+/+`), where the sign pattern is one of `+/+`,
`+/-`, `-/+`, `-/-` (TF2G / R2G regulation directions). "Direct" uses direct TF-motif
annotations; "extended" uses extended/indirect annotations.

One tidy (long) table of regulons and their unique target genes:
  - A regulon is identified by TF + sign pattern (direct/extended normalized away), e.g.
    `Arnt_+/+`. `Arnt_+/+` and `Arnt_-/+` are distinct regulons.
  - Only regulons with a positive R2G (second) sign are kept (+/+, -/+).
  - If a regulon exists in both direct and extended, keep its DIRECT targets; otherwise keep its
    EXTENDED targets.
  - The same gene appears once per enhancer region, so target pairs are deduped to unique
    (regulon, Gene).

Reads:
  data/raw/regulons/yoo25_L2_3/L2_3_direct_e_regulon_metadata.csv
  data/raw/regulons/yoo25_L2_3/L2_3_extended_e_regulon_metadata.csv
Outputs:
  data/input/03.yoo25_l23_regulon_targets.tsv
"""

import os
import pandas as pd

from l23.regulon import merge_direct_extended

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_DIRECT   = os.path.join(REPO_ROOT, 'data', 'raw', 'regulons', 'yoo25_L2_3', 'L2_3_direct_e_regulon_metadata.csv')
INPUT_EXTENDED = os.path.join(REPO_ROOT, 'data', 'raw', 'regulons', 'yoo25_L2_3', 'L2_3_extended_e_regulon_metadata.csv')
OUT_TSV        = os.path.join(REPO_ROOT, 'data', 'input', '03.yoo25_l23_regulon_targets.tsv')

KEEP_DIRECTIONS = {'+/+', '-/+'}   # positive R2G (second) sign


def load_and_annotate(csv_path, source):
    """Load an eRegulon metadata CSV and add normalized regulon-identity columns."""
    df = pd.read_csv(csv_path, index_col=0)
    # normalized regulon key: drop the direct/extended tag, e.g. Arnt_direct_+/+ -> Arnt_+/+
    df['regulon'] = (df['eRegulon_name']
                     .str.replace('_direct_', '_', regex=False)
                     .str.replace('_extended_', '_', regex=False))
    df['regulation_direction'] = df['eRegulon_name'].str.rsplit('_', n=1).str[-1]
    df['source'] = source
    if not set(df['regulation_direction'].unique()) <= {'+/+', '+/-', '-/+', '-/-'}:
        raise ValueError(f"Unexpected sign pattern(s) in {csv_path}: "
                         f"{sorted(df['regulation_direction'].unique())}")
    return df[df['regulation_direction'].isin(KEEP_DIRECTIONS)].copy()


# multiple-region rows collapse to unique (regulon, target gene) pairs
table = merge_direct_extended(pd.concat([load_and_annotate(INPUT_DIRECT, 'direct'),
                                         load_and_annotate(INPUT_EXTENDED, 'extended')], ignore_index=True))
n_direct = table.loc[table['source'] == 'direct', 'regulon'].nunique()

os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
table.to_csv(OUT_TSV, sep='\t', index=False)

print(f"regulons total         : {table['regulon'].nunique()}")
print(f"  from direct          : {n_direct}")
print(f"  from extended (only) : {table.loc[table['source'] == 'extended', 'regulon'].nunique()}")
print(f"  +/+ (activating)     : {table.loc[table['regulation_direction'] == '+/+', 'regulon'].nunique()}")
print(f"  -/+ (repressing)     : {table.loc[table['regulation_direction'] == '-/+', 'regulon'].nunique()}")
print(f"unique target genes    : {table['Gene'].nunique()}")
print(f"target pairs (rows)    : {len(table)}")
print(f'wrote {OUT_TSV}')
