"""Regulon x archetype enrichment (expression-stratified null) and regulon target activity.

Enrichment: does an archetype's marker set M overlap a regulon's target set T more than expected?
Universe U = the tested expression genes of the adult L2/3 cells (shared expressed / nonzero-variance
genes across cheng22 + yoo25; 3_regulon/01). A plain Fisher test assumes M and T are exchangeable draws
from U; they are not -- both are biased toward highly expressed genes. The null here draws a pseudo-marker set M*
by sampling, within each expression stratum b, exactly m_b = |M & b| genes uniformly without
replacement from the n_b = |U & b| genes in that stratum. With t_b = |T & b|:

    |M* & T| = sum_b |M*_b & T|,   |M*_b & T| ~ Hypergeometric(n_b, t_b, m_b), independent

so the null is a convolution of per-stratum hypergeometrics, available in closed form:

    exp   = sum_b m_b t_b / n_b
    var   = sum_b m_b (t_b/n_b) (1 - t_b/n_b) (n_b - m_b) / (n_b - 1)
    pmf   = conv_b hypergeom.pmf(0..min(m_b,t_b); n_b, t_b, m_b)
    p     = sum_{k >= x} pmf[k]

The Fisher test is the one-stratum special case; its `log2_or` / `pval` / `fdr` are carried along.
About half of the pairs have overlap 0, where the exact one-sided p is 1.0; `p_mid` (mid-p) is
therefore the BH input (`fdr_strat`), and the strict `p_strat` is kept.

Target activity (`load_expression`, `regulon_score`): per-cell log2(CP10k+1) of a regulon's
targets, each gene scaled to [0, 1] by a high percentile, averaged over the targets.
"""

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import scipy.stats
from statsmodels.stats.multitest import multipletests

from .expr import log2_cp10k, pctile_score

MIN_REGULON_GENES = 5    # drop regulons with fewer than this many targets in-universe
# 20 quantile strata of mean log2-CP10k expression
N_BINS = 20
BIN_COVARIATE = 'mean'
# the closed form is checked against brute-force resampling on one pair; a unit test for the
# convolution indexing, not a way of computing anything
MC_CHECK_DRAWS = 20000
MC_CHECK_RTOL_MEAN = 0.05
MC_CHECK_RTOL_TAIL = 0.10


def log2_odds_ratio(x, M, T, N):
    """Haldane-Anscombe corrected log2 odds ratio for the 2x2 [[x,M-x],[T-x,N-M-T+x]]."""
    a, b, c, d = x + 0.5, (M - x) + 0.5, (T - x) + 0.5, (N - M - T + x) + 0.5
    return float(np.log2((a * d) / (b * c)))


def make_strata(stats, universe_order, n_bins, covariate):
    """Quantile strata of `covariate` over the universe, as an int code per gene.

    Ranked before qcut: the raw values have heavy ties at the low end, on which qcut
    produces unequal and sometimes empty bins.
    """
    v = stats.loc[universe_order, covariate].rank(method='first')
    bins = pd.qcut(v, n_bins, labels=False).values.astype(np.int64)
    counts = np.bincount(bins, minlength=n_bins)
    assert (counts > 0).all(), \
        f'empty stratum with n_bins={n_bins} on {covariate}: counts={counts.tolist()}'
    return bins, counts


def stratified_null(mb, tb, nb, need_pmf=True):
    """Exact null of |M* & T| under stratum-wise resampling of the marker side.

    Returns (mean, sd, pmf). `pmf` is the convolution of the per-stratum hypergeometrics;
    strata with no markers or no targets contribute a point mass at 0 and are skipped.
    """
    assert (mb <= nb).all() and (tb <= nb).all(), 'stratum count exceeds stratum size'
    frac = tb / nb
    mean = float(np.sum(mb * frac))
    var = float(np.sum(mb * frac * (1.0 - frac) * (nb - mb) / np.maximum(nb - 1, 1)))
    if not need_pmf:
        return mean, np.sqrt(var), None

    pmf = np.array([1.0], dtype=np.float64)
    for b in np.nonzero((mb > 0) & (tb > 0))[0]:
        k = np.arange(min(mb[b], tb[b]) + 1)
        pmf = np.convolve(pmf, scipy.stats.hypergeom.pmf(k, nb[b], tb[b], mb[b]))
    assert abs(pmf.sum() - 1.0) < 1e-9, \
        f'convolved pmf lost mass: sum={pmf.sum():.12f} (repeated convolution underflow?)'
    return mean, np.sqrt(var), pmf


def tail_and_midp(pmf, x):
    """One-sided P(X >= x) and its mid-p variant P(X > x) + 0.5 P(X == x)."""
    if x >= len(pmf):
        return 0.0, 0.0
    upper = float(pmf[x:].sum())
    return upper, float(upper - 0.5 * pmf[x])


def score_pairs(layer, bins, counts, M_idx, T_idx, arch_labels):
    """All (archetype, regulon) pairs under one stratification."""
    n_bins = len(counts)
    mb_by_a = {a: np.bincount(bins[idx], minlength=n_bins) for a, idx in M_idx.items()}
    rows = []
    for r, t_idx in T_idx.items():
        tb = np.bincount(bins[t_idx], minlength=n_bins)
        t_set = set(t_idx.tolist())
        for a in arch_labels:
            x = len(t_set & set(M_idx[a].tolist()))
            mean, sd, pmf = stratified_null(mb_by_a[a], tb, counts)
            row = dict(layer=layer, archetype=a, regulon=r, overlap=x, exp=mean, null_sd=sd,
                       log2_enr=float(np.log2((x + 0.5) / (mean + 0.5))),
                       z=float((x - mean) / sd) if sd > 0 else np.nan)
            row['p_strat'], row['p_mid'] = tail_and_midp(pmf, x)
            rows.append(row)
    return pd.DataFrame(rows)


def mc_check(bins, counts, m_idx, t_idx, rng):
    """Brute-force the stratified resampling for one pair; guards the convolution indexing."""
    n_bins = len(counts)
    mb = np.bincount(bins[m_idx], minlength=n_bins)
    tb = np.bincount(bins[t_idx], minlength=n_bins)
    mean, sd, pmf = stratified_null(mb, tb, counts)
    order = np.argsort(bins, kind='stable')
    starts = np.searchsorted(bins[order], np.arange(n_bins + 1))
    t_set = set(t_idx.tolist())
    x_ref = len(t_set & set(m_idx.tolist()))
    draws = np.empty(MC_CHECK_DRAWS, dtype=np.int64)
    for i in range(MC_CHECK_DRAWS):
        tot = 0
        for b in range(n_bins):
            if mb[b] == 0:
                continue
            pool = order[starts[b]:starts[b + 1]]
            tot += len(t_set.intersection(rng.choice(pool, mb[b], replace=False).tolist()))
        draws[i] = tot
    tail_exact, _ = tail_and_midp(pmf, x_ref)
    tail_mc = float((draws >= x_ref).mean())
    print(f'    MC check (B={MC_CHECK_DRAWS}): mean {mean:.4f} vs {draws.mean():.4f}, '
          f'sd {sd:.4f} vs {draws.std():.4f}, P(X>={x_ref}) {tail_exact:.5f} vs {tail_mc:.5f}')
    assert abs(mean - draws.mean()) <= MC_CHECK_RTOL_MEAN * max(mean, 1e-6), \
        f'closed-form mean {mean:.4f} disagrees with resampling {draws.mean():.4f}'
    assert abs(tail_exact - tail_mc) <= MC_CHECK_RTOL_TAIL * max(tail_exact, 1e-3) + 0.005, \
        f'closed-form tail {tail_exact:.5f} disagrees with resampling {tail_mc:.5f}'


def stratified_enrichment(layer, universe, stats, markers, regulons, arch_labels, arch_relabel, rng):
    """Every (archetype, regulon) pair: stratified null + Fisher statistics, BH-FDR, labels.

    markers  : DataFrame with gene / archetype (archetype in arch_labels)
    regulons : DataFrame with regulon / TF / regulation_direction / Gene
    Returns (long table, in-universe marker sets per archetype).
    """
    U = np.array(sorted(universe))
    gi = {g: i for i, g in enumerate(U)}
    stats = stats.loc[U]
    N = len(U)

    M_raw = {a: set(markers.loc[markers['archetype'] == a, 'gene']) & universe for a in arch_labels}
    M_idx = {a: np.array(sorted(gi[g] for g in s), dtype=np.int64) for a, s in M_raw.items()}

    reg_meta = regulons.drop_duplicates('regulon').set_index('regulon')[['TF', 'regulation_direction']]
    T_raw = {r: set(g['Gene']) & universe for r, g in regulons.groupby('regulon')}
    n_total = len(T_raw)
    T_raw = {r: t for r, t in T_raw.items() if len(t) >= MIN_REGULON_GENES}
    T_idx = {r: np.array(sorted(gi[g] for g in s), dtype=np.int64) for r, s in T_raw.items()}
    print(f'    regulons: {len(T_idx)} kept (>= {MIN_REGULON_GENES} in-universe targets) '
          f'of {n_total}; universe N={N}')

    bins, counts = make_strata(stats, U, N_BINS, BIN_COVARIATE)
    long = score_pairs(layer, bins, counts, M_idx, T_idx, arch_labels)

    long['n_markers'] = long['archetype'].map({a: len(s) for a, s in M_raw.items()})
    long['n_targets'] = long['regulon'].map({r: len(s) for r, s in T_raw.items()})
    long['universe'] = N
    long['n_bins'] = N_BINS
    long['bin_covariate'] = BIN_COVARIATE
    long['TF'] = long['regulon'].map(reg_meta['TF'])
    long['regulation_direction'] = long['regulon'].map(reg_meta['regulation_direction'])
    # Fisher (single-stratum) statistics, for side-by-side comparison
    long['log2_or'] = [log2_odds_ratio(x, m, t, N) for x, m, t in
                       zip(long['overlap'], long['n_markers'], long['n_targets'])]
    long['pval'] = [scipy.stats.fisher_exact([[x, m - x], [t - x, N - m - t + x]],
                                             alternative='greater')[1]
                    for x, m, t in zip(long['overlap'], long['n_markers'], long['n_targets'])]
    long['fdr'] = multipletests(long['pval'].values, method='fdr_bh')[1]
    long['fdr_strat'] = multipletests(long['p_mid'].values, method='fdr_bh')[1]
    long['arch_letter'] = long['archetype'].map(arch_relabel)

    cols = ['layer', 'archetype', 'arch_letter', 'regulon', 'TF', 'regulation_direction',
            'overlap', 'n_markers', 'n_targets', 'universe', 'n_bins', 'bin_covariate',
            'exp', 'null_sd', 'log2_enr', 'z', 'p_strat', 'p_mid', 'fdr_strat',
            'log2_or', 'pval', 'fdr']
    long = long[cols]

    # one convolution-vs-resampling check, on the pair with the largest overlap
    top = long.sort_values('overlap', ascending=False).iloc[0]
    mc_check(bins, counts, M_idx[top['archetype']], T_idx[top['regulon']], rng)
    return long, M_raw


# --------------------------------------------------------------------------- regulon tables

REGULON_COLS = ['regulon', 'TF', 'regulation_direction', 'source', 'Gene']


def merge_direct_extended(df):
    """SCENIC+ (regulon, target) rows from both annotations -> one table: a regulon keeps its DIRECT
    targets if it has any, else its EXTENDED ones; unique (regulon, Gene) pairs, sorted."""
    direct_keys = set(df.loc[df['source'] == 'direct', 'regulon'])
    keep = (df['source'] == 'direct') | ~df['regulon'].isin(direct_keys)
    table = (df.loc[keep, REGULON_COLS]
             .drop_duplicates(subset=['regulon', 'Gene'])
             .sort_values(['regulon', 'Gene'])
             .reset_index(drop=True))
    n_direct = table.loc[table['source'] == 'direct', 'regulon'].nunique()
    if n_direct != len(direct_keys):
        raise ValueError(f'direct regulon count mismatch: {n_direct} != {len(direct_keys)}')
    return table


# --------------------------------------------------------------------------- target activity

NORM_MAX_PCTILE = 99


def load_expression(h5ad_path, genes):
    """log2(CP10k+1) of `genes` in every Yoo25 L2/3 cell, plus (Age, Sample, day, condition).

    `h5ad_path` is the 0_prep/01 extract (raw counts in .X, library size = sum over all genes in obs lib).
    """
    src = ad.read_h5ad(h5ad_path)
    present = [g for g in genes if g in src.var_names]
    print(f'{len(present)} of {len(genes)} target genes present')
    X = sp.csr_matrix(src.X)[:, src.var_names.get_indexer(present)].toarray()
    L = log2_cp10k(X, src.obs['lib'].to_numpy())
    cells = (src.obs[['Age', 'Sample', 'day', 'condition']]
             .astype({'Age': str, 'Sample': str, 'condition': str}).reset_index(drop=True))
    print(f'{len(cells)} L2/3 cells; per Age: {cells["Age"].value_counts().to_dict()}')
    return pd.DataFrame(L, columns=present), cells


def regulon_score(tf, genes, L):
    """Per-cell score: each target / its NORM_MAX_PCTILE percentile over L's cells, capped at 1,
    mean over targets. Returns (score, n_targets_used)."""
    used = [g for g in genes if g in L.columns]
    if not used:
        raise ValueError(f'{tf}: none of its {len(genes)} targets is in the data')
    return pctile_score(L[used].to_numpy(), NORM_MAX_PCTILE), len(used)
