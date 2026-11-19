"""Mouse vs human L2/3 axes: PCA + varimax per species, and closed-form CCA on ortholog loadings.

Varimax: HVGs by raw variance of the log-normalized matrix -> StandardScaler -> PCA -> Kaiser varimax,
components reordered by descending score variance; per-component variance partition (partial R^2 per
factor) for choosing the cell-type (Gate-A) components.

CCA: an "axis" is a linear mixture of a species' varimax (VX) components; with genes as observations
and VX components as features (rows paired by 1-to-1 orthology), maximizing corr(Xa, Yb) is exactly
CCA, and the canonical correlations are the cosines of the principal angles between the two
gene-loading subspaces. Per pair:
  1. Closed-form CCA spectrum, sum cos^2 (subspace overlap), per-block condition number.
  2. Cross-validation over genes, folds BLOCKED BY CO-EXPRESSION MODULE (random folds leak through
     co-expressed genes and inflate CV); blocked and random, over CV_R seeds.
  3. Bootstrap weight stability; weights interpretable only when the median |cos| >= STABLE_THRESH.
  4. Permutation null (mouse gene labels); z is a standardized effect size on bounded support, NOT a
     normal deviate, and is subordinate to CV because co-expression violates gene exchangeability.
"""

import numpy as np
import pandas as pd
from numpy.linalg import LinAlgError
from scipy.linalg import cholesky, solve_triangular
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.stats import zscore
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression

from .expr import log2_cp10k

# --- CCA parameters ---
N_PERM          = 4000     # permutation null replicates (z quoted to +-0.2)
N_BOOT          = 500      # gene bootstraps for weight stability
CV_K            = 5        # cross-validation folds over genes
CV_R            = 10       # module-assignment seeds (spread of blocked/random CV)
MODULE_DIVISOR  = 8        # n_modules ~ n_genes / MODULE_DIVISOR
BOOT_NCOMP      = 2        # bootstrap stability computed for components 1..2
STABLE_THRESH   = 0.9      # bootstrap median |cos| gate for interpreting weights
SEED            = 0


# ======================================================================================
# PCA + varimax
# ======================================================================================

def varimax(L, gamma=1.0, max_iter=1000, tol=1e-6):
    """Kaiser varimax rotation of loading matrix L (n_vars x n_factors)."""
    n, p = L.shape
    R = np.eye(p)
    for _ in range(max_iter):
        R_old = R.copy()
        for i in range(p - 1):
            for j in range(i + 1, p):
                Lr = L @ R
                u = Lr[:, i] ** 2 - Lr[:, j] ** 2
                v = 2 * Lr[:, i] * Lr[:, j]
                A = u.sum()
                B = v.sum()
                C = (u ** 2 - v ** 2).sum()
                D = 2 * (u * v).sum()
                theta = 0.25 * np.arctan2(
                    D - gamma * 2 * A * B / n,
                    C - gamma * (A ** 2 - B ** 2) / n,
                )
                c, s = np.cos(theta), np.sin(theta)
                Rij = np.eye(p)
                Rij[i, i] = Rij[j, j] = c
                Rij[i, j] = -s
                Rij[j, i] = s
                R = R @ Rij
        if np.max(np.abs(R - R_old)) < tol:
            break
    return R


def pca_varimax(X_norm, gene_names, n_hvg, n_pcs, random_state=0):
    """Top-n_hvg genes by variance of X_norm -> StandardScaler -> PCA -> varimax, VX by descending variance.

    Returns (hvg_names, pc_scores, vx_scores, vx_loadings); loadings are genes x components.
    """
    gene_var = X_norm.var(axis=0)
    hvg_idx = np.argsort(gene_var)[::-1][:n_hvg]
    X_scaled = StandardScaler().fit_transform(X_norm[:, hvg_idx])
    pca = PCA(n_components=n_pcs, random_state=random_state)
    scores = pca.fit_transform(X_scaled)   # (n_cells, n_pcs)
    L = pca.components_.T                  # (n_hvg, n_pcs)
    print(f'  Explained variance ratio: {pca.explained_variance_ratio_.round(3)}')
    R = varimax(L)
    vx_scores = scores @ R
    vx_loadings = L @ R
    order = np.argsort(vx_scores.var(axis=0))[::-1]
    return np.asarray(gene_names)[hvg_idx], scores, vx_scores[:, order], vx_loadings[:, order]


def r2(X, y):
    reg = LinearRegression().fit(X, y)
    ss_res = np.sum((y - reg.predict(X)) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - ss_res / ss_tot


def dummies(values, name):
    """One-hot design matrix; fails fast on all-zero dummy columns (stale categorical levels)."""
    X = pd.get_dummies(values).values.astype(float)
    if (X.sum(axis=0) == 0).any():
        raise ValueError(f'{name} design matrix has {(X.sum(axis=0) == 0).sum()} all-zero dummy column(s) — '
                         f'stale categorical levels survived subsetting')
    return X


def library_size_factor(libsize):
    return zscore(libsize).reshape(-1, 1)


def variance_partition(scores, comp_cols, factors, label):
    """Partial R^2 per component for each factor in `factors` (dict name -> design matrix)."""
    X_full = np.hstack(list(factors.values()))
    print(f'\nVariance partitioning ({label}):')
    rows = []
    for i, col in enumerate(comp_cols):
        y = scores[:, i]
        r2_full = r2(X_full, y)
        partial = {}
        for name in factors:
            X_others = np.hstack([v for k, v in factors.items() if k != name])
            partial[name] = max(r2_full - r2(X_others, y), 0.0)
        residual = max(1.0 - r2_full, 0.0)
        rows.append({'component': col, **partial, 'residual': residual})
        print(f'  {col}  R²_full={r2_full:.3f}  '
              + '  '.join(f'{k}={v:.3f}' for k, v in partial.items())
              + f'  resid={residual:.3f}')
    return pd.DataFrame(rows).set_index('component')


def mouse_log2_cp10k(adata):
    """log2(CP10k + 1) from the integer counts in .raw (cheng22 .X is already normalized)."""
    X_raw = adata.raw[:, adata.var_names].X.toarray().astype(np.float32)
    depths = X_raw.sum(axis=1)
    depths[depths == 0] = 1
    return log2_cp10k(X_raw, depths)


# ======================================================================================
# Orthologs
# ======================================================================================

def load_orthologs(path):
    """1-to-1 human/mouse ortholog pairs (the order of the two drop_duplicates calls matters)."""
    return (pd.read_csv(path, sep='\t')
            .drop_duplicates('human_symbol').drop_duplicates('mouse_symbol'))


# ======================================================================================
# Core CCA (closed form)
# ======================================================================================

def orthonormal(M):
    """Orthonormal basis (QR) for the column span of centered M — genes x k."""
    return np.linalg.qr(M - M.mean(axis=0))[0]


def cca_spectrum(X, Y):
    """Canonical correlations = cos(principal angles) between spans of centered X, Y."""
    s = np.linalg.svd(orthonormal(X).T @ orthonormal(Y), compute_uv=False)
    return np.clip(s, 0.0, 1.0)


def cca_fit(X, Y):
    """Full closed-form CCA. Returns (corrs, A, B): canonical correlations and the weight
    matrices mapping centered X, Y to canonical variates (Xc @ A, Yc @ B)."""
    n = X.shape[0]
    Xc = X - X.mean(axis=0)
    Yc = Y - Y.mean(axis=0)
    Sxx = Xc.T @ Xc / (n - 1)
    Syy = Yc.T @ Yc / (n - 1)
    Sxy = Xc.T @ Yc / (n - 1)
    Lx = cholesky(Sxx, lower=True)                        # Sxx = Lx Lx'
    Ly = cholesky(Syy, lower=True)
    K = solve_triangular(Lx, Sxy, lower=True)             # Lx^-1 Sxy
    K = solve_triangular(Ly, K.T, lower=True).T           # Lx^-1 Sxy Ly^-T
    U, D, Vt = np.linalg.svd(K)
    A = solve_triangular(Lx.T, U, lower=False)            # Lx^-T U
    B = solve_triangular(Ly.T, Vt.T, lower=False)         # Ly^-T V
    return np.clip(D, 0.0, 1.0), A, B


def sumcos2(s):
    """Subspace overlap sum cos^2 = sum s_i^2."""
    return float(np.dot(s, s))


def cond_zscored(M):
    """Condition number of the column-standardized block."""
    Mz = (M - M.mean(axis=0)) / M.std(axis=0)
    return float(np.linalg.cond(Mz))


# Permuting rows of Y commutes with orthonormalization (Qy_perm = P Qy), so each replicate is one
# k x k SVD of Qx' (P Qy). The full spectrum is retained per replicate.
def permutation_null(X, Y, n_perm, seed):
    Qx, Qy = orthonormal(X), orthonormal(Y)
    n, k = Qy.shape[0], min(X.shape[1], Y.shape[1])
    rng = np.random.default_rng(seed)
    null_spec = np.empty((n_perm, k))
    for i in range(n_perm):
        Qyp = Qy[rng.permutation(n)]
        s = np.clip(np.linalg.svd(Qx.T @ Qyp, compute_uv=False), 0.0, 1.0)
        null_spec[i] = s
    null_sumcos2 = (null_spec ** 2).sum(axis=1)
    return null_spec, null_sumcos2


# ======================================================================================
# Co-expression modules and fold assignment for cross-validation
# ======================================================================================

def coexpression_modules(X):
    """Ward-cluster the human loading matrix (genes x kx) into ~ n/MODULE_DIVISOR modules."""
    n = X.shape[0]
    n_modules = max(2, round(n / MODULE_DIVISOR))
    Z = linkage(X, method='ward')
    return fcluster(Z, t=n_modules, criterion='maxclust')


def module_folds(modules, K, seed):
    """Assign whole modules to K folds (genes in a module always share a fold)."""
    rng = np.random.default_rng(seed)
    mods = np.unique(modules)
    perm = rng.permutation(len(mods))
    fold_of_mod = {mods[perm[i]]: i % K for i in range(len(mods))}
    return np.array([fold_of_mod[m] for m in modules])


def random_folds(n, K, seed):
    rng = np.random.default_rng(seed)
    idx = rng.permutation(n)
    f = np.empty(n, dtype=int)
    for i, g in enumerate(idx):
        f[g] = i % K
    return f


def cv_one_assignment(X, Y, fold_of_gene, r, rng):
    """One K-fold pass. Returns (comp_folds K x r signed, sub_folds K, sub_null_folds K)."""
    K = int(fold_of_gene.max()) + 1
    comp_folds, sub_folds, sub_null = [], [], []
    for f in range(K):
        te = fold_of_gene == f
        tr = ~te
        if te.sum() < 3 or tr.sum() <= X.shape[1]:
            continue
        _, A, B = cca_fit(X[tr], Y[tr])
        Xt, Yt = X[te] @ A, Y[te] @ B                     # Pearson centers internally
        cc = np.array([np.corrcoef(Xt[:, i], Yt[:, i])[0, 1] for i in range(r)])
        comp_folds.append(cc)
        Qxt = np.linalg.qr(Xt - Xt.mean(axis=0))[0]
        Qyt = np.linalg.qr(Yt - Yt.mean(axis=0))[0]
        s = np.clip(np.linalg.svd(Qxt.T @ Qyt, compute_uv=False), 0.0, 1.0)
        sub_folds.append(sumcos2(s))
        Qyp = Qyt[rng.permutation(Qyt.shape[0])]          # destroy pairing on held-out genes
        sp = np.clip(np.linalg.svd(Qxt.T @ Qyp, compute_uv=False), 0.0, 1.0)
        sub_null.append(sumcos2(sp))
    return np.array(comp_folds), np.array(sub_folds), np.array(sub_null)


def cross_validate(X, Y, modules, r, seed):
    """Blocked (by module) and random K-fold CV over CV_R seeds.

    Returns dict with component-wise signed means (blocked/random), subspace sum cos^2
    means (blocked/random) and its blocked null, plus seed-0 blocked per-fold comp values.
    """
    n = X.shape[0]
    rng = np.random.default_rng(seed)
    blk_comp, rnd_comp = [], []
    blk_sub, rnd_sub, blk_sub_null, rnd_sub_null = [], [], [], []
    comp_folds0 = None
    for s in range(CV_R):
        cf_b, sf_b, sn_b = cv_one_assignment(X, Y, module_folds(modules, CV_K, seed + s), r, rng)
        cf_r, sf_r, sn_r = cv_one_assignment(X, Y, random_folds(n, CV_K, seed + 1000 + s), r, rng)
        if cf_b.size:
            blk_comp.append(cf_b.mean(axis=0))
            blk_sub.append(sf_b.mean())
            blk_sub_null.append(sn_b.mean())
        if cf_r.size:
            rnd_comp.append(cf_r.mean(axis=0))
            rnd_sub.append(sf_r.mean())
            rnd_sub_null.append(sn_r.mean())
        if s == 0:
            comp_folds0 = cf_b
    sub_blk = float(np.mean(blk_sub))
    sub_blk_null = float(np.mean(blk_sub_null))
    sub_rnd = float(np.mean(rnd_sub))
    sub_rnd_null = float(np.mean(rnd_sub_null))
    return {
        'comp_blocked': np.mean(blk_comp, axis=0),
        'comp_random':  np.mean(rnd_comp, axis=0),
        # raw held-out sum cos^2 carries a random-subspace baseline (~ r^2/n_test); the
        # generalization signal is the excess over the permuted-pairing null, which is the
        # quantity that returns ~0 under a destroyed pairing.
        'sub_blocked_raw': sub_blk, 'sub_blocked_null': sub_blk_null,
        'sub_random_raw':  sub_rnd, 'sub_random_null':  sub_rnd_null,
        'sub_blocked': sub_blk - sub_blk_null,   # excess (headline)
        'sub_random':  sub_rnd - sub_rnd_null,
        'comp_folds0':  comp_folds0,
    }


# ======================================================================================
# Bootstrap weight stability
# ======================================================================================

def _abscos(a, b):
    return abs(float(a @ b) / (np.linalg.norm(a) * np.linalg.norm(b)))


def bootstrap_weights(X, Y, A_full, B_full, n_boot, seed):
    """|cos| between resampled and full-data canonical weight vectors, components 1..BOOT_NCOMP."""
    n = X.shape[0]
    ncomp = min(BOOT_NCOMP, A_full.shape[1])
    rng = np.random.default_rng(seed)
    cos_h = {c: [] for c in range(ncomp)}
    cos_m = {c: [] for c in range(ncomp)}
    n_skip = 0
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        try:
            _, A, B = cca_fit(X[idx], Y[idx])
        except LinAlgError:
            n_skip += 1
            continue
        for c in range(ncomp):
            cos_h[c].append(_abscos(A[:, c], A_full[:, c]))
            cos_m[c].append(_abscos(B[:, c], B_full[:, c]))

    def summarize(store):
        out = {}
        for c, vals in store.items():
            v = np.array(vals)
            out[c] = (float(np.median(v)), float(np.percentile(v, 5)), float((v < STABLE_THRESH).mean()))
        return out

    return summarize(cos_h), summarize(cos_m), n_skip


# ======================================================================================
# Per-pair analysis and output tables
# ======================================================================================

def analyze_pair(X, Y, seed):
    """Steps 1, 2, 4 for one (X, Y) pair. Modules from the human block X."""
    kx, ky = X.shape[1], Y.shape[1]
    r = min(kx, ky)
    spec = cca_spectrum(X, Y)
    obs_sc2 = sumcos2(spec)

    null_spec, null_sc2 = permutation_null(X, Y, N_PERM, seed)
    null_mean = null_spec.mean(axis=0)
    null_sd = null_spec.std(axis=0)
    z = (spec - null_mean) / null_sd
    sc2_mean, sc2_sd = float(null_sc2.mean()), float(null_sc2.std())
    sc2_z = (obs_sc2 - sc2_mean) / sc2_sd
    sc2_p_emp = (1 + int((null_sc2 >= obs_sc2).sum())) / (N_PERM + 1)

    modules = coexpression_modules(X)
    cv = cross_validate(X, Y, modules, r, seed)

    return {
        'n': X.shape[0], 'kx': kx, 'ky': ky, 'r': r,
        'spec': spec, 'sumcos2': obs_sc2, 'sumcos2_frac': obs_sc2 / r,
        'null_mean': null_mean, 'null_sd': null_sd, 'z': z,
        'sc2_null_mean': sc2_mean, 'sc2_null_sd': sc2_sd, 'sc2_z': float(sc2_z),
        'sc2_p_emp': sc2_p_emp,
        'cond_x': cond_zscored(X), 'cond_y': cond_zscored(Y),
        'cv': cv, 'modules': modules,
    }


def sign_fix(a, b):
    """Flip both weight vectors by the sign of a's largest-|.| entry (keeps corr sign)."""
    s = np.sign(a[np.argmax(np.abs(a))])
    if s == 0:
        s = 1.0
    return a * s, b * s


def run_cca(X, Y):
    """Steps 1-4 for one pair; prints a summary. Returns (res, A, B, boot_h, boot_m)."""
    res = analyze_pair(X, Y, SEED)
    print(f'  n={res["n"]}, kx={res["kx"]}, ky={res["ky"]}')
    print(f'  spectrum : {np.round(res["spec"], 4)}')
    print(f'  Σcos²θ   : {res["sumcos2"]:.3f}  (frac {res["sumcos2_frac"]:.3f}, '
          f'z {res["sc2_z"]:.1f}, null {res["sc2_null_mean"]:.3f})')
    print(f'  z/comp   : {np.round(res["z"], 1)}')
    print(f'  CV comp blocked: {np.round(res["cv"]["comp_blocked"], 3)}  '
          f'random: {np.round(res["cv"]["comp_random"], 3)}')
    print(f'  CV Σcos²θ excess blocked {res["cv"]["sub_blocked"]:.3f} '
          f'(raw {res["cv"]["sub_blocked_raw"]:.3f} − null {res["cv"]["sub_blocked_null"]:.3f}) | '
          f'random {res["cv"]["sub_random"]:.3f}')
    _, A, B = cca_fit(X, Y)
    boot_h, boot_m, n_skip = bootstrap_weights(X, Y, A, B, N_BOOT, SEED)
    print(f'  bootstrap ({N_BOOT}, {n_skip} skipped):')
    for c in boot_h:
        print(f'    comp{c+1}: human med {boot_h[c][0]:.3f} p5 {boot_h[c][1]:.3f} '
              f'frac<0.9 {boot_h[c][2]:.2f} | mouse med {boot_m[c][0]:.3f}')
    return res, A, B, boot_h, boot_m


def write_spectrum(res, out_path):
    cv = res['cv']
    r = res['r']
    rows = []
    for i in range(len(res['spec'])):
        folds0 = cv['comp_folds0']
        fold_str = ';'.join(f'{v:.4f}' for v in folds0[:, i]) if (folds0 is not None and i < r) else ''
        rows.append({
            'component': f'CCA{i+1}', 'r': res['spec'][i], 'frac': np.nan,
            'null_mean': res['null_mean'][i], 'null_sd': res['null_sd'][i], 'z': res['z'][i],
            'r_cv_blocked': cv['comp_blocked'][i] if i < r else np.nan,
            'r_cv_random': cv['comp_random'][i] if i < r else np.nan,
            'r_cv_folds': fold_str,
            'cond_x': res['cond_x'], 'cond_y': res['cond_y'],
        })
    rows.append({
        'component': 'subspace', 'r': res['sumcos2'], 'frac': res['sumcos2_frac'],
        'null_mean': res['sc2_null_mean'], 'null_sd': res['sc2_null_sd'], 'z': res['sc2_z'],
        # r_cv_blocked/random are the excess-over-permuted-null generalization signal;
        # raw held-out sum cos^2 and its baseline null are recorded in r_cv_folds.
        'r_cv_blocked': cv['sub_blocked'], 'r_cv_random': cv['sub_random'],
        'r_cv_folds': (f"blk_raw={cv['sub_blocked_raw']:.4f};blk_null={cv['sub_blocked_null']:.4f};"
                       f"rnd_raw={cv['sub_random_raw']:.4f};rnd_null={cv['sub_random_null']:.4f};"
                       f"p_emp={res['sc2_p_emp']:.2e}"),
        'cond_x': res['cond_x'], 'cond_y': res['cond_y'],
    })
    pd.DataFrame(rows).to_csv(out_path, sep='\t', index=False)
    print(f'  saved {out_path}')


def write_weights(res, A, B, boot_h, boot_m, human_vx, mouse_vx, out_human, out_mouse):
    """Unit-norm canonical weights per species (rows CCA1..CCAr), with bootstrap stability."""
    r = res['r']
    for side, vx, boot, out in [('human', human_vx, boot_h, out_human), ('mouse', mouse_vx, boot_m, out_mouse)]:
        rows = []
        for i in range(r):
            a, b = sign_fix(A[:, i].copy(), B[:, i].copy())
            w = a if side == 'human' else b
            w = w / np.linalg.norm(w)
            row = {'component': f'CCA{i+1}', 'canonical_r': res['spec'][i]}
            row.update({v: w[j] for j, v in enumerate(vx)})
            if i in boot:
                med, p5, frac = boot[i]
                row.update({'boot_cos_median': med, 'boot_cos_p5': p5, 'boot_frac_below': frac})
            else:
                row.update({'boot_cos_median': np.nan, 'boot_cos_p5': np.nan, 'boot_frac_below': np.nan})
            row['stable'] = (i in boot_h and i in boot_m
                             and boot_h[i][0] >= STABLE_THRESH and boot_m[i][0] >= STABLE_THRESH)
            rows.append(row)
        pd.DataFrame(rows).to_csv(out, sep='\t', index=False)
        print(f'  saved {out}')
