# Converted from 31.multiome_basics_pcs_global_hvgs.ipynb (outputs dropped; shell "!" and IPython "%" lines commented out)
# Copied from scomics/scripts/revision/imported/analysis_multiome_py/31.multiome_basics_pcs_global_hvgs.py
# and run as is, with minimal changes to make it run here (all marked "# CHANGED"):
#   - inputs: L23_allmultiome_raw.h5ad -> results/0_prep/01.yoo25_l23.h5ad (same 38,539 cells);
#     nCount_RNA -> lib; var['features'] -> var_names
#   - scores_l23abc.csv -> results/0_prep/02.yoo25_l23_abc_scores.tsv (same values to float32
#     precision; score_* renamed to scores_*)
#   - gene annotation csv (L23-ABC-genes, loaded but unused) dropped
#   - scanpy, umap, scroutines (not installed here) dropped; GeneModules cell dropped;
#     powerplots.savefig_autodate -> fig.savefig
#   - Agg backend; plt.show() -> show(), which closes every open figure (not saved)
#   - notebook cell outputs of key tables wrapped in print()

# %%
import numpy as np
import pandas as pd
from scipy.stats import zscore
from scipy import sparse 
from scipy.spatial import ConvexHull
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns

from py_pcha import PCHA

# %%
import os
import numpy as np
import pandas as pd
import anndata 
import seaborn as sns
from scipy.stats import zscore
import matplotlib.pyplot as plt
import collections
from natsort import natsorted
import re

from scipy import stats
from scipy import sparse
from sklearn.decomposition import PCA
# from umap import UMAP  # CHANGED: not installed, unused

from matplotlib.colors import LinearSegmentedColormap

# CHANGED: scroutines not installed here
# from scroutines.config_plots import *
# from scroutines import powerplots # .config_plots import *
# from scroutines import pnmf
# from scroutines import basicu
# from scroutines.gene_modules import GeneModules
import matplotlib
matplotlib.use('Agg')

# %%
from scipy import stats
from matplotlib.ticker import MaxNLocator

def plot(x, y, aspect_equal=False, density=False, hue='type'):
    n = len(cases)
    fig, axs = plt.subplots(1,n,figsize=(4*n,4*1), sharex=True, sharey=True)
    for i, (ax, cond) in enumerate(zip(axs.flat, cases)):
        ax.set_title(cond)
        sns.scatterplot(data=res, 
                        x=x, y=y, 
                        c='lightgray',
                        s=1, edgecolor='none', 
                        legend=False,
                        ax=ax,
                        rasterized=True,
                       )
        if hue == 'type':
            sns.scatterplot(data=res[res['cond']==cond].sample(frac=1, replace=False),
                            x=x, y=y, 
                            hue='type',
                            hue_order=list(palette_types.keys()),
                            palette=palette_types,
                            s=3, edgecolor='none', 
                            legend=False,
                            ax=ax,
                            rasterized=True,
                           )
        else:
            sns.scatterplot(data=res[res['cond']==cond].sample(frac=1, replace=False),
                            x=x, y=y, 
                            hue='rep',
                            # hue_order=list(palette_types.keys()),
                            # palette=palette_types,
                            s=3, edgecolor='none', 
                            legend=False,
                            ax=ax,
                            rasterized=True,
                           )
            
        if density:
            sns.histplot(data=res[res['cond']==cond],
                            x=x, y=y, 
                            legend=False,
                            ax=ax,
                            rasterized=True,
                           )
        sns.despine(ax=ax)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
        if aspect_equal:
            ax.set_aspect('equal')
        if i > 0:
            ax.set_xlabel('')
            ax.set_ylabel('')
    return fig
    # plt.show()
    
def plot2(x, y, hue=None, aspect_equal=False, s=10, vmin=-2.5, vmax=2.5, vminp=None, vmaxp=None, cmap='coolwarm'):
    n = len(cases)
    fig, axs = plt.subplots(1,n,figsize=(4*n,4*1), sharex=True, sharey=True)
    fig.suptitle(hue, x=0, ha='left')
    for i, (ax, cond) in enumerate(zip(axs.flat, cases)):
        ax.set_title(cond)
        show = res[res['cond']==cond]
        
        if vminp is not None:
            vmin = np.percentile(show[hue], vminp)
        if vmaxp is not None:
            vmax = np.percentile(show[hue], vmaxp)
            
        if hue:
            # order = np.argsort(show[hue].values)
            ax.scatter(
                       # show[x].iloc[order], show[y].iloc[order], c=show[hue].iloc[order], 
                       show[x], show[y], c=show[hue], 
                       cmap=cmap,
                       vmin=vmin, vmax=vmax,
                       s=s, 
                       edgecolor='none', 
                       rasterized=True,
                      )
        else:
            r, p = stats.spearmanr(show[x], show[y])
            ax.scatter(show[x], show[y],  
                       s=s, 
                       edgecolor='none', 
                       rasterized=True,
                      )
            ax.set_title(f'{cond}\n r={r:.2f}')
        sns.despine(ax=ax)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
        if aspect_equal:
            ax.set_aspect('equal')
        if i == 0: 
            ax.set_xlabel(x)
            ax.set_ylabel(y)
        else:
            ax.set_xlabel('')
            ax.set_ylabel('')
        ax.grid(False)
    fig.tight_layout()
    
    return fig
    # plt.show()

# %%
import numpy as np
import pandas as pd
from scipy.stats import zscore
from scipy import sparse 
from scipy.spatial import ConvexHull
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns

from py_pcha import PCHA

from scipy.stats import gaussian_kde

# %%
def norm(x, depths):
    """
    Arguments: 
        x - cell by gene count matrix
        depths - sequencing depth per cell
        
    Output:
        xn - normalized count matrix

    This function takes raw counts as the input, and does the following steps sequencially.
         1. size normalization (CP10k) 
         2. log1p normalization (base 2 - log2(1+CP10k))
         3. zscore per gene  
    """

    xn = x/depths.reshape(-1,1)*1e4
    xn = np.log2(1+xn)
    xn = zscore(xn, axis=0)

    if np.any(np.isnan(xn)):
        print('Warning: the normalized matrix contains nan values. Check input.')
        xn = np.nan_to_num(xn, 0)

    return xn

# %%
def proj(x_norm, ndim, method='PCA'):
    """
    Arguments: 
        x_norm - normalized cell by gene feature matrix
        ndim   - number of dimensions

    Output:
        x_proj - a low-dimensional representation of `x_norm` 

    Here we only implemented PCA - a common projection method widely used, including by
    Adler et al. 2019 and Xie et al. 2024 for the Achetypal Analysis of scRNA-seq data.

    In principle, one can also choose to use other projection methods as needed.
    """

    if method == 'PCA':
        x_proj = PCA(n_components=ndim).fit_transform(x_norm)
    else:
        raise ValueError('methods other than PCA are not implemented...')

    return x_proj

# %%
def pcha(X, noc=3, delta=0, **kwargs):
    """
    """
    XC, S, C, SSE, varexpl = PCHA(X, noc=noc, delta=delta, **kwargs)
    XC = np.array(XC)
    XC = XC[:,np.argsort(XC[0])] # assign an order according to x-axis 
    return XC

# %%
def downsamp(x, which='cell', p=0.8, seed=None):
    """
    Arguments:
        x - cell by gene matrix
        which - downsample cells (rows) or genes (columns)
        p - fraction of cells/genes to keep - should be a value between ~ [0,1]
    """
    n0, n1 = x.shape
    
    rng = np.random.default_rng(seed=seed)
    
    if which in [0, 'cell', 'row']:
        return x[rng.random(n0)<p, :]
    elif which in [1, 'gene', 'col', 'column']:
        return x[:, rng.random(n1)<p]
    else:
        raise ValueError('choose from cell or gene')

# %%
def shuffle_rows_per_col(x, seed=None):
    """
    Arguments:
       x - cell by gene matrix
       seed - a random seed for reproducibility
    
    shuffles entries across rows (cells) independently for each col (gene)
    """
    rng = np.random.default_rng(seed=seed)
    x_shuff = rng.permuted(x, axis=0)
    return x_shuff

# %%
def plot_archetype(ax, aa, fmt='--o', color='k', **kwargs):
    """
    """
    ax.plot(aa[0].tolist()+[aa[0,0]], aa[1].tolist()+[aa[1,0]], fmt, color=color, **kwargs)

# %%
def get_t_ratio(xp, aa):
    """
    Arguments:
     xp -- projected matrix (cell by 2)
     aa -- inferred archetypes (2 by noc)
     note that this function only works for 2-dimensional space only
     
    Return: 
     t-ratio - ratio of areas (convex hull vs PCH)
     
    """
    assert xp.shape[1] == aa.shape[0] == 2
    
    ch_area  = ConvexHull(xp).volume
   
    x = aa[0]
    y = aa[1]
    pch_area = 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

    return ch_area/pch_area

# %%
class SingleCellArchetype():
    """
    """
    def __init__(self, x, depths, types):
        """
        Arguments: 
            x - cell by gene count matrix
            depths - sequencing depth per cell
            types  - cell type labels per cell
        
        Initiate the SingleCellArchetype object

        """
        
        # input
        self.x = x
        self.depths = depths
        self.types = types
        
        # cell type label
        types_idx, types_lbl = pd.factorize(types, sort=True)
        
        self.types_idx = types_idx
        self.types_lbl = types_lbl 
        
        # normalize
        self.xn = norm(self.x, self.depths)
        
        # feature matrix 
        self.xf = None 
        return 
        
    def setup_feature_matrix(self, method='data'):
        """
        """
        if method == 'data': 
            self.xf = self.xn
            print('use data')
            return  
        
        elif method == 'gshuff':
            # shuffle gene expression globally across all cells
            self.xf = shuffle_rows_per_col(self.xn)
            print('use shuffled data')
            return
            
        elif method == 'tshuff':
            # shuff each gene across cells independently - internally for each type A,B,C
            xn = self.xn
            xn_tshuff = xn.copy()
            
            types_lbl = self.types_lbl
            types_idx = self.types_idx
            for i in range(len(types_lbl)):
                xn_tshuff[types_idx==i] = shuffle_rows_per_col(xn[types_idx==i])
            self.xf = xn_tshuff
            print('use per-type shuffled data')
            return
        else:
            raise ValueError('choose from (data, gshuff, tshuff)')
    
    def proj_and_pcha(self, ndim, noc, **kwargs):
        """
        """
        xp = proj(self.xf, ndim)
        aa = pcha(xp.T, noc=noc, **kwargs)
        
        self.xp = xp
        self.aa = aa
        return (xp, aa)
        
    def downsamp_proj_pcha(self, ndim, noc, nrepeats=10, which='cell', p=0.8, **kwargs): 
        """
        """
        aa_dsamps = []
        for i in range(nrepeats):
            xn_dsamp = downsamp(self.xf, which=which, p=p)
            xp_dsamp = proj(xn_dsamp, ndim)
            aa_dsamp = pcha(xp_dsamp.T, noc=noc, **kwargs)
            aa_dsamps.append(aa_dsamp)
            
        return aa_dsamps
    
    def t_ratio_test(self, ndim, noc, nrepeats=10, **kwargs): 
        """
        this only work for 2-dimensional space for now
        """
        assert ndim == 2
        
        self.setup_feature_matrix(method='data')
        xp, aa = self.proj_and_pcha(ndim, noc)
        t_ratio = get_t_ratio(xp, aa)
        
        t_ratios_shuff = []
        for i in range(nrepeats):
            self.setup_feature_matrix(method='gshuff')
            xp_shuff, aa_shuff = self.proj_and_pcha(ndim, noc)
            t_ratio_shuff = get_t_ratio(xp_shuff, aa_shuff)
            t_ratios_shuff.append(t_ratio_shuff)
            
        pvalue = (np.sum(t_ratio > t_ratios_shuff)+1)/nrepeats
        
        return t_ratio, t_ratios_shuff, pvalue

# %%
# CHANGED: paths
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
INPUT_SCORES = os.path.join(REPO_ROOT, 'results', '0_prep', '02.yoo25_l23_abc_scores.tsv')
INPUT_ADATA  = os.path.join(REPO_ROOT, 'results', '0_prep', '01.yoo25_l23.h5ad')
OUTDIR_FIG   = os.path.join(REPO_ROOT, 'figures', '2_development', '07.multiome_basics_pcs_global_hvgs')
outdirfig = OUTDIR_FIG
os.makedirs(outdirfig, exist_ok=True)

def show():
    """CHANGED: stands in for plt.show(); closes every open figure."""
    plt.close('all')

# %% [markdown]
# # load gene annotation

# %%
# CHANGED: GeneModules (scroutines) not available
# gene_modules = GeneModules()
# anno, color, gene_styled = gene_modules.check_genes('Cdh13')
# print("\t".join(anno))
# print("\t".join(color))
# print("\t".join(gene_styled))

# %%
scores_abc = pd.read_csv(INPUT_SCORES, sep='\t', index_col=0) # CHANGED
scores_abc = scores_abc.rename(columns={'score_a': 'scores_a', 'score_b': 'scores_b', 'score_c': 'scores_c'})
scores_abc.index.name = None # reset_index() below then names the cell column 'index', as with the csv
print(scores_abc)

# %% [markdown]
# # load data

# %%
# f_anndata_in  = "../../data/v1_multiome/superdupermegaRNA_hasraw_multiome_l23.h5ad"# _multiome_l23.h5ad"#  L23_allmultiome_raw.h5ad"
f_anndata_in  = INPUT_ADATA # CHANGED
adata = anndata.read_h5ad(f_anndata_in)
print(adata)

# %%
sample_labels = ["-".join(cell.split(' ')[0].split('-')[2:]).replace('-2023', '') for cell in adata.obs.index]
time_labels = [s[:-1].replace('DR', '') for s in sample_labels]

adata.obs['n_counts'] = adata.obs['lib'] # CHANGED: was nCount_RNA
adata.obs['sample'] = sample_labels
adata.obs['time']   = time_labels

uniq_samples = natsorted(np.unique(sample_labels))
uniq_times = natsorted(np.unique(time_labels))

nr_samples = [s for s in uniq_samples if "DR" not in s]
dr_samples = [s for s in uniq_samples if "DR" in s]
print(uniq_times)
print(nr_samples)
print(dr_samples)

# adata.obs['sample'] = sample_labels

# %%
# select samples
adata.obs['cond'] = adata.obs['sample'].apply(lambda x: x[:-1]) # .unique()

# remove mitocondria genes
adata = adata[:,~adata.var_names.str.contains(r'^mt-')] # CHANGED: var['features'] -> var_names
adata = adata[:,~adata.var_names.str.contains(r'Xist')]

# select
print(adata.obs['sample'].unique(), adata.obs['cond'].unique())

# %%

print(adata.obs['Type'].unique())

# %%
# filter genes
cond = np.ravel((adata.X>0).sum(axis=0)) > 10 # expressed in more than 10 cells
adata = adata[:,cond]
genes = adata.var.index.values

# counts
x = sparse.csr_matrix(adata.X) # CHANGED: ensure sparse
cov = adata.obs['n_counts'].values

# CP10k
# xn = x/cov.reshape(x.shape[0], -1)*1e4
xn = (sparse.diags(1/cov).dot(x))*1e4

# log2(CP10k+1)
xln = xn.copy()
xln.data = np.log2(xln.data+1)

adata.layers[    'norm'] = np.array(xn.todense())
adata.layers[ 'lognorm'] = np.array(xln.todense())
adata.layers['zlognorm'] = zscore(np.array(xln.todense()), axis=0)

# %%
# select HVGs with mean and var
nbin = 20
qth = 0.3

# min
gm = np.ravel(xn.mean(axis=0))

# var
tmp = xn.copy()
tmp.data = np.power(tmp.data, 2)
gv = np.ravel(tmp.mean(axis=0))-gm**2

# cut 
lbl = pd.qcut(gm, nbin, labels=np.arange(nbin))
gres = pd.DataFrame()
gres['name'] = genes
gres['lbl'] = lbl
gres['mean'] = gm
gres['var'] = gv
gres['ratio']= gv/gm

# select
gres_sel = gres.groupby('lbl')['ratio'].nlargest(int(qth*(len(gm)/nbin))) #.reset_index()
gsel_idx = np.sort(gres_sel.index.get_level_values(1).values)
assert np.all(gsel_idx != -1)

# %%
adata_hvg = adata[:,gsel_idx]
print(adata_hvg.shape)
genes_hvg = adata_hvg.var.index.values

adata_hvg_early = adata_hvg[adata_hvg.obs['cond'].isin(['P8'])]
adata_hvg_later = adata_hvg[adata_hvg.obs['cond'].isin(['P21'])]

# %%
genes_hvg # to File
# np.savetxt('/u/home/f/f7xiesnm/v1_multiome/l23_alltime_hvgs_n4940.txt', genes_hvg, fmt='%s')

# %%
# learn on parts, project for everything
pca_early = PCA(n_components=5)
pca_early.fit(adata_hvg_early.layers['lognorm'][...]) #
pcs_early = pca_early.transform(adata_hvg.layers['lognorm'][...])

pcs_early[:,0] = -pcs_early[:,0] # flip PC1
# pcs_early[:,2] = -pcs_early[:,2] # flip PC3

pca_later = PCA(n_components=5)
pca_later.fit(adata_hvg_later.layers['lognorm'][...]) #
pcs_later = pca_later.transform(adata_hvg.layers['lognorm'][...])

pcs_later[:,0] = -pcs_later[:,0] # flip PC1
pcs_later[:,1] = -pcs_later[:,1] # flip PC2

# %%
pca_all = PCA(n_components=20)
pca_all.fit(adata_hvg.layers['lognorm'][...]) #
pcs_all = pca_all.transform(adata_hvg.layers['lognorm'][...])

pcs_all[:,0] = -pcs_all[:,0] # flip PC1
pcs_all[:,1] = -pcs_all[:,1] # flip PC1

# %%
pca_later_shuff = PCA(n_components=5)
pca_later_shuff.fit(shuffle_rows_per_col(adata_hvg_later.layers['lognorm'][...])) # random direction
pcs_later_shuff = pca_later_shuff.transform(adata_hvg.layers['lognorm'][...]) # original data

# %%
# %%time
# ucs_all = UMAP(n_components=2, n_neighbors=30).fit_transform(pcs_all) #

# %%
vt = pca_early.components_
topgenes_early = genes_hvg[np.flip(np.argsort(np.abs(vt), axis=1), axis=1)]
print(topgenes_early[:,:8])

vt = pca_later.components_
topgenes_later = genes_hvg[np.flip(np.argsort(np.abs(vt), axis=1), axis=1)]
print(topgenes_later[:,:8])

# %%
adata.obsm['pca_early'] = pcs_early
adata.obsm['pca_later'] = pcs_later
adata.obsm['pca_all']   = pcs_all

# %% [markdown]
# # set up for plotting

# %%
res0 = pd.DataFrame(adata.layers['zlognorm'][...], columns=genes)
res0['cond'] = adata.obs['cond'].values
res0['type'] = np.char.add('c', adata.obs['Type'].values.astype(str))
res0['samp'] = adata.obs['sample'].values
res0['rep']  = res0['samp'].apply(lambda x: x[-1])

res1 = pd.DataFrame(pcs_early, columns=np.char.add("pcs_early", ((1+np.arange(pcs_early.shape[1])).astype(str))))
res2 = pd.DataFrame(pcs_later, columns=np.char.add("pcs_later", ((1+np.arange(pcs_later.shape[1])).astype(str))))
res3 = pd.DataFrame(pcs_all  , columns=np.char.add("pcs_all",   ((1+np.arange(pcs_all.shape[1])).astype(str))))

res4 = pd.DataFrame(pcs_later_shuff, columns=np.char.add("pcs_later_shuff", ((1+np.arange(pcs_later_shuff.shape[1])).astype(str))))
res5 = scores_abc.reindex(adata.obs.index.values).reset_index().drop('cond', axis=1)
res = pd.concat([res0, res1, res2, res3, res4, res5], axis=1)

# res4 = pd.DataFrame(ucs_all, columns=np.char.add("ucs_all", ((1+np.arange(ucs_all.shape[1])).astype(str))))
# res = pd.concat([res0, res1, res2, res3, res4], axis=1)

# %%
allcolors = sns.color_palette('tab20b', 20)
allcolors

# %%
allcolors2 = sns.color_palette('tab10', 20)
allcolors2

# %%
palette = collections.OrderedDict({
     "P6": allcolors[1],
     "P8": allcolors[0],
    "P10": allcolors[4+2],
    "P12": allcolors[4+1],
    "P14": allcolors[4+0],
    
    "P17": allcolors[8+2],
    "P21": allcolors[8+0],
    
    "P12DR": allcolors[8+2],
    "P14DR": allcolors[8+0],
    "P17DR": allcolors[8+0],
    "P21DR": allcolors[8+0],
    
})
cases = np.array(list(palette.keys()))

cond_order_dict = {
    'P6':  0,
    'P8':  1,
    'P10': 2,
    'P12': 3,
    'P14': 4,
    'P17': 5,
    'P21': 6,
    
    'P12DR': 7,
    'P14DR': 8,
    'P17DR': 9,
    'P21DR': 10,
}
unq_conds = np.array(list(cond_order_dict.keys()))
adata.obs['cond_order'] = adata.obs['cond'].apply(lambda x: cond_order_dict[x])

palette_types = collections.OrderedDict({
    'L2/3_A': allcolors2[0],
    'L2/3_B': allcolors2[1],
    'L2/3_C': allcolors2[2],
})             

palette_types = {
    'c14': 'C0', 
    'c18': 'C1',
    'c16': 'C2', 
    
    'c13': 'C0', 
    'c15': 'C1', 
    'c17': 'C2',
}
type_order = [key for key, val in palette_types.items()]
type_order

# %%
palette_time = sns.cubehelix_palette(n_colors=7, start=.5, rot=-.5)
palette_time

# %% [markdown]
# # Plot A vs C genes aligning cells along early vs late PCs

# %%
# plot('pcs_all1', 'pcs_all2', aspect_equal=True)
# plt.show()
# plot('pcs_early1', 'pcs_early2', aspect_equal=True)
# plt.show()
# plot('pcs_later1', 'pcs_later2', aspect_equal=True)
# plt.show()

# %%


# %% [markdown]
# # Archetypes

# %%
todo_conds = [
    'P12DR', 'P14DR', 'P17DR', 'P21DR',
    'P6', 'P8', 'P10', 'P12', 'P14', 'P17', 'P21', 
]
todo_samps = [
    'P12DRa', 'P12DRb',
    'P14DRa', 'P14DRb',
    'P17DRa', 'P17DRb',
    'P21DRa', 'P21DRb',
    'P6a', 'P6b', 'P6c', 
    'P8a', 'P8b', 'P8c', 
    'P10a', 'P10b', 
    'P12a', 'P12b', 'P12c', 
    'P14a', 'P14b',
    'P17a', 'P17b', 
    'P21a', 'P21b', 
]

# %%
res['scores_c-a'] = res['scores_c'] - res['scores_a']

# %%
ndim = 2
noc = 3
delta = 0.2 # 0

aas = []

for cond in todo_conds:
    print(cond)
    xp = res[res['cond']==cond][['pcs_later1', 'pcs_later2']].values # pcs_later[:,:2]
    aa = pcha(xp.T, noc=noc, delta=delta)
    aas.append(aa)

# %%
# plot
n = len(todo_conds)
fig, axs = plt.subplots(1,n,figsize=(8*n,6), sharex=True, sharey=True)
for i, cond in enumerate(todo_conds):
    ax = axs[i]
    aa = aas[-1]
    xp = res[res['cond']==cond][['pcs_later1', 'pcs_later2']].values # pcs_later[:,:2]
    
    ax.scatter(xp[:,0], xp[:,1], s=2) #  c=types_colorvec, s=2)
    plot_archetype(ax, aa, fmt='-o', color='k', zorder=0)
    ax.set_title(cond)

    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_aspect('equal')
    sns.despine(ax=ax)
    ax.grid(False)
    
show()

# %%
xbins = np.linspace(-6,6,2*12+1)
ybins = np.linspace(-5,5,2*10+1)

xbins, ybins

# %%
from matplotlib.colors import LinearSegmentedColormap

colors_a = [(0.0, 'white'), (1.0, 'C0')]      
colors_b = [(0.0, 'white'), (1.0, 'C1')]      
colors_c = [(0.0, 'white'), (1.0, 'C2')]      
colors_n = [(0.0, 'white'), (1.0, 'lightgray')]      
cmap_a = LinearSegmentedColormap.from_list('cmap_a', colors_a)
cmap_b = LinearSegmentedColormap.from_list('cmap_b', colors_b)
cmap_c = LinearSegmentedColormap.from_list('cmap_c', colors_c)
cmap_n = LinearSegmentedColormap.from_list('cmap_n', colors_n)

colors_ac = [
    np.array(cmap_a(1.0)),
    0.5*np.array(cmap_a(1.0))+0.5*np.array(cmap_n(1.0)),
    np.array(cmap_n(1.0)),
    0.5*np.array(cmap_c(1.0))+0.5*np.array(cmap_n(1.0)),
    np.array(cmap_c(1.0)),
]
cmap_ac = LinearSegmentedColormap.from_list('cmap_ac', colors_ac)

colors_b2 = [
    'lightgray',
    cmap_b(0.2),
    cmap_b(1.0),
]
cmap_b2 = LinearSegmentedColormap.from_list('cmap_b2', colors_b2)

# %%

genes_viz = ['Nkain2', 'Nkain3', 'Meis2', 'Foxp1', 'Nfib', 'Rfx3', 'Sox5', 'Zbtb20', 'Tox', 'Tshz3']
for gn in genes_viz:
    fig = plot2('pcs_later1', 'pcs_later2', hue=gn, aspect_equal=True, s=20)
    axs = fig.get_axes()
    for ax in axs:
        # plot_archetype(ax, aa, fmt='--o', color='gray', zorder=2)
        plot_archetype(ax, aa, fmt='--', color='gray', zorder=2)
        # ax.scatter(aa[0,0], aa[1,0], color='C0', zorder=2)
        # ax.scatter(aa[0,1], aa[1,1], color='C1', zorder=2)
        # ax.scatter(aa[0,2], aa[1,2], color='C2', zorder=2)
    
    # output = os.path.join(outdirfig, f'pc12_{gn}.pdf')
    # powerplots.savefig_autodate(fig, output)
    show()
    
    break

# %%
genes_viz = ['scores_b','scores_c-a']
cmaps = [
    cmap_b2,
    cmap_ac,
]
for gn, cmap in zip(genes_viz, cmaps):
    
    vmin = np.percentile(res[gn],  5)
    vmax = np.percentile(res[gn], 95)
    
    fig = plot2('pcs_later1', 'pcs_later2', hue=gn, aspect_equal=True, s=20, vmin=vmin, vmax=vmax,cmap=cmap) #vminp=0, vmaxp=98)
    axs = fig.get_axes()
    for ax in axs:
        plot_archetype(ax, aa, fmt='--', color='gray', zorder=2)
        # ax.scatter(aa[0,0], aa[1,0], color='C0', zorder=2)
        # ax.scatter(aa[0,1], aa[1,1], color='C1', zorder=2)
        # ax.scatter(aa[0,2], aa[1,2], color='C2', zorder=2)
    
    output = os.path.join(outdirfig, f'pc12_{gn}.pdf')
    fig.savefig(output) # CHANGED: was powerplots.savefig_autodate
    show()
    
    # break

# %%
# genes_viz = ['scores_a', 'scores_b', 'scores_c', 'scores_c-a']
# for gn in genes_viz:
fig = plot2('scores_c-a', 'scores_b', aspect_equal=True, s=20) # , vmin=None, vmax=None, vminp=0, vmaxp=98)
show()
    
    # break

# %% [markdown]
# # variance along A-C axis vs B-axis

# %%
# define A-C (t1) and B (t2) axis

t1 = aa[:,2]-aa[:,0]
t1 = t1/np.sqrt(t1.dot(t1))
t2 = np.array([t1[1], -t1[0]])

print(np.sum(t1**2), t1.dot(t2))

# %%
fjg, ax = plt.subplots(figsize=(2,2))
plot_archetype(ax, aa, fmt='--', color='gray', zorder=2)
ax.plot([aa[0,0],aa[0,0]+t1[0]], [aa[1,0], aa[1,0]+t1[1]])
ax.plot([aa[0,0],aa[0,0]+t2[0]], [aa[1,0], aa[1,0]+t2[1]])
ax.set_aspect('equal')

# %%
res['pcs_later_t1'] = res[['pcs_later1', 'pcs_later2']].values.dot(t1)
res['pcs_later_t2'] = res[['pcs_later1', 'pcs_later2']].values.dot(t2)
res['pcs_later_shuff_t1'] = res[['pcs_later_shuff1', 'pcs_later_shuff2']].values.dot(t1)
res['pcs_later_shuff_t2'] = res[['pcs_later_shuff1', 'pcs_later_shuff2']].values.dot(t2)

# %%
var_mtx = res[['pcs_later_t1', 'pcs_later_t2', 
               'pcs_later_shuff_t1', 'pcs_later_shuff_t2', 
               'pcs_early1',
               'cond', 'samp']].groupby(['cond']).std(numeric_only=True) 
var_mtx = np.power(var_mtx, 2).reindex(todo_conds)
print(var_mtx)

# %%
var_mtx_rep = res[['pcs_later_t1', 'pcs_later_t2', 
                   'pcs_later_shuff_t1', 'pcs_later_shuff_t2', 
                   'pcs_early1',
                   'cond', 'samp']].groupby(['samp']).std(numeric_only=True)
var_mtx_rep = np.power(var_mtx_rep, 2).reindex(todo_samps)
print(var_mtx_rep)

# %%
times = np.array([6,8,10,12,14,17,21])
todo_conds_t = np.array([int(re.sub(r'[a-zA-Z]', '', a)) for a in todo_conds])
todo_samps_t = np.array([int(re.sub(r'[a-zA-Z]', '', a)) for a in todo_samps])
print(times)
print(todo_conds_t)
print(todo_samps_t)

# %%
scale1 = var_mtx['pcs_later_t1'].max()
scale2 = var_mtx['pcs_later_t2'].max()


fig, axs = plt.subplots(1,3,figsize=(3*3,1*3), sharex=True, sharey=True)
ax = axs[0]
ax.plot(todo_conds_t[4:], var_mtx['pcs_early1'].values[4:]/scale1,       color='C1')
ax.plot(todo_conds_t[:4], var_mtx['pcs_early1'].values[:4]/scale1,       color='k')

ax.plot(todo_samps_t[8:], var_mtx_rep['pcs_early1'].values[8:]/scale1,       'o', markersize=5, fillstyle='none', color='C1')
ax.plot(todo_samps_t[:8], var_mtx_rep['pcs_early1'].values[:8]/scale1,       's', markersize=5, fillstyle='none', color='k')
sns.despine(ax=ax)
ax.grid(False)
ax.set_xlabel('Postnatal day (P)')
ax.set_ylabel('norm. variance')
ax.set_title('P8 PC1 (early A-C axis)')

ax = axs[1]
ax.plot(todo_conds_t[4:], var_mtx['pcs_later_t1'].values[4:]/scale1,       color='C1')
ax.plot(todo_conds_t[:4], var_mtx['pcs_later_t1'].values[:4]/scale1,       color='k')
ax.plot(todo_conds_t[4:], var_mtx['pcs_later_shuff_t1'].values[4:]/scale1, color='gray', linestyle='--')
ax.plot(todo_conds_t[:4], var_mtx['pcs_later_shuff_t1'].values[:4]/scale1, color='lightgray', linestyle='--')

ax.plot(todo_samps_t[8:], var_mtx_rep['pcs_later_t1'].values[8:]/scale1,       'o', markersize=5, fillstyle='none', color='C1')
ax.plot(todo_samps_t[:8], var_mtx_rep['pcs_later_t1'].values[:8]/scale1,       's', markersize=5, fillstyle='none', color='k')
ax.plot(todo_samps_t[8:], var_mtx_rep['pcs_later_shuff_t1'].values[8:]/scale1, 'o', markersize=5, fillstyle='none', color='gray')
ax.plot(todo_samps_t[:8], var_mtx_rep['pcs_later_shuff_t1'].values[:8]/scale1, 's', markersize=5, fillstyle='none', color='lightgray')

sns.despine(ax=ax)
ax.grid(False)
ax.set_xlabel('Postnatal day (P)')
ax.set_ylabel('norm. variance')
ax.set_title('A-C axis')

ax = axs[2]
ax.plot(todo_conds_t[4:], var_mtx['pcs_later_t2'].values[4:]/scale2,       color='C1', label='NR')
ax.plot(todo_conds_t[:4], var_mtx['pcs_later_t2'].values[:4]/scale2,       color='k',  label='DR')
ax.plot(todo_conds_t[4:], var_mtx['pcs_later_shuff_t2'].values[4:]/scale2, color='gray', linestyle='--', label='shuffled NR')
ax.plot(todo_conds_t[:4], var_mtx['pcs_later_shuff_t2'].values[:4]/scale2, color='lightgray', linestyle='--', label='shuffled DR')

ax.plot(todo_samps_t[8:], var_mtx_rep['pcs_later_t2'].values[8:]/scale2,       'o', markersize=5, fillstyle='none', color='C1')
ax.plot(todo_samps_t[:8], var_mtx_rep['pcs_later_t2'].values[:8]/scale2,       's', markersize=5, fillstyle='none', color='k')
ax.plot(todo_samps_t[8:], var_mtx_rep['pcs_later_shuff_t2'].values[8:]/scale2, 'o', markersize=5, fillstyle='none', color='gray')
ax.plot(todo_samps_t[:8], var_mtx_rep['pcs_later_shuff_t2'].values[:8]/scale2, 's', markersize=5, fillstyle='none', color='lightgray')

sns.despine(ax=ax)
ax.grid(False)
ax.set_ylim(bottom=0)
ax.set_xticks(times)
ax.legend(bbox_to_anchor=(1,1))
ax.set_title('B axis')

output = os.path.join(outdirfig, "ACvsB_var_vs_time.pdf")
fig.savefig(output) # CHANGED: was powerplots.savefig_autodate
show()

# %% [markdown]
# # PC1 (for every sample)

# %%
res['pc1_this'] = 0

pc1_percond = {}

for cond in todo_conds:
    print(cond)
    adatasub = adata[adata.obs['cond']==cond]
    res_sub = res[res['cond']==cond]
    assert np.all(adatasub.obs.index.values == res_sub['index'].values)
    
    pc1_this = PCA(n_components=1).fit_transform(adatasub.layers['lognorm'][...])[:,0] #
    a_scores = res_sub['scores_a'].values
    
    r, _ = stats.spearmanr(pc1_this, a_scores)
    
    if r > 0:
        pc1_this = -pc1_this
        
    pc1_percond[cond]= pc1_this
    # break

# %%
def mean_bin(x, y, bins=None):
    df = pd.DataFrame() # 
    df['x'] = x
    df['y'] = y 
    
    if bins is None:
        binned, bins = pd.qcut(df['x'], 10, retbins=True)
        df['xbin'] = binned
    
    df_mean = df.groupby('xbin').mean()
    mid_bins = np.mean(np.vstack([bins[:-1], bins[1:]]).T, axis=1)
    
    df_mean['mid'] = mid_bins
    return df_mean

# %%
fig1, axs1 = plt.subplots(1,11, figsize=(3*11,3*1), sharex=True, sharey=True)
fig2, axs2 = plt.subplots(1,11, figsize=(3*11,3*1), sharex=True, sharey=True)
fig3, axs3 = plt.subplots(1,11, figsize=(3*11,3*1), sharex=True, sharey=True)

for i, cond in enumerate(todo_conds):
    print(cond)
    res_sub = res[res['cond']==cond]
    pc1_this = res_sub['pcs_later_t1']
    a_scores = res_sub['scores_a'].values
    b_scores = res_sub['scores_b'].values
    c_scores = res_sub['scores_c'].values
    
    ax = axs1[i]
    ax.set_title(cond)
    ax.scatter(pc1_this, a_scores, s=1, color='C0', rasterized=True)
    means = mean_bin(pc1_this, a_scores)
    ax.plot(means['x'], means['y'], color='k',)# , facecolor='none')
    
    ax = axs2[i]
    ax.scatter(pc1_this, b_scores, s=1, color='C1', rasterized=True)
    means = mean_bin(pc1_this, b_scores)
    ax.plot(means['x'], means['y'], color='k',)# , facecolor='none')
    
    ax = axs3[i]
    ax.scatter(pc1_this, c_scores, s=1, color='C2', rasterized=True)
    means = mean_bin(pc1_this, c_scores)
    ax.plot(means['x'], means['y'], color='k',)# , facecolor='none')
    
    # break
show()

# %%
all_vars = np.zeros((11,3))

for i, cond in enumerate(todo_conds):
    print(cond)
    pc1_this = pc1_percond[cond]
    res_sub = res[res['cond']==cond]
    coords_ac = res_sub['pcs_later_t1'].values
    coords_b  = res_sub['pcs_later_t2'].values
    
    var_pc1 = np.std(pc1_this)**2
    var_ac  = np.std(coords_ac)**2
    var_b   = np.std(coords_b)**2
    
    all_vars[i,0] = var_pc1
    all_vars[i,1] = var_ac
    all_vars[i,2] = var_b
                        
    # a_scores = res_sub['scores_a'].values
    # b_scores = res_sub['scores_b'].values
    # c_scores = res_sub['scores_c'].values

# %%
fig, axs = plt.subplots(1,2,figsize=(2*5,1*4))
ax = axs[0]
ax.plot(todo_conds_t[4:], all_vars[4:,1]/all_vars[4:,0], '-o', color='C1')
ax.plot(todo_conds_t[:4], all_vars[:4,1]/all_vars[:4,0], '-o', color='k')
ax.set_ylim(ymin=0)
ax.set_xticks(times)

ax = axs[1]
ax.plot(todo_conds_t[4:], all_vars[4:,2]/all_vars[4:,0], '-o', color='C1')
ax.plot(todo_conds_t[:4], all_vars[:4,2]/all_vars[:4,0], '-o', color='k')
ax.set_ylim(ymin=0)
ax.set_xticks(times)
print(pd.DataFrame(all_vars, index=todo_conds, columns=['var_pc1_this', 'var_ac', 'var_b']))
show()

# %% [markdown]
# # mean logFC of ABC archetype vs L2/3 mean in the space of HVGs.

# %%
num_archetypal_cells = 100
offset = 1
num_shuff = 10

abc_mean_l2fc = np.zeros((11,3))
abc_mean_l2fc_shuff = np.zeros((11,num_shuff))

for i, cond in enumerate(todo_conds):
    print(cond)
    adatasub = adata_hvg[adata_hvg.obs['cond']==cond]
    res_sub = res[res['cond']==cond]
    n_cells = len(adatasub)
    
    assert np.all(adatasub.obs.index.values == res_sub['index'].values)
    
    ranks_ac = res_sub['scores_c-a'].rank()
    ranks_b  = res_sub['scores_b'].rank()
    
    precond_a = ranks_ac <= num_archetypal_cells
    precond_c = ranks_ac > adatasub.shape[0] - num_archetypal_cells
    precond_b = ranks_b  > adatasub.shape[0] - num_archetypal_cells
    
    cond_a = np.all([ precond_a, ~precond_b, ~precond_c], axis=0)
    cond_b = np.all([~precond_a,  precond_b, ~precond_c], axis=0)
    cond_c = np.all([~precond_a, ~precond_b,  precond_c], axis=0)
    
        
    adatasub_a = adatasub[cond_a]
    adatasub_b = adatasub[cond_b]
    adatasub_c = adatasub[cond_c]
    
    mat_a = adatasub_a.layers['norm'][...]
    mat_b = adatasub_b.layers['norm'][...]
    mat_c = adatasub_c.layers['norm'][...]
    mat = adatasub.layers['norm'][...]
    
    for j, mat_j in enumerate([mat_a, mat_b, mat_c]):
        l2fc_j = np.mean(np.abs(np.log2(np.mean(mat_j, axis=0)+offset) - np.log2(np.mean(mat, axis=0)+offset))) # log2FC (CP10k as raw counts)
        abc_mean_l2fc[i,j] = l2fc_j
        
    for j in range(num_shuff):
        cond_n = np.random.choice(n_cells, size=num_archetypal_cells, replace=False)
        adatasub_n = adatasub[cond_n]
        mat_n = adatasub_n.layers['norm'][...]
        
        l2fc_j = np.mean(np.abs(np.log2(np.mean(mat_n, axis=0)+offset) - np.log2(np.mean(mat, axis=0)+offset))) # log2FC (CP10k as raw counts)
        abc_mean_l2fc_shuff[i,j] = l2fc_j

# %%
abc_mean_l2fc_shuff_mean = np.mean(abc_mean_l2fc_shuff, axis=1) # .shape
abc_mean_l2fc_shuff_err = 1.96*np.std(abc_mean_l2fc_shuff, axis=1)/np.sqrt(num_shuff) # .shape
print(pd.DataFrame(np.hstack([abc_mean_l2fc, abc_mean_l2fc_shuff_mean[:,None]]), index=todo_conds,
                   columns=['A', 'B', 'C', 'shuff_mean']))

# %%
fig, axs = plt.subplots(1,3,figsize=(3*4,1*4), sharey=True)
ax = axs[0]
ax.plot(todo_conds_t[4:], abc_mean_l2fc[4:,0])
ax.plot(todo_conds_t[:4], abc_mean_l2fc[:4,0])

ax.fill_between(todo_conds_t[4:], 
                 abc_mean_l2fc_shuff_mean[4:]-abc_mean_l2fc_shuff_err[4:],
                 abc_mean_l2fc_shuff_mean[4:]+abc_mean_l2fc_shuff_err[4:],
                )
ax.fill_between(todo_conds_t[:4], 
                 abc_mean_l2fc_shuff_mean[:4]-abc_mean_l2fc_shuff_err[:4],
                 abc_mean_l2fc_shuff_mean[:4]+abc_mean_l2fc_shuff_err[:4],
                )

ax = axs[1]
ax.plot(todo_conds_t[4:], abc_mean_l2fc[4:,1])
ax.plot(todo_conds_t[:4], abc_mean_l2fc[:4,1])

ax.fill_between(todo_conds_t[4:], 
                 abc_mean_l2fc_shuff_mean[4:]-abc_mean_l2fc_shuff_err[4:],
                 abc_mean_l2fc_shuff_mean[4:]+abc_mean_l2fc_shuff_err[4:],
                )
ax.fill_between(todo_conds_t[:4], 
                 abc_mean_l2fc_shuff_mean[:4]-abc_mean_l2fc_shuff_err[:4],
                 abc_mean_l2fc_shuff_mean[:4]+abc_mean_l2fc_shuff_err[:4],
                )

ax = axs[2]
ax.plot(todo_conds_t[4:], abc_mean_l2fc[4:,2])
ax.plot(todo_conds_t[:4], abc_mean_l2fc[:4,2])

ax.fill_between(todo_conds_t[4:], 
                 abc_mean_l2fc_shuff_mean[4:]-abc_mean_l2fc_shuff_err[4:],
                 abc_mean_l2fc_shuff_mean[4:]+abc_mean_l2fc_shuff_err[4:],
                )
ax.fill_between(todo_conds_t[:4], 
                 abc_mean_l2fc_shuff_mean[:4]-abc_mean_l2fc_shuff_err[:4],
                 abc_mean_l2fc_shuff_mean[:4]+abc_mean_l2fc_shuff_err[:4],
                )

show()

# %%

fig, axs = plt.subplots(1,3,figsize=(3*4,1*4), sharey=True)
ax = axs[0]
ax.plot(todo_conds_t[4:], abc_mean_l2fc[4:,0]-abc_mean_l2fc_shuff_mean[4:])
ax.plot(todo_conds_t[:4], abc_mean_l2fc[:4,0]-abc_mean_l2fc_shuff_mean[:4])


ax = axs[1]
ax.plot(todo_conds_t[4:], abc_mean_l2fc[4:,1]-abc_mean_l2fc_shuff_mean[4:])
ax.plot(todo_conds_t[:4], abc_mean_l2fc[:4,1]-abc_mean_l2fc_shuff_mean[:4])


ax = axs[2]
ax.plot(todo_conds_t[4:], abc_mean_l2fc[4:,2]-abc_mean_l2fc_shuff_mean[4:])
ax.plot(todo_conds_t[:4], abc_mean_l2fc[:4,2]-abc_mean_l2fc_shuff_mean[:4])


show()

# %% [markdown]
# # dev delay 
#
# - x P21/P10 FC (mean across all L2/3) for all HVGs
# - y DR / NR FC (mean across dev time)

# %%

