# v1-multiome-continua

Code to reproduce the L2/3 archetype analyses of mouse primary visual cortex (V1). L2/3 neurons
vary along continua between archetypes A', B' and C'. The archetype genes are rebuilt from raw
counts (`2_development/01`, identical to Table S4) and are also used by results 3-5. Five results:

| # | Result | Folder |
|---|---|---|
| 1 | Mouse vs human L2/3 | `scripts/1_mouse_human` |
| 2 | Development of mouse L2/3 under normal rearing (NR, P6-P21) | `scripts/2_development` |
| 3 | Regulon enrichment of the L2/3 archetype genes | `scripts/3_regulon` |
| 4 | Effect of dark rearing (DR, P12-P21) | `scripts/4_dr` |
| 5 | EASI-FISH projection labeling along L2/3 depth | `scripts/5_easifish` |

## Quick start
```
conda env create -f environment.yml
conda activate v1-multiome-continua
pip install --no-deps -e .     # shared code src/l23/, importable by every script
# put the inputs in data/raw/ (see "Inputs")
bash run_all.sh
```
`run_all.sh` runs every script in number order, folder by folder (0_prep, then 1_mouse_human to
5_easifish); each folder reads only from itself and lower-numbered folders. It writes tables to `results/<folder>/`, figures to
`figures/<folder>/` and one log per script to `results/logs/`. The whole run takes about 3 minutes.
Each output file name starts with the number of the script that writes it. Within a folder, the
computing scripts come first and write tables; the figure scripts (`*figure*.py`, and those listed
with figures below) come after them and only read tables. Numerical notes on the pinned
environment (numpy, OpenBLAS, DESeq2) are in [REPRODUCIBILITY.md](REPRODUCIBILITY.md).

## Scripts and figures

### 0. Shared prep
| Script | What it does |
|---|---|
| `0_prep/01.extract_l23.py` | L2/3 neurons with raw counts: Yoo25 of every age (NR P6-P21, DR P12-P21) and Cheng22 P28 NR; every later step reads these extracts |
| `0_prep/02.extract_merfish_l23.py` | L2/3 MERFISH cells, P28 NR (for 1/06) and P21 NR (for 5/02); later steps read these extracts |

### 1. Mouse vs human L2/3
| Script | What it does | Figures |
|---|---|---|
| `01` | Mouse (Cheng22 P28 NR) and human (Jorstad23) L2/3: PCA + varimax | |
| `02` | Mouse L2/3 archetypes, markers and scores | |
| `03` | Mouse archetype scores on human cells | |
| `04` | Cross-species CCA of the varimax axes | |
| `05` | A'/B'/C' scores along CCA1/CCA2, human vs mouse: binned profiles and human-mouse r; per-cell CCA coordinates | |
| `06` | Conserved A vs C genes in situ, mouse P28 MERFISH vs human MERFISH (Fang22): relative L2/3 depth, profiles | |
| `07` | Figure of 05 | `07.l23_abc_normalized_compact.pdf` |
| `08` | Figures of 06 | `08.*.pdf` |
| `09` | Mouse archetype scores (B', C' − A') of human and mouse cells on the CCA1 x CCA2 embedding, mouse archetype triangle (from 05) | `09.l23_abc_scores_cca_embedding.pdf` |

### 2. Development (NR)
| Script | What it does | Figures |
|---|---|---|
| `01` | Per-cell A/B/C scores, archetype cells and archetype genes per age (NR and DR); checks the exact match to Table S4 | |
| `02` | Archetype cell groups and pseudobulk counts, P6-P21 | |
| `03` | Pseudobulk DESeq2: DR effect and its change with age (P12-P21) | |
| `04` | Archetype gene modules M1-M5 over age, NR vs DR | |
| `05` | Archetype chromatin (ATAC) region modules over age, NR vs DR | |
| `06` | Cells in P8 / P21 PC space of the global HVGs; variance per age along the P21 archetype A-C and B axes | |
| `07` | Number of archetype genes per age, NR vs DR | `07.*.pdf` |
| `08` | Figures of 04 (gene modules) and 05 (region modules) | `08.*.pdf` |
| `09` | Figures of 06 | `09.*.pdf` |

### 3. Regulon enrichment
| Script | What it does | Figures |
|---|---|---|
| `01` | Yoo25 L2/3 SCENIC+ regulon -> target table | |
| `02` | Gene universe: genes expressed in adult L2/3 (Cheng22 + Yoo25) and their expression stats | |
| `03` | Expression-stratified regulon x archetype enrichment; enriched TFs | |
| `04` | TF -> target network: edges and layout | |
| `05` | Human (Wang25) regulon -> target table | |
| `06` | AP-1/EGR regulon target activity over development, human (Wang25) and mouse (B' TFs) | |
| `07`, `08` | Dot plot of enriched regulons; archetype-gene overlap bars | `07.*.pdf`, `08.*.pdf` |
| `09` | Figure of 04 | `09.regulon_network_free.pdf` |
| `10` | Figure of 06, mouse vs human | `10.*.pdf` |

### 4. Dark rearing
Results 2's `03`, `04` and `05` also show NR vs DR.

| Script | What it does | Figures |
|---|---|---|
| `01` | Target activity of archetype-enriched regulons, NR vs DR | |
| `02` | Selected genes over age, NR vs DR, per archetype cell group | `02.*.pdf` |
| `03` | Figures of 01 | `03.*.pdf` |

### 5. EASI-FISH
| Script | What it does | Figures |
|---|---|---|
| `01` | Per-cell A'-C' axis, Yoo25 P21 NR | |
| `02` | Archetype scores along the A'-C' axis (snRNA) and depth (P21 MERFISH); EASI-FISH LM vs RL depth profiles (3 samples) | |
| `03` | Figure of 02 | `03.*.pdf` |
| `04` | EASI-FISH LM vs RL, 3 samples: depth profiles, x-y projection bias, per x-y bin mean-depth difference vs a label-shuffle null | |
| `05` | Figures of 04: cell maps, depth profiles, x-y maps, null histograms | `05.*.pdf` |

## Inputs
All inputs go in `data/raw/`, which is not tracked:

| File in `data/raw/` | Content | Used by | Source |
|---|---|---|---|
| `superdupermegaRNA_hasraw.h5ad` | Cheng22 + Yoo25 mouse V1 snRNA-seq, raw counts | 0_prep/01 | TBD |
| `jorstad23_human_WithinArea_L23IT.h5ad` | Human L2/3 IT (Jorstad et al. 2023) | 1/01, 1/03 | TBD |
| `wang25_human_multiome_gex_v1_l23.h5ad` | Human developing V1 L2/3, multiome GEX (Wang et al. 2025) | 3/06 | TBD |
| `regulons/yoo25_L2_3/L2_3_{direct,extended}_e_regulon_metadata.csv` | Yoo25 L2/3 SCENIC+ eRegulon metadata | 3/01 | TBD |
| `regulons/human_regulons_SuppTable13.xlsx` | Wang25 human eRegulons (Supp. Table 13) | 3/05 | TBD |
| `atac/pmat_l23concensus_{age}.h5ad` | Yoo25 L2/3 ATAC peak matrices, 11 conditions | 2/05 | TBD |
| `P28NRDR_v1glut_rna_merfish_240729.h5ad` | Mouse P28 NR/DR V1 MERFISH; only the NR sections' L2/3 labels are used | 0_prep/02 | TBD |
| `P28NR{a,b}_{ant,pos}_l2_v1_240723.h5ad` | Mouse P28 NR MERFISH sections (500 genes) | 0_prep/02 | TBD |
| `fang22_human_merfish/*.csv` | Human cortex MERFISH, 6 sections (Fang et al. 2022) | 1/06 | TBD |
| `P8P14P21NR_v1l23glut_rna_merfish_250411.h5ad` | Mouse P8/P14/P21 NR V1 L2/3 MERFISH; only the P21 NR sections are used | 0_prep/02 | TBD |
| `easifish/easifish_sample{1,2,3}_res.csv` | EASI-FISH retrograde labels (LM, RL), 3 samples | 5/02, 5/04, 5/05 | TBD |
| `Cheng22_Xie25_L23_ABC_genes_n286.csv` | Adult L2/3 archetype genes, the seed of the archetype genes | 2/01 | TBD |
| `gene_sets_l23{a,b,c}_250423.csv` | Module (M1-M5) of each archetype gene | 2/04 | TBD |
| `atac_modules/L23{A,B,C}_regions_M{1-5}.bed` | Archetype chromatin region modules | 2/05 | TBD |
| `human_mouse_orthologs.tsv` | Human-mouse ortholog pairs | 1/03, 1/04, 1/06 | TBD |
| `TableS4_L23_archetypal_genes.tsv` | Published Table S4; exact-match check only | 2/01 | TBD |

## Layout
```
v1-multiome-continua/
├── README.md
├── REPRODUCIBILITY.md      numerical notes on the pinned environment
├── run_all.sh              runs every script in dependency order
├── environment.yml         pinned conda env (see REPRODUCIBILITY.md)
├── data/raw/               inputs (not tracked; see "Inputs")
├── src/l23/                shared code, imported by the scripts
│   ├── constants.py        Yoo25 ages, NR vs DR test days, colours
│   ├── expr.py             log2(CP10k+1), percentile-scaled gene-set scores
│   ├── timecourse.py       NR vs DR OLS test and time-course panel
│   ├── archetype.py        archetype analysis (PCA + PCHA)
│   ├── utils.py            helpers of archetype.py: normalization, projection, PCHA
│   ├── crossspecies.py     cross-species varimax / CCA
│   ├── regulon.py          regulon enrichment, regulon tables, target activity
│   └── regulon_plot.py     shared encoding of the enrichment figures
├── scripts/
│   ├── 0_prep/             NN.*.py, numbered in run order
│   ├── 1_mouse_human/
│   ├── 2_development/
│   ├── 3_regulon/
│   ├── 4_dr/
│   └── 5_easifish/
├── results/<folder>/       tables (not tracked); results/logs/ = one log per script
├── figures/<folder>/       figures (not tracked)
└── plan/                   development plans
```
