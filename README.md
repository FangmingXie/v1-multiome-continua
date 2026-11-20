# v1-multiome-continua

Code to reproduce the L2/3 archetype analyses of mouse primary visual cortex (V1). L2/3 neurons
vary along continua between archetypes A, B and C. Five results:

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
bash download_inputs.sh        # data/input/ from the input bundles (see "Data")
bash run_all.sh                # the analysis, folders 1-5
```
`download_inputs.sh` downloads the input bundles from Google Drive with `curl` (about 1.9 GB), checks
their sha256 and unpacks them into `data/input/`.

`run_all.sh` runs the analysis scripts in number order, folder by folder (1_mouse_human to 5_easifish),
in about 4 minutes. Each folder reads only from itself, lower-numbered folders and `data/input/`.
Scripts write tables to `results/<folder>/`, figures to `figures/<folder>/` and one log each to
`results/logs/`; each output file name starts with the number of the script that writes it. Within a
folder, the computing scripts come first and write tables; the figure scripts come after them and only
read tables.

Numerical notes on the pinned environment (numpy, OpenBLAS, DESeq2) are in
[REPRODUCIBILITY.md](REPRODUCIBILITY.md).

## Data
`download_inputs.sh` installs `data/input/` from two archives (version 2026-11-20) on Google Drive:

| Bundle | Content | Size (packed) | Needed by |
|---|---|---|---|
| [core](https://drive.google.com/file/d/10kRY471qCVIJj3j6GRAzAmrYdjXbCWNS/view?usp=sharing) | every input except the three external human datasets | 1.7 GB (741 MB) | folders 2-5 |
| [external](https://drive.google.com/file/d/1yhkWYA0xS0GC7xOS-CQS1K3ZN45ALzog/view?usp=sharing) | Jorstad23, Wang25 and Fang22 human datasets | 4 GB (1.1 GB) | folder 1, 3_regulon/03 and its figure 07 |

The script checks each archive against the sha256 recorded in `download_inputs.sh`, and that every file
listed in the archive's `MANIFEST_<bundle>.tsv` (path, size and sha256 of each file) is unpacked with its
size. `--core-only` skips the external bundle; `--from <dir>` installs archives already downloaded (e.g.
from the links above in a browser).

Each file's content, the scripts that read it, its source and licence: [DATA_ATTRIBUTION.md](DATA_ATTRIBUTION.md).
Rebuilding `data/input/` from `data/raw/` (0_prep) and publishing new bundles (maintainers):
[INPUT_PREP.md](INPUT_PREP.md).

## Scripts and figures
Each figure script's row links the PDFs it draws, so this section doubles as the figure report.
`run_all.sh` ends with `check_readme_figures.py`, which warns if a PDF in `figures/` is not linked here
or a link points to a missing file.

### 1. Mouse vs human L2/3
| Script | What it does | Figures |
|---|---|---|
| `01` | PCA + varimax of L2/3 neurons, mouse (Cheng22, P28 NR) and human (Jorstad23) | |
| `02` | Mouse archetypes (PCHA on the varimax axes), their marker genes, and archetype scores of mouse and human cells (via orthologs) | |
| `03` | CCA aligning the mouse and human varimax axes; archetype scores along CCA1 and CCA2 in both species, with the human-mouse correlation | |
| `04` | L2/3 depth profiles of conserved archetype genes, mouse P28 MERFISH vs human MERFISH (Fang22) | |
| `05` | Figures of 03: score profiles along CCA1 and CCA2; cells on the CCA embedding | [05.archetype_scores_along_cca.pdf](figures/1_mouse_human/05.archetype_scores_along_cca.pdf)<br>[05.archetype_scores_cca_embedding.pdf](figures/1_mouse_human/05.archetype_scores_cca_embedding.pdf) |
| `06` | Figures of 04: depth profiles and in situ maps | [06.merfish_gene_in_situ_maps.pdf](figures/1_mouse_human/06.merfish_gene_in_situ_maps.pdf)<br>[06.merfish_gene_depth_profiles.pdf](figures/1_mouse_human/06.merfish_gene_depth_profiles.pdf)<br>[06.merfish_l23_top_fit.pdf](figures/1_mouse_human/06.merfish_l23_top_fit.pdf) |

### 2. Development (NR)
| Script | What it does | Figures |
|---|---|---|
| `01` | Archetype scores per cell and archetype genes per age, NR and DR (reproduces Table S4) | |
| `02` | PCA of the L2/3 cells of every age (fit at P8 and P21); spread of each age along the A-C and B axes | |
| `03` | Archetype cell groups (top 20% of cells per age by each archetype score) and their pseudobulk counts per sample | |
| `04` | Pseudobulk DESeq2: effect of dark rearing and its change with age (P12-P21) | |
| `05` | Archetype gene modules (M1-M5) over age, NR vs DR | |
| `06` | Archetype chromatin (ATAC) region modules over age, NR vs DR | |
| `07` | Figures of 01 and 02: archetype gene counts per age; cells in P21 PC space colored by score; spread along the axes over age | [07.archetype_gene_counts.pdf](figures/2_development/07.archetype_gene_counts.pdf)<br>[07.axis_variance_over_age.pdf](figures/2_development/07.axis_variance_over_age.pdf)<br>[07.pc_space_archetype_scores.pdf](figures/2_development/07.pc_space_archetype_scores.pdf) |
| `08` | Figures of 05 and 06 | [08.atac_modules_by_cell_group.pdf](figures/2_development/08.atac_modules_by_cell_group.pdf)<br>[08.atac_modules_by_module.pdf](figures/2_development/08.atac_modules_by_module.pdf)<br>[08.gene_modules_by_cell_group.pdf](figures/2_development/08.gene_modules_by_cell_group.pdf)<br>[08.gene_modules_by_module.pdf](figures/2_development/08.gene_modules_by_module.pdf) |

### 3. Regulon enrichment
| Script | What it does | Figures |
|---|---|---|
| `01` | Which Yoo25 L2/3 regulons (TF -> targets) are enriched in each archetype's genes, against an expression-matched null; the enriched TFs | |
| `02` | Network of the enriched TFs and their archetype-gene targets | |
| `03` | Target activity of the B-enriched AP-1 and EGR regulons over development, mouse (Yoo25) vs human (Wang25) | |
| `04`, `05` | Figures of 01: enrichment dot plot; archetype-gene overlap bars | [04.regulon_enrichment_dotplot.pdf](figures/3_regulon/04.regulon_enrichment_dotplot.pdf)<br>[05.regulon_archetype_gene_overlap.pdf](figures/3_regulon/05.regulon_archetype_gene_overlap.pdf) |
| `06` | Figure of 02: the network | [06.regulon_network.pdf](figures/3_regulon/06.regulon_network.pdf) |
| `07` | Figure of 03: mouse vs human | [07.b_regulon_activity_mouse_vs_human.pdf](figures/3_regulon/07.b_regulon_activity_mouse_vs_human.pdf) |

### 4. Dark rearing
Folder 2's `04`-`06` also compare NR and DR, for genes and modules.

| Script | What it does | Figures |
|---|---|---|
| `01` | Target activity of the archetype-enriched regulons (3_regulon/01) over age, NR vs DR | |
| `02` | Selected genes over age, NR vs DR, per archetype cell group | [02.selected_genes_nr_vs_dr.pdf](figures/4_dr/02.selected_genes_nr_vs_dr.pdf) |
| `03` | Figures of 01 | [03.regulon_activity_nr_vs_dr.pdf](figures/4_dr/03.regulon_activity_nr_vs_dr.pdf)<br>[03.regulon_activity_p21_dr_change.pdf](figures/4_dr/03.regulon_activity_p21_dr_change.pdf) |

### 5. EASI-FISH
| Script | What it does | Figures |
|---|---|---|
| `01` | Archetype scores along the A-C axis (P21 snRNA) and along L2/3 depth (P21 MERFISH); EASI-FISH LM- vs RL-projecting cells: depth profiles and x-y distribution, against a label-shuffle null | |
| `02` | Figure of 01: archetype scores next to the LM and RL depth profiles | [02.archetype_scores_and_projection_depth.pdf](figures/5_easifish/02.archetype_scores_and_projection_depth.pdf) |
| `03` | Figures of 01: EASI-FISH cell maps, depth profiles and x-y analyses | [03.easifish_cell_maps.pdf](figures/5_easifish/03.easifish_cell_maps.pdf)<br>[03.easifish_depth_profiles.pdf](figures/5_easifish/03.easifish_depth_profiles.pdf)<br>[03.easifish_xy_depth_difference_null.pdf](figures/5_easifish/03.easifish_xy_depth_difference_null.pdf)<br>[03.easifish_xy_projection_bias.pdf](figures/5_easifish/03.easifish_xy_projection_bias.pdf) |

## Layout
```
v1-multiome-continua/
├── README.md
├── REPRODUCIBILITY.md      numerical notes on the pinned environment
├── DATA_ATTRIBUTION.md     every input file: content, scripts that read it, source, licence
├── INPUT_PREP.md           (maintainers) 0_prep, data/raw/, publishing bundles
├── download_inputs.sh      installs data/input/ from the input bundles
├── make_input_bundle.sh    (maintainers) packs data/input/ into the input bundles
├── run_input_prep.sh       (maintainers) runs 0_prep (data/raw/ -> data/input/)
├── run_all.sh              runs the analysis, folders 1-5, in dependency order
├── check_readme_figures.py checks that README.md links every figure (run by run_all.sh)
├── environment.yml         pinned conda env (see REPRODUCIBILITY.md)
├── data/                   not tracked
│   ├── input/              inputs of scripts 1-5 (DATA_ATTRIBUTION.md)
│   └── raw/                (maintainers) sources of 0_prep (INPUT_PREP.md)
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
│   ├── 0_prep/             (maintainers) data/raw/ -> data/input/ (INPUT_PREP.md)
│   ├── 1_mouse_human/
│   ├── 2_development/
│   ├── 3_regulon/
│   ├── 4_dr/
│   └── 5_easifish/
├── results/<folder>/       tables (not tracked); results/logs/ = one log per script
├── figures/<folder>/       figures (tracked)
└── plan/                   development plans (not tracked)
```
