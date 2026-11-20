# Data attribution

Sources and terms of the files in `data/input/`. These are the files that are shipped in the input bundles
(`download_inputs.sh`). A copy of this file is in each bundle, as `ATTRIBUTION.md`.

## This study and our lab's earlier publications (core bundle)

Unless stated otherwise below, our lab's data in the bundles are released under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Please cite this study (Yoo et al., in
preparation) when you use them.

| Files | Data | Source | Terms |
|---|---|---|---|
| `01.yoo25_l23.h5ad` | Yoo25 mouse V1 L2/3 snRNA-seq, every age (NR P6-P21, DR P12-P21), raw counts | this study (unpublished) | CC BY 4.0 |
| `atac/pmat_l23concensus_*.h5ad`, `atac_modules/*.bed` | Yoo25 L2/3 ATAC peak matrices; archetype chromatin region modules | this study (unpublished) | CC BY 4.0 |
| `03.yoo25_l23_regulon_targets.tsv` | Yoo25 L2/3 SCENIC+ regulon -> target table | this study (unpublished) | CC BY 4.0 |
| `02.merfish_p21nr_l23.h5ad` | mouse P21 NR V1 MERFISH, L2/3 cells | this study (unpublished) | CC BY 4.0 |
| `easifish/easifish_sample{1,2,3}_res.csv` | EASI-FISH retrograde labels (LM, RL) | this study (unpublished) | CC BY 4.0 |
| `TableS4_L23_archetypal_genes.tsv` | L2/3 archetype genes per condition (Table S4) | this study (unpublished) | CC BY 4.0 |
| `Cheng22_Xie25_L23_ABC_genes_n286.csv`, `gene_sets_l23{a,b,c}_250423.csv` | adult L2/3 archetype genes; their modules | our lab | CC BY 4.0 |
| `01.cheng22_p28nr_l23.h5ad` | Cheng22 mouse V1 L2/3 snRNA-seq, P28 NR, raw counts | Cheng S, Butrus S, et al. Vision-dependent specification of cell types and function in the developing cortex. *Cell* 185, 311-327 (2022). [doi:10.1016/j.cell.2021.12.022](https://doi.org/10.1016/j.cell.2021.12.022). Data: GEO [GSE190940](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190940) | Our lab's data. GEO places no restrictions on use or distribution. Cite Cheng et al. 2022 |
| `02.merfish_p28nr_l23.h5ad` | mouse P28 NR V1 MERFISH, L2/3 cells (500 genes) | Xie F, Jain S, Xu R, et al. Spatial profiling of the interplay between cell type- and vision-dependent transcriptomic programs in the visual cortex. *PNAS* 122 (2025). [doi:10.1073/pnas.2421022122](https://doi.org/10.1073/pnas.2421022122). Data: Zenodo [doi:10.5281/zenodo.13916878](https://doi.org/10.5281/zenodo.13916878) | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Changes: kept the P28 NR sections and the L2/3 cells |
| `human_mouse_orthologs.tsv` | 1-to-1 human-mouse ortholog pairs | upstream source not recorded | Cite the source once it is identified |

## External datasets (external bundle)

| Files | Data | Source | Terms |
|---|---|---|---|
| `jorstad23_human_WithinArea_L23IT.h5ad` | human V1 L2/3 IT neurons, snRNA-seq | Jorstad NL, et al. Transcriptomic cytoarchitecture reveals principles of human neocortex organization. *Science* 382, eadf6812 (2023). [doi:10.1126/science.adf6812](https://doi.org/10.1126/science.adf6812). Dataset "Dissection: Primary visual cortex(V1)", CZ CELLxGENE Discover, [dataset c81e9267-8822-4eab-aad9-dcecccca3fa8](https://datasets.cellxgene.cziscience.com/c81e9267-8822-4eab-aad9-dcecccca3fa8.h5ad), [collection d17249d2-0e6e-4500-abb8-e6c93fa1ac6f](https://cellxgene.cziscience.com/collections/d17249d2-0e6e-4500-abb8-e6c93fa1ac6f) | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Changes: kept the L2/3 IT cells |
| `wang25_human_multiome_gex_v1_l23.h5ad` | human developing V1 L2/3, snMultiome gene expression | Wang L, et al. Molecular and cellular dynamics of the developing human neocortex. *Nature* (2025). [doi:10.1038/s41586-024-08351-7](https://doi.org/10.1038/s41586-024-08351-7). Dataset "snMultiome data of the developing human neocortex", CZ CELLxGENE Discover, [dataset a4310202-4dc8-4e1b-a96d-d9675f5b14d1](https://datasets.cellxgene.cziscience.com/a4310202-4dc8-4e1b-a96d-d9675f5b14d1.h5ad) | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Changes: kept the V1 L2/3 cells |
| `fang22_human_merfish/*.csv` | human cortex MERFISH, 6 sections (features, genes, matrix), unmodified | Fang R, et al. Conservation and divergence of cortical cell organization in human and mouse revealed by MERFISH. *Science* 377, 56-62 (2022). [doi:10.1126/science.abm1741](https://doi.org/10.1126/science.abm1741). Data: Dryad [doi:10.5061/dryad.x3ffbg7mw](https://doi.org/10.5061/dryad.x3ffbg7mw) | [CC0](https://creativecommons.org/publicdomain/zero/1.0/) (Dryad). Please cite Fang et al. 2022 |
| `04.human_wang25_regulon_targets.tsv` (core bundle) | Wang25 human eRegulon -> target table | derived by 0_prep/04 from Supplementary Table 13b of Wang et al. 2025 (above) | Please cite Wang et al. 2025. Changes: direct and extended eRegulons merged (direct kept where both exist), one row per regulon-target pair |

CC BY 4.0 data on CZ CELLxGENE Discover: per its Data Submission Policy, users "can download, share, and use
the data without restriction beyond providing attribution". GEO: NCBI "places no restrictions on the use or
distribution of the GEO data".
