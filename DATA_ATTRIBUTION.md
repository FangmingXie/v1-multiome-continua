# Data

Every file in `data/input/`, the inputs of the analysis (folders 1-5). For each file: its content, the
scripts that read it, its source and its terms. `download_inputs.sh` installs these files from two
bundles; a copy of this file is in each bundle, as `ATTRIBUTION.md`.

- **Used by:** `2/01` = `scripts/2_development/01.*.py`.
- **Number prefix:** a file name with a number prefix was written by that 0_prep script from `data/raw/`
  (see [INPUT_PREP.md](INPUT_PREP.md)).
- **Source:** bracketed numbers refer to [References](#references).

## Core bundle

| File | Content | Used by | Source | Terms |
|---|---|---|---|---|
| `01.yoo25_l23.h5ad` | Yoo25 mouse V1 L2/3 snRNA-seq, every age (NR P6-P21, DR P12-P21), raw counts | 2/01, 2/02, 2/03, 3/01, 3/03, 4/01, 5/01 | this study [7] | CC BY 4.0 |
| `01.cheng22_p28nr_l23.h5ad` | Cheng22 mouse V1 L2/3 snRNA-seq, P28 NR, raw counts | 1/01, 1/02, 3/01 | our lab: Cheng et al. 2022 [1], GEO GSE190940 | GEO places no restrictions on use or distribution; cite [1] |
| `02.merfish_p28nr_l23.h5ad` | mouse P28 NR V1 MERFISH, L2/3 cells (500 genes) | 1/04 | our lab: Xie et al. 2025 [2], Zenodo 13916878 | CC BY 4.0. Changes: kept the P28 NR sections and the L2/3 cells |
| `02.merfish_p21nr_l23.h5ad` | mouse P21 NR V1 MERFISH, L2/3 cells (361 genes) | 5/01 | this study [7] | CC BY 4.0 |
| `03.yoo25_l23_regulon_targets.tsv` | Yoo25 L2/3 SCENIC+ regulon -> target table | 3/01, 3/02, 3/03, 4/01 | this study [7] | CC BY 4.0 |
| `04.human_wang25_regulon_targets.tsv` | Wang25 human eRegulon -> target table | 3/03 | Wang et al. 2025 [4], Supplementary Table 13b | cite [4]. Changes: direct and extended eRegulons merged (direct kept where both exist), one row per regulon-target pair |
| `atac/pmat_l23concensus_{age}.h5ad` | Yoo25 L2/3 ATAC peak matrices, 11 conditions | 2/06 | this study [7] | CC BY 4.0 |
| `atac_modules/L23{A,B,C}_regions_M{1-5}.bed` | archetype chromatin region modules | 2/06 | this study [7] | CC BY 4.0 |
| `easifish/easifish_sample{1,2,3}_res.csv` | EASI-FISH retrograde labels (LM, RL), 3 samples | 5/01, 5/03 | this study [7] | CC BY 4.0 |
| `Cheng22_Xie25_L23_ABC_genes_n286.csv` | adult L2/3 archetype genes, the seed of the archetype genes | 2/01 | our lab: Xie et al. 2025 [2], from the Cheng22 data [1] | CC BY 4.0 |
| `TableS4_L23_archetypal_genes.tsv` | L2/3 archetype genes per condition (Table S4 of this study); exact-match check only | 2/01 | this study [7] | CC BY 4.0 |
| `gene_sets_l23{a,b,c}_250423.csv` | module (M1-M5) of each archetype gene | 2/05 | this study [7] | CC BY 4.0 |
| `human_mouse_orthologs.tsv` | human-mouse ortholog pairs, every orthology type (the scripts keep 1-to-1 pairs) | 1/02, 1/03, 1/04 | Ensembl release 115 [6], BioMart: human genes with their mouse homologs (gene names) | no restrictions; cite [6] |

## External bundle

| File | Content | Used by | Source | Terms |
|---|---|---|---|---|
| `jorstad23_human_WithinArea_L23IT.h5ad` | human V1 L2/3 IT neurons, snRNA-seq | 1/01, 1/02 | Jorstad et al. 2023 [3], CZ CELLxGENE Discover | CC BY 4.0. Changes: kept the L2/3 IT cells |
| `wang25_human_multiome_gex_v1_l23.h5ad` | human developing V1 L2/3, snMultiome gene expression | 3/03 | Wang et al. 2025 [4], CZ CELLxGENE Discover | CC BY 4.0. Changes: kept the V1 L2/3 cells |
| `fang22_human_merfish/*.csv` | human cortex MERFISH, 6 sections (features, genes, matrix), unmodified | 1/04 | Fang et al. 2022 [5], Dryad | CC0; please cite [5] |

## Terms

- **Our lab's data** (this study [7] and [1], [2], unless stated otherwise above) are released in the
  bundles under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Please cite this study [7] when
  you use them.
- **CZ CELLxGENE Discover** datasets are CC BY 4.0: per its Data Submission Policy, users "can download,
  share, and use the data without restriction beyond providing attribution".
- **GEO:** NCBI "places no restrictions on the use or distribution of the GEO data".
- **Dryad** publishes every dataset under [CC0](https://creativecommons.org/publicdomain/zero/1.0/).
- **Ensembl** data are free to use, with no restrictions.

## References

1. Cheng S, Butrus S, et al. Vision-dependent specification of cell types and function in the developing
   cortex. *Cell* 185, 311-327 (2022). [doi:10.1016/j.cell.2021.12.022](https://doi.org/10.1016/j.cell.2021.12.022).
   Data: GEO [GSE190940](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190940).
2. Xie F, Jain S, Xu R, et al. Spatial profiling of the interplay between cell type- and vision-dependent
   transcriptomic programs in the visual cortex. *PNAS* 122 (2025).
   [doi:10.1073/pnas.2421022122](https://doi.org/10.1073/pnas.2421022122). Data: Zenodo
   [doi:10.5281/zenodo.13916878](https://doi.org/10.5281/zenodo.13916878).
3. Jorstad NL, et al. Transcriptomic cytoarchitecture reveals principles of human neocortex organization.
   *Science* 382, eadf6812 (2023). [doi:10.1126/science.adf6812](https://doi.org/10.1126/science.adf6812).
   Data: CZ CELLxGENE Discover, "Dissection: Primary visual cortex(V1)",
   [dataset c81e9267-8822-4eab-aad9-dcecccca3fa8](https://datasets.cellxgene.cziscience.com/c81e9267-8822-4eab-aad9-dcecccca3fa8.h5ad),
   [collection d17249d2-0e6e-4500-abb8-e6c93fa1ac6f](https://cellxgene.cziscience.com/collections/d17249d2-0e6e-4500-abb8-e6c93fa1ac6f).
4. Wang L, et al. Molecular and cellular dynamics of the developing human neocortex. *Nature* (2025).
   [doi:10.1038/s41586-024-08351-7](https://doi.org/10.1038/s41586-024-08351-7). Data: CZ CELLxGENE
   Discover, "snMultiome data of the developing human neocortex",
   [dataset a4310202-4dc8-4e1b-a96d-d9675f5b14d1](https://datasets.cellxgene.cziscience.com/a4310202-4dc8-4e1b-a96d-d9675f5b14d1.h5ad).
5. Fang R, et al. Conservation and divergence of cortical cell organization in human and mouse revealed by
   MERFISH. *Science* 377, 56-62 (2022). [doi:10.1126/science.abm1741](https://doi.org/10.1126/science.abm1741).
   Data: Dryad [doi:10.5061/dryad.x3ffbg7mw](https://doi.org/10.5061/dryad.x3ffbg7mw).
6. Dyer SC, et al. Ensembl 2025. *Nucleic Acids Research* 53, D948-D957 (2025).
   [doi:10.1093/nar/gkae1071](https://doi.org/10.1093/nar/gkae1071). Data: Ensembl release 115 (September
   2025), [BioMart](https://sep2025.archive.ensembl.org/biomart/martview): dataset Human genes, attributes
   "Gene name" and "Mouse gene name", unique rows with a mouse homolog; header renamed to `human_symbol`,
   `mouse_symbol`.
7. Yoo J, Xie F, Butrus S, et al. Establishing a continuum of cell types in the visual cortex. *bioRxiv*
   (2025). [doi:10.1101/2025.09.22.677893](https://doi.org/10.1101/2025.09.22.677893). This study; a
   revised version is in preparation.
