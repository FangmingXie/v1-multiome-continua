# Streamlining plan: v1-multiome-continua

Goal: every kept table stays the same, there are fewer and simpler scripts, each step either computes or plots, shared work runs once, and the figures you no longer need are gone.

## 1. Audit: is computation separate from plotting?

| Step | Today | Action |
|---|---|---|
| 0_prep/01-02 | compute only | ok |
| 1/01-05, 07-10, 12 | compute only | ok |
| 1/06 | CCA coords + binned profiles + permutation r **and** 3 figures | split: compute -> TSV; figure script (1 PDF) |
| 1/11 | MERFISH loading/flattening/profiles **and** 5 figures | split: compute -> TSVs; figure script (3 PDFs) |
| 1/13 | human stage stats **and** figure | becomes compute only (figure dropped) |
| 1/14 | mouse score **and** figure | mouse compute joins 13; 14 = figure only |
| 2/01 | Table S4 check **and** figure | move the S4 check into 0_prep/02 (where the genes are made); 2/01 = figure only |
| 2/02-04, 06 | compute | ok |
| 2/05, 07 | figure | ok (but the two are near copies; see 2h) |
| 3/01-04 | compute | ok |
| 3/05, 06 | figure | ok |
| 3/07 | star-rule re-selection + edges + layout **and** figures | split: compute (edges, layout positions) -> TSV; figure script |
| 4/01 | figure (light log-CPM transform inline) | ok |
| 4/02 | scores + OLS **and** 2 figures | split: compute; figure script |
| 5/01 | figure | delete (as requested) |
| 5/02 | compute | ok |
| 5/03 | snRNA/MERFISH/EASI-FISH binning **and** figure | split: compute; figure script |

## 2. Repeated computation: things to do once

a. **Yoo25 L2/3 raw counts are pulled from the 20+ GB `superdupermegaRNA_hasraw.h5ad` five times**
   (0_prep/01, 0_prep/02, 2/02, 1/14 and 4/02 via `regulon.load_expression`), each with its own copy of
   the Study/Subclass/Age filter and the AGES dict.
   -> 0_prep/01 extracts **all 11 ages** of Yoo25 L2/3 (38,539 cells) once to `results/0_prep/01.yoo25_l23.h5ad`
   (raw counts; obs Age, Sample, day, condition, lib computed in float64). Every later step reads that file.
   5/02-03 read it too, filtering to P21 NR. The P21-only extract goes away.

b. **"Percentile-scaled mean" score** (each gene divided by its p-th percentile, capped at 1, then averaged over
   genes) is written out 7 times: 1/11, 1/13, 2/02, 5/02, 5/03, `regulon.regulon_score`, plus the min-max
   2-98 variant in 1/03 and 1/04. -> `l23.expr.pctile_score()` and `l23.expr.minmax_score()`.

c. **log2(CP10k+1) from counts** is written out 9 times (0_prep/02, 2/02, 1/13, 5/02, 5/03,
   `regulon.load_expression`, `crossspecies.mouse_log2_cp10k`, 1/09 `load_mouse` (an exact copy of it),
   3/03 = `regulon.reconstruct_gene_universe`). -> `l23.expr.log2_cp10k()`.

d. **Adult L2/3 two-dataset expression** (cheng22 + yoo25, shared genes, mito out, log2 CP10k) is built
   three times: in 3/02 (counts), 3/03, and `regulon.reconstruct_gene_universe` in 3/04 (the last two are the same code).
   -> one loader `l23.expr.load_adult_l23()`. 3/02 also writes the universe gene stats (mean/var/detect),
   so 3/04 no longer loads either h5ad.

e. **PCA + varimax**: 1/01, 02, 07, 08 are four copies that differ only in species and N_HVG; 1/09 then
   reloads both h5ads and renormalizes to extend the 4000-HVG loadings.
   -> one varimax script that loads each species once and fits each N_HVG (and extends the
   loadings while the matrix is still in memory, if the hvg4000 branch is kept; see Q1).

f. **Archetype marker test** (conservative one-vs-each Wilcoxon) is duplicated in 1/03 and 3/03.
   -> `l23.archetype.archetype_markers()`, using vectorized `scipy.stats.ranksums(axis=0)`
   instead of a per-gene Python loop (3/03 currently takes about 30 s).

g. **The NR vs DR OLS model** `score ~ age + DR + DR*(day-16)` with BH, plus the D0/TEST_DAYS checks, is copied in
   2/04, 2/06 and 4/02. -> `l23.timecourse.dr_ols()`.

h. **The NR vs DR time-course panel** (lines + sample dots + bracket + stars) is copied in 2/05, 2/07, 4/01 and 4/02;
   2/05 and 2/07 are ~95% identical. -> `l23.timecourse.draw_nr_dr_panel()`, `stars()`, `draw_page()`.

i. **Star-rule TF selection** happens twice (3/04 `select()`, 3/07 `select_tfs()`) and gives the same table.
   -> the 3/07 compute step reads `04.enriched_regulon_selection.tsv`.

j. **SCENIC+ regulon dedup** (direct over extended, unique regulon x gene pairs) is in both 1/12 and 3/01 -> a small shared helper.

k. **Repeated constants** (AGES map, TEST_DAYS/D0, PCHA <-> A'/B'/C' relabel, ARCH_COLORS, COND_COLORS, STARS,
   Gate-A VX lists in 05/06/10/11) -> `l23/constants.py`.

New and changed modules in `src/l23/`: `expr.py` (new), `timecourse.py` (new), `constants.py` (new);
`archetype.py`, `regulon.py` and `crossspecies.py` slimmed so they use these.

## 3. Figure trimming (as requested) and the code it orphans

| Step | Keep | Drop | Also removed |
|---|---|---|---|
| 1/06 | `06.l23_abc_normalized_compact.pdf` | `..._ac_cca1_b_cca2_normalized{,_compact}.pdf` | the two-row figure code; stats now computed per ABC panel |
| 1/11 | `in_situ_all_sections`, `rel_depth_profiles`, `top_fit_all_sections` | `ac_scores_mouse_human`, `in_situ_mouse_human` | the A-vs-C score, `SELECTED`, `11.l23_selected_genes_ac_scores_cells.tsv`, the CCA1 per-gene z labels (only shown in the dropped figure) |
| 1/13-14 | `14.l23_b_regulon_targets_mouse_vs_human.pdf` | `13.wang25_..._over_stage.pdf` | 13 figure code (its TSV stays as input) |
| 3/07 | `07.regulon_network_free.pdf` | `*_free.png`, `*_alltargets.{pdf,png}` | `07.regulon_network_alltargets_edges.tsv`; this cuts 3/07 from about 80 s to about 40 s |
| 5 | `03.*.pdf` | `01.easifish_proj_depth.pdf` | script 5/01 |

**Orphan cascade (see Q1):** once 11 stops drawing CCA1 labels, nothing reads 1/07-10 (hvg4000 varimax,
loading extension, hvg4000 CCA) or 3/03 (adult PCHA markers, used only for 11's CCA1 orientation check).
The 3/02 Harmony coords are still needed (3/04 gets the gene-universe cells from them).

## 4. Proposed layout (recommended option: orphans removed)

Scripts are renumbered in run order, with compute before figure, following the existing 2_development 04/05 convention.

```
0_prep        01 extract_yoo25_l23 (all ages)      02 archetype_genes (+ Table S4 check)
3_regulon     01 regulon_targets   02 harmony_embed (+ universe stats)   03 regulon_enrichment
              04 regulon_network (edges + layout)  | figs: 05 dotplot  06 overlap_bar  07 network
1_mouse_human 01 varimax (mouse+human)  02 mouse_archetypes  03 scores_on_human  04 axis_cca
              05 abc_scores_along_cca  06 merfish_sections_depth  07 human_wang25_regulon_targets
              08 regulon_targets_over_development (human+mouse tables)
              | figs: 09 abc_along_cca  10 merfish_depth (3 PDFs)  11 b_regulon_targets_mouse_vs_human
2_development 01 cell_groups_pseudobulk  02 dr_slope_deseq2  03 gene_modules  04 atac_region_modules
              | figs: 05 archetype_gene_counts  06 gene_modules  07 atac_region_modules
4_dr          01 regulon_target_activity_nr_vs_dr | figs: 02 genes_time_course  03 regulon_target_activity (2 PDFs)
5_easifish    01 yoo25_p21nr_ac_axis  02 scores_and_proj_depth_bins | figs: 03 scores_and_proj_depth
```
`run_all.sh` keeps its folder order; README tables, docstrings and "Reads:" lines are updated to match.

## 5. Behavior preservation and verification

- Snapshot the current `results/` and `figures/` first (baseline from the pinned env).
- Refactor one folder at a time, and after each run that folder and diff every kept TSV against the baseline
  (exact match; the one exception is the ~1e-6 float32 marker noise of 3_regulon/03 already noted in the README).
  Render the kept PDFs with `pdftoppm` and pixel-diff them against the baseline.
- Points to watch: the extract must keep the hasraw cell order and float64 library sizes (0_prep/02, 2/02 and the
  regulon scores depend on them); vectorized `ranksums` must reproduce the per-gene p-values exactly.
- Finish with a full `run_all.sh` and a final diff of every table.

## 6. Issues found that this plan does NOT change (flagged only)

- 3/03 divides log2FC by ln 2 even though its matrix is already log2 (1/03 explicitly doesn't). This inflates its
  log2FC by 1.44x against the threshold. The point is moot if 3/03 is removed.
- The 1/14 mouse regulon score scales each gene over NR cells only; 4/02 scales over NR+DR cells. So the per-cell
  scores can't be shared without changing one figure. Kept as is; only the loader is shared.
- 0_prep/02 writes `abc_scores`, `abc_cells` and `abc_gene_stats.parquet`, which nothing reads; kept as provenance tables.

## 7. Decisions (2026-11-18)

- Orphans removed: 1_mouse_human 07-10 and 3_regulon/03 are deleted, along with their outputs.
- Follow-on: once 3/03 is gone, 3/04 uses 3/02 only for its list of cells (all L2/3 cells of cheng22 + yoo25, in
  file order), not for the Harmony embedding. 3/02 therefore becomes the adult L2/3 gene-universe step
  (shared genes, log2 CP10k, per-gene mean/var/detect) with no PCA or Harmony. The 3/03 enrichment table must stay identical.
- Numbering: scripts are renumbered in run order, compute before figures (layout in section 4).
- This plan lives in `plan/` and is not committed.

## 8. Outcome (2026-11-18)

- Done as planned. Verified against a fresh run of the original pipeline: all 45 kept tables are byte-identical
  (`1_mouse_human/08.l23_b_regulon_targets_nr_mouse.tsv` gains a `family` column; otherwise equal), and 16 of 17 kept
  figures are pixel-identical. The network figure differs only in its title: the empty "0 non-marker" count is gone.
- Bug caught by the pixel check and fixed: once 1/09 read p-values from a TSV, an exact float comparison made
  "p < 0.001" print as "p = 0.000999"; it now uses `np.isclose`.
- Full run: about 2.7 min (was about 5 min).
- `results/0_prep/01.yoo25_l23.h5ad` (all ages) is 1.0 GB, replacing the 0.18 GB P21-only extract.
- `harmonypy` was dropped from environment.yml (3_regulon/02 no longer runs Harmony); `src/l23/hvg.py` was deleted.
