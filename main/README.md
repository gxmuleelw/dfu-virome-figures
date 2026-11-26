[English](README.md) | [简体中文](README.zh-CN.md)

# Main figure methods

The figures compare patients with Wagner grade 4 diabetic foot ulcers with those with grades 1–3. The analysis scripts calculate the quantities below. The drawing scripts use prepared summaries, figure settings or existing artwork. Commands assume the repository root as the working directory. Supply the study inputs separately and use a separate writable directory for outputs. [Repository overview](../README.md)

| Figure | Analysis or descriptive summary | Drawing |
| --- | --- | --- |
| **1b–f** | [figure1_catalogue_analysis.py](figure1_catalogue_analysis.py) · `summarize()` | [figure1_reflow.py](figure1_reflow.py) · `draw_panels()`, with `bars()` from [main_1.py](main_1.py) |
| **2a–e** | [figure2_lifestyle_analysis.R](figure2_lifestyle_analysis.R) · `lifestyle_analysis()`, `same_genus_coverage()`, `figure2_sensitivity()` | [figure2_render.py](figure2_render.py) · `main_figure()` |
| **3a–c** | [figure3_paired_analysis.py](figure3_paired_analysis.py) · `analyze()` | [figure3_render.py](figure3_render.py) · `host_overview()` |
| **4a–e** | [figure4_host_analysis.R](figure4_host_analysis.R) · `paired_ranks()`, `prior_mag_score()`, `run_host_models()` | [figure4_render.py](figure4_render.py) · `focal_hosts()` |
| **5a–c** | [figure5_analysis.py](figure5_analysis.py) · `marker_span()`, `marker_consistency()`, `tree_neighbourhood()` | [figure5_reflow.py](figure5_reflow.py) · command-line script |

The Figure 2 and 4 analyses share `fit_hc3()` and `hc3_result()` for regression inference. [hc3.R](hc3.R)

The renderers for Figures 2–4 each take a directory of prepared TSVs, a figure-specific JSON configuration and an output directory. For example, run `python main/figure2_render.py DATA_DIR CONFIG_JSON OUTPUT_DIR`. Use `main/figure3_render.py` or `main/figure4_render.py` with the same argument order to draw Figures 3 or 4. Each writes PDF, SVG and PNG.

## Figure 1. Viral catalogue and analysis populations

**1a** is the study-design diagram drawn in PowerPoint. **1b–d** count viral operational taxonomic units (vOTUs) by CheckV quality, viral cluster membership and deepest assigned taxonomic rank. **1e** counts the patients in whom each vOTU was detected, requiring a positive mapped count and at least 50% breadth. A patient with observed measurements and no confirmed detection remains in the detection denominator. Their relative viral composition is undefined. **1f** gives the number of patients in each Wagner group for detection, composition and ANCOM-BC2 analyses.

`python main/figure1_catalogue_analysis.py --help` lists the input tables and required fields. Pass the resulting JSON to `python main/figure1_reflow.py INPUT_JSON OUTPUT_DIR` to generate `Figure1b-f.pdf`.

## Figure 2. Viral lifestyle and ulcer severity

Load the analysis with `source("main/figure2_lifestyle_analysis.R")`. For **2a–c**, `lifestyle_analysis(rpkm, category, meta)` expresses each lifestyle category as a percentage of the patient's total retained viral RPKM. It estimates the age- and sex-adjusted difference in percentage points for Wagner grade 4 minus grades 1–3. HC3 intervals use residual degrees of freedom. The three category tests form one Benjamini–Hochberg (BH) family per catalogue scenario.

For **2d**, `same_genus_coverage()` counts vOTUs whose best hit to a metagenome-assembled genome (MAG) matches the predicted host genus and covers **strictly more than 50%** of the query. The denominator includes every vOTU in the lifestyle category. For **2e**, `figure2_sensitivity(rpkm, catalog, meta)` returns five contrasts for the predicted lytic percentage. These use the original catalogue, exclude vOTUs assigned to *Escherichia*, broaden the lytic label, adjust for the *E. coli* covariate, or exclude vOTUs with more host than viral genes. Catalogue exclusions change the denominator. The broader-label analysis retains the full denominator. The broader-label and bacterial-adjusted fits report raw P values without BH q values. Supply the bacterial covariate standardized in the full cohort.

## Figure 3. Viral and bacterial co-abundance

Import the analysis function with `from main.figure3_paired_analysis import analyze`. Supply breadth-filtered viral RPKM and counts, prepared host assignments, MAG counts and genus annotations, and both the permutation and bootstrap seeds. Read mapping and host predictions are external inputs.

**3a** correlates the fraction of viral abundance assigned to each eligible host genus with the bacterial fraction for that genus. **3b** repeats the comparison for individual vOTUs meeting the detection threshold. Both use Spearman ranks across patients with defined viral composition and two-sided tests that permute bacterial patient labels. Only the genus-level 95% intervals resample paired patients. **3c** counts associations with q < 0.05 for predicted-host and nonmatched pairs. Genus, member and full-matrix tests have separate BH families. Within each matrix, matched and nonmatched cells share one correction. Pass the excluded vOTU IDs as `excluded` to recalculate viral fractions and reselect eligible sets for the exclusion matrix, retaining the same patients and bacterial denominator. The function returns the genus, member, matrix and cross-pairing summary tables.

## Figure 4. Focal host genera

Load the analysis with `source("main/figure4_host_analysis.R")`. For **4a–b**, call `paired_ranks()` separately for *Escherichia* and *Collinsella*. It ranks viral and bacterial abundances across both Wagner groups together, retaining valid zeros. The line describes the ranks. The correlation, interval and q value come from the genus-level analysis in Figure 3. The bacterial abundance here sums MAGs assigned to the genus.

For **4c–e**, `run_host_models()` tests detection by Wagner group with two-sided Fisher exact tests (4c). Among patients with positive viral abundance it fits HC3 linear models of the log2 viral fraction to estimate the severity contrast (4d) and the association with the targeted bacterial score (4e). BH correction uses a family of eight tests and the available P values. Each positive-abundance model requires at least five patients with complete data in each Wagner group, at least twenty overall and at least ten residual degrees of freedom.

Construct the targeted score with `prior_mag_score()` before selecting positive patients. Express each of the three MAG counts as a percentage of all MAG counts in that patient, add half the smallest positive percentage across the three MAGs and all patients, then take log10. Standardize each MAG across the full cohort. The score is the *E. coli* z value minus the mean z value of the two *Collinsella* MAGs.

Run `Rscript main/figure4_host_analysis.R EXTERNAL_WORKDIR` with five prepared TSVs in that directory (`rpkm.tsv`, `counts.tsv`, `metadata.tsv`, `bacteria.tsv` and `hosts.tsv`). The command writes test results, membership, patient metrics, diagnostics and fitted model objects there.

## Figure 5. Evidence for the focal viral representative

`python -m main.figure5_analysis --help` lists three commands. **5a** uses `span` to summarize the coordinates of a contiguous marker range in an existing annotation table. **5b** uses `markers` to compare the first-ranked supported split containing the target across markers and count members by category when those target sides agree. **5c** uses `neighbourhood` to select tips for display by distance on a finished tree with branch lengths. Each command prints a JSON summary. The full-length MAG match, tree, branch supports and original artwork come from existing results.

To lay out the genome and tree artwork, run `python main/figure5_reflow.py SOURCE_SVG LABEL_JSON OUTPUT_DIR` with the compatible source SVG and labels. The script preserves the tree paths and writes SVG, PDF and PNG files.
