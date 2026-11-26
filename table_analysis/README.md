[English](README.md) | [简体中文](README.zh-CN.md)

# Table analysis methods

[Repository overview](../README.md) · [Main figure methods](../main/README.md)

These scripts work with externally prepared count matrices, annotations, model objects or trees. Where a table shares its method with a main figure, the entry below points to that figure's analysis code.

## Main Tables 1 and 2

| Table | Source of the reported values |
| --- | --- |
| Table 1. Primary tests for the host-associated viral sets | `run_host_models()` in the [Figure 4 analysis](../main/figure4_host_analysis.R): two-sided Fisher's exact tests for detection, and HC3 linear models of log2 viral relative abundance among patients with positive abundance, with BH correction in a family of eight tests. |
| Table 2. Evidence and taxonomic assignments for the two focal viral representatives | Existing annotations, MAG alignments and marker trees. The comparison of target-containing groups across markers and the nearest formal reference by patristic distance come from `marker_consistency()` and `tree_neighbourhood()` in the [Figure 5 analysis](../main/figure5_analysis.py). |

## Sample and catalogue summaries

| Script and functions | Inputs, calculation and interpretation |
| --- | --- |
| [Alignment quality control](alignment_qc.py) · `summarize()` | Pass aligned vOTU-by-sample arrays of mapped counts, RPKM and coverage breadth (a fraction), per-sample clean read totals and sample labels in column order. It totals counts before and after retaining entries with breadth ≥0.5, sums retained RPKM and checks raw mapped counts against clean read totals. A patient has an available viral composition only when retained counts and RPKM both have a positive total. |
| [Prevalence](prevalence.R) · `prevalence_filter()` | Supply the count matrix, minimum count for detection, prevalence fraction and scenario label. For each vOTU, the function counts samples meeting the detection threshold and records whether that count reaches `ceiling(number of input samples × prevalence fraction)`. The denominator is the number of columns supplied. It returns a decision and reason for each vOTU. |
| [Alpha diversity](diversity.R) · `alpha_values()`, `diversity_values()` | Supply relative abundances summing to one per sample and mapped counts with matching feature and sample names. The functions return Shannon diversity (natural logarithm), inverse Simpson diversity and the number of vOTUs with positive counts. `alpha_values()` requires positive total counts in each supplied sample. |
| [Predicted host assignment](host_assignment.py) · `assign()`, `number15()`, `text_bool()` | Supply catalogue vOTU IDs and iPHoP records with `Virus`, `Host genus` and `Confidence score`. Scores ≥90 qualify. A single genus at the highest qualifying score is assigned even when it appears more than once. Ties across genera yield `Ambiguous`, and absence of a qualifying record yields `Unassigned`. The two formatting helpers record numerical scores and Boolean flags. The assignments summarize computational predictions. |

## Models and correlations

| Code | Inputs, method and output |
| --- | --- |
| [Figure 2 analysis](../main/figure2_lifestyle_analysis.R) · `lifestyle_analysis()` | In R, source the script and call `lifestyle_analysis()` with breadth-filtered RPKM, named vOTU categories and patient metadata. It calculates within-patient percentages and uses `safe_lifestyle_fit()` for the age- and sex-adjusted severity contrast. HC3 intervals use residual degrees of freedom. The three categories form a BH family per catalogue scenario. `figure2_sensitivity()` assembles five contrasts, and `same_genus_coverage()` summarizes existing best MAG hits with >50% query coverage in the same genus. |
| [Figure 3 analysis](../main/figure3_paired_analysis.py) · `analyze()` | Five helpers support the correlations. `unit_ranks()` averages tied ranks and normalizes the rank vectors, `rho()` computes Spearman correlation, `permutation()` permutes bacterial patient labels for two-sided tests, `boot_ci()` resamples paired patients, and `bh()` takes the family size as an argument. For the full analysis, call `analyze()` with patient-aligned viral counts and RPKM, host assignments, MAG counts and genera, and both random seeds. Only genus-level correlations receive bootstrap intervals. |
| [Figure 4 analysis](../main/figure4_host_analysis.R) · `run_host_models()` | `prior_mag_score()` constructs the bacterial composite score from full-cohort MAG counts. `paired_ranks()` ranks both abundances within one genus across severity groups, retaining valid zeros. `run_host_models()` reads five TSVs and tests detection by two-sided Fisher exact test. It fits HC3 linear models of log2 viral fraction among patients with positive abundance. Estimable P values are adjusted in a family of eight tests. |
| [Collinearity diagnostics](ecology_diagnostics.R) · `vif_table()` | Supply a fitted model with an intercept and multiple predictor columns, plus a model label. For each non-intercept model-matrix column, the function regresses it on the others and returns its variance inflation factor, `1/(1 − R²)`. A factor represented by several columns has a value for each column. |

From the repository root, run `Rscript main/figure4_host_analysis.R EXTERNAL_WORKDIR` for the Figure 4 models. In a separate writable directory, supply `rpkm.tsv`, `counts.tsv`, `metadata.tsv`, `bacteria.tsv` and `hosts.tsv`. Results, set membership, sample metrics, diagnostics and fitted model objects are written there. Figures 2 and 4 use the [shared HC3 functions](../main/hc3.R) and the R package `sandwich`. The correlation implementation requires NumPy, pandas and SciPy.

## Tree and coordinate summaries

[Tree summaries](tree_summary.py) provides two functions. `neighbourhood()` takes an existing tree with branch lengths, a target tip, tip categories and the requested counts of neighbours and ICTV references. It returns the selected tips and a row for every tip giving its distance to the target and reason for display. Identifiers break distance ties for display only. The [Figure 5 analysis](../main/figure5_analysis.py) calls this function for its local tree view. `union_bases()` takes inclusive coordinate pairs, merges overlapping or adjacent intervals and returns the covered base count with the merged intervals. Tree inference and branch support estimation are upstream steps.
