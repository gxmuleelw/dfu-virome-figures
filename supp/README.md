[English](README.md) | [简体中文](README.zh-CN.md)

# Supplementary figures S1–S14

[Back to the repository overview](../README.md)

Each figure module exposes `draw(...)`. Supply the tables, estimates or tree objects listed below. Python returns a Matplotlib figure, which you can save with `fig.savefig("output.pdf")`. In R, source the selected script, call `draw(...)` and save the returned plot with `ggplot2::ggsave("output.pdf", plot = fig)`. Choose the output filename when saving.

| Figure | What it shows | Module and required inputs |
| --- | --- | --- |
| S1 | Catalogue mapping, detection, sequence lengths and CheckV gene counts | [S01.py](S01.py) takes sample-level QC, catalogue lengths and viral/host gene counts, and the annotation text, target coordinates and label position. |
| S2 | Effect of omitting each predicted-lytic member, and the 12 largest contributions | [S02.py](S02.py) takes the leave-one-vOTU-out estimates and intervals, member contribution table and original model coefficient. |
| S3 | Lifestyle differences across sensitivity analyses | [S03.py](S03.py) takes estimates and intervals by scenario, lifestyle category and covariate adjustment, with P and q values and four ordered scenario IDs and labels. |
| S4 | Shannon diversity, Bray–Curtis and robust Aitchison ordinations, and variation associated with severity | [S04.py](S04.py) takes sample Shannon values, the adjusted model summary row and group labels. For each distance, supply a tuple of ordination coordinates, PERMANOVA results, dispersion test results and ellipse geometry. |
| S5 | Differential abundance of individual vOTUs and gene counts for exploratory candidates (q < 0.25) | [S05.py](S05.py) takes all differential results, the selected candidate table with intervals, predicted host genus and CheckV gene counts, and the model sample size. |
| S6 | Associations between clinical phenotypes and virome features, with intervals for the ten smallest raw P values | [S06.R](S06.R) takes phenotype–feature estimates, 95% confidence limits, raw P, BH q and effective sample size. |
| S7 | Virus–bacterium co-abundance matrix with predicted-host matches | [S07.py](S07.py) takes one completed test row per bacterium–vOTU pair (Spearman ρ, matrix-wide BH q and host-match annotation) and the column indices to label. |
| S8 | *Escherichia* severity contrasts across scenarios and sensitivity of host co-abundance to member omission | [S08.py](S08.py) takes scenario estimates, correlations and intervals, the change in ρ after each member omission, four scenario keys and labels, and the IDs to highlight. |
| S9 | Summary of marker trees and node support | [S09.py](S09.py) takes one summary row per marker with the job identifier, tip count, trimmed amino-acid alignment length, selected model, SH-aLRT and UFBoot support. |
| S10 | Six full marker trees | [S10_S12_S13.py](S10_S12_S13.py) draws one supplied tree per page. Use `view="full"` and leave `chosen=None` to show all tips. |
| S11 | Gene organisation along viral representatives | [S11.R](S11.R) takes gene coordinates, orientations and feature classes, backbone lengths, label coordinates, panel titles and a colour palette. Supply the overlay tables and text for lysis intervals and notes as applicable. |
| S12 | Local neighbourhood on the terminase large subunit (TerL) tree | [S10_S12_S13.py](S10_S12_S13.py) takes the existing tree and a set of tip names in `chosen`. Use `view="local"` to display the selected neighbourhood and summarise omitted branches. |
| S13 | Three full marker trees | [S10_S12_S13.py](S10_S12_S13.py) draws each supplied tree with `view="full"` and `chosen=None`. |
| S14 | Abundance of viral carriers of resistance, auxiliary metabolic and virulence genes, and contrasts by auxiliary metabolic gene category | [S14.py](S14.py) takes the functional model table, participant-level carrier abundances with severity group and panel codes, and group labels. |

Python requires Matplotlib, NumPy and pandas. Make the `supp/` directory available on the Python import path so the figure modules can load [supp_style.py](supp_style.py), which sets the shared fonts and colours. Its `canvas(height_mm)` creates a 175 mm wide figure, and `panel(ax, letter, title, ...)` adds a panel letter and title. S10, S12 and S13 also require Biopython. Their shared renderer takes an existing tree with branch lengths and node support, the target tip, annotations for every tip, the order of the two support metrics and the page title. S6 uses ggplot2 and patchwork. S11 uses ggplot2 and gggenes.

This directory contains drawing code. Study records, fitted estimates, reference trees and finished panels are held separately. The caller supplies study-specific memberships, counts, result text and artwork identifiers where the drawing functions expect them. Sequence inspection, alignment, tree inference and model fitting are outside the scope of these modules.
