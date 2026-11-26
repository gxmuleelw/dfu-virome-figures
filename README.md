[English](README.md) | [简体中文](README.zh-CN.md)

# Gut virome and diabetic foot ulcer severity

Code for the manuscript “Lower relative abundance of lytic gut phages in patients with Wagner grade 4 diabetic foot ulcers”.

We compared the gut virome of patients with Wagner grade 4 diabetic foot ulcers with that of patients with grades 1–3. For each figure, the analysis link opens the code that calculates the quantities shown; the drawing link opens the code that displays them. The repository contains code and documentation. Study data, fitted results and original artwork are held separately.

| Figure and research question | Analysis or descriptive summary | Drawing |
| --- | --- | --- |
| **Figure 1b–f.** What does the viral catalogue contain, and which patients contributed to each analysis? | [Count catalogue categories, detections and eligible patients](main/figure1_catalogue_analysis.py) | [Catalogue and cohort panels](main/figure1_reflow.py) |
| **Figure 2.** Do predicted phage lifestyle proportions differ by Wagner grade? | [Calculate within-patient percentages, fit adjusted group contrasts and repeat the sensitivity comparisons](main/figure2_lifestyle_analysis.R) | [Distributions and intervals](main/figure2_render.py) |
| **Figure 3.** How does viral abundance covary with the predicted host genus compared with other bacterial genera? | [Form genus and member pairs, permute patient labels, compute genus-level bootstrap intervals and compare matched with nonmatched pairs](main/figure3_paired_analysis.py) | [Correlations and counts](main/figure3_render.py) |
| **Figure 4.** How do detection and abundance of the *Escherichia* and *Collinsella* viral collections relate to severity and bacterial abundance? | [Rank paired abundances, test detection, construct the bacterial score and model abundance among positive samples](main/figure4_host_analysis.R) | [Detection and abundance panels](main/figure4_render.py) |
| **Figure 5.** What evidence do existing annotations and trees provide for the focal viral representative? | [Summarize marker coordinates, compare reviewed splits and select neighbouring tips on a finished tree](main/figure5_analysis.py) | [Layout from existing vector artwork](main/figure5_reflow.py) |

Figure 1a is a study-design diagram drawn in PowerPoint.

## Directory guides

- [Main figures](main/README.md): panel-level methods, inputs, function names and drawing commands.
- [Supplementary figures](supp/README.md): plotting code for Figures S1–S14 and shared drawing functions.
- [Tables](table_analysis/README.md): sources of the values in Tables 1 and 2, and the calculations behind table metrics.

## Data and software

Raw metagenomic sequencing data are available in the NCBI Sequence Read Archive under accession PRJNA1471479. Statistical analyses used R 4.5.0 and Python 3.12.10. Package versions are listed in the Key resources table (Supplementary Methods, section 11).
