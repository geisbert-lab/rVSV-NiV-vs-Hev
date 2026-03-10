# Transcriptomic analyses for _A single-cycle, recombinant VSV platform Nipah vaccine cross-protects against Hendra virus in nonhuman primates_

> Paper citation TBD; manuscript in preparation

## Methods

The expression of ~779 host mRNAs was quantified via the Nanostring NHP Immunology v2 panel according to the manufacturer’s instructions. Raw RCC files were loaded into [nSolver v4.0](https://nanostring.com/products/ncounter-analysis-system/ncounter-analysis-solutions/), and background thresholding was performed using the default parameters. Thresholded count matrices were exported from nSolver and analyzed with `[limma v3.65.3](https://doi.org/10.1093/nar/gkv007)` (PMID [25605792](https://pubmed.ncbi.nlm.nih.gov/25605792/) via `[edgeR v4.7.3](https://doi.org/10.1093/bioinformatics/btp616)` (PMID [19910308](https://pubmed.ncbi.nlm.nih.gov/19910308/)) in R v4.5.0.

Three differential expression (DE) analyses were performed. First, the host response in unvaccinated controls was evaluated by comparing 4 DPI and terminal (6-7 DPI) to baseline (0 DPI), binning by days postinfection. Second, the host response in vaccinated animals was evaluated by comparing 4, 7, and 10 DPI samples to baseline. Finally, vaccinated animals were compared to controls at baseline, 4 DPI, and 7 DPI/terminal to quantify diverging host response patterns. For each comparison, significantly DE genes were determined by an FDR-adjusted p-value < 0.05 and a log2 fold change >1 or <-1. 

Digital cell type deconvolution was performed on all Nanostring-profiled samples via `[CIBERSORTx](https://doi.org/10.1038/s41587-019-0114-2)` (PMID [31061481](https://pubmed.ncbi.nlm.nih.gov/31061481/)) run in absolute mode using the `LM22` signature matrix, B-mode batch correction, 1000 permutations, and the CPM matrix as input. Significance comparisons in control and vaccinated time-series analyses were performed for each cell type via Friedman rank-sum test with post-hoc Dunn's test. FDR-adjusted p-values < 0.05 were considered significant. Significance comparisons between treatment groups (e.g., controls at 0 DPI vs. vaccinated at 0 DPI) were performed via Mann-Whitney U-test for each cell type, with p-values < 0.05 considered significant. Statistical comparisons were performed via `[rstatix v0.7.2](https://rpkgs.datanovia.com/rstatix/)` in R v4.5.0. 

## Data availability

The raw RCC files are available via NCBI GEO (_accession TBD_). The unformatted output figures are in [`analysis`](analysis/). The code is in [`analysis.r`](analysis.r). 

