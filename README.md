
# Reanalysis of GSE227205: infant vs adult human jejunal enteroids

Exploratory reanalysis of public bulk RNA-seq data (3 infant-derived and 3 adult-derived jejunal enteroid lines;
Adeniyi-Ipadeola et al., mBio 2024) with DESeq2, enrichment analysis (ORA and GSEA) and signed WGCNA.
The contrast is **developmental stage (infant vs adult)**, not chronological aging. With n = 3 per group all
network results are exploratory.

## Run
1. Open `enteroid_reanalysis.Rproj` in RStudio (developed with R 4.6.1 / Bioconductor 3.23).
2. `source("scripts/00_setup.R")` (once; installs packages).
3. `source("run_all.R")`, or run `scripts/01` to `07` one by one with `source("scripts/0X_name.R")`.

Internet is needed for the GEO download and for KEGG. Runtime is about 30-60 min, mostly WGCNA.
RAM: `MAX_BLOCK <- 5000` in `R/config.R` suits 8 GB; set it to 20058 with >= 32 GB for a single-block network.

## Outputs
- `results/key_numbers.csv`: every number used in the manuscript (`earlier_value` = value from the first analysis, `status` = OK/DIFFERS).
- `results/tables/`: Tables 1-4, DEG list, enrichment results, module and hub tables.
- `results/figures/`: Figures 1-9 (PNG 300 dpi and PDF).
- `results/sessionInfo.txt`: exact software versions.

## Scripts
| script | content |
|---|---|
| 01_data.R | GEO download, sample groups from GEO metadata, gene annotation |
| 02_deseq2.R | filter (total count >= 10), DESeq2 adult vs infant, VST, PCA |
| 03_enrichment.R | ORA (adult-up / infant-up) and GSEA on the Wald statistic |
| 04_wgcna.R | signed WGCNA (power 12), module-trait tests (Student, exact permutation, BH), hubs, sensitivity analyses |
| 05_robustness.R | leave-one-out DESeq2 |
| 06_figures.R | figures |
| 07_session.R | software versions |

## Notes
- The GEO `GeneSymbol` column is corrupted by spreadsheet date conversion (e.g. "1-Dec"); symbols are re-derived from Ensembl IDs.
- The filter that reproduces 20,058 genes is total count >= 10 over the six samples.
- With 3 vs 3 samples only 20 balanced label permutations exist, so the smallest exact two-sided permutation p is 0.1.
- No WGCNA module is significant after Benjamini-Hochberg correction; the two main modules largely restate the DEG lists.
- The random seed is fixed (`SEED` in `R/config.R`).

## Data and license
Data: GEO accession GSE227205 (downloaded by the scripts, not stored here). Code: MIT (see LICENSE).


## Status
Manuscript in preparation. All numbers in `results/key_numbers.csv` are reproduced by `run_all.R`.
