
# Central parameters. All scripts start with source("R/config.R").
if (!file.exists("R/config.R")) stop("Open enteroid_reanalysis.Rproj first (working directory must be the project root).")

GEO_ID          <- "GSE227205"
PADJ_CUT        <- 0.05
LFC_CUT         <- 1
MIN_TOTAL_COUNT <- 10       # gene kept if total count over all six samples >= 10
SOFT_POWER      <- 12       # power used in the original analysis
SOFT_POWER_SENS <- 18       # sensitivity analysis (WGCNA FAQ suggests 18 for signed, n < 20)
MIN_MODULE      <- 30
MERGE_CUT       <- 0.25
DEEP_SPLIT      <- 2
MAX_BLOCK       <- 5000     # 8 GB RAM -> 5000. With >= 32 GB RAM use 20058 (single block).
WGCNA_THREADS   <- 1        # keep 1 on small machines
SEED            <- 1

for (d in c("data/raw", "data/processed", "results/tables", "results/figures", "results/rds"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

source("R/helpers.R")

