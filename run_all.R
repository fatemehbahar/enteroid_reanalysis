
# Runs the whole pipeline. Open enteroid_reanalysis.Rproj first.
stopifnot(file.exists("R/config.R"))
unlink("results/key_numbers.csv")
steps <- c("scripts/01_data.R", "scripts/02_deseq2.R", "scripts/03_enrichment.R", "scripts/04_wgcna.R",
           "scripts/05_robustness.R", "scripts/06_figures.R", "scripts/07_session.R")
for (s in steps) {
  message("\n=== ", s, " ===")
  sys.source(s, envir = new.env(parent = globalenv()))
}
kn <- read.csv("results/key_numbers.csv")
print(kn, row.names = FALSE)
message("\nItems with status DIFFERS: ", sum(kn$status == "DIFFERS"), " (see results/key_numbers.csv)")

