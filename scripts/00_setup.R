
# Run once. Installs all packages. Takes 10-30 min the first time.
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran <- c("ggplot2", "ggrepel", "pheatmap", "cowplot", "msigdbr", "dplyr")
bioc <- c("GEOquery", "Biobase", "DESeq2", "SummarizedExperiment", "AnnotationDbi",
          "org.Hs.eg.db", "clusterProfiler", "enrichplot", "WGCNA")

install.packages(setdiff(cran, rownames(installed.packages())))
BiocManager::install(setdiff(bioc, rownames(installed.packages())), update = FALSE, ask = FALSE)

# msigdbr >= 10 may ask for the msigdbdf data package. If scripts/03 complains, run:
# install.packages("msigdbdf", repos = c("https://igordot.r-universe.dev", "https://cloud.r-project.org"))

print(sapply(c(cran, bioc), requireNamespace, quietly = TRUE))   # all must be TRUE

