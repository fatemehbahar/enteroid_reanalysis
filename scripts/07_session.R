
# 07: record software versions (needed for the Methods section)
source("R/config.R")
writeLines(utils::capture.output(utils::sessionInfo()), "results/sessionInfo.txt")
pk <- c("DESeq2", "WGCNA", "clusterProfiler", "org.Hs.eg.db", "msigdbr", "GEOquery", "enrichplot")
for (p in pk) record(paste0("version_", p), as.character(utils::packageVersion(p)))
record("version_R", R.version.string)
record("version_Bioconductor", as.character(BiocManager::version()))

