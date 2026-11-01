
# 05: leave-one-out DESeq2 (each sample removed in turn, model refitted)
source("R/config.R")
need("DESeq2")
inp <- readRDS("data/processed/input.rds")
counts <- inp$counts; coldata <- inp$coldata; annot <- inp$annot
res <- readRDS("results/rds/deseq_results.rds")$res
full_deg <- res$ensembl[res$sig != "NS"]

genes <- c(CXCL8 = "CXCL8", NFKBIA = "NFKBIA", IL1B = "IL1B", CLDN2 = "CLDN2")
gid <- annot$ensembl[match(genes, annot$symbol)]; stopifnot(!anyNA(gid))
names(gid) <- names(genes)

loo <- do.call(rbind, lapply(colnames(counts), function(s) {
  m  <- counts[, colnames(counts) != s]; m <- m[rowSums(m) >= MIN_TOTAL_COUNT, ]
  dd <- DESeq2::DESeq(DESeq2::DESeqDataSetFromMatrix(m, coldata[colnames(m), ], design = ~ group), quiet = TRUE)
  r  <- as.data.frame(DESeq2::results(dd, contrast = c("group", "adult", "infant")))
  isdeg <- function(id) !is.na(r[id, "padj"]) && r[id, "padj"] < PADJ_CUT && abs(r[id, "log2FoldChange"]) > LFC_CUT
  ids <- rownames(r)[which(r$padj < PADJ_CUT & abs(r$log2FoldChange) > LFC_CUT)]
  out <- data.frame(removed = s, n_DEG = length(ids), kept_of_1202 = sum(full_deg %in% ids))
  for (g in names(gid)) out[[g]] <- ifelse(isdeg(gid[[g]]), "DEG", "ns")
  out
}))
print(loo)
utils::write.csv(loo, "results/tables/Table4_leave_one_out.csv", row.names = FALSE)
expected <- c(J1005 = 695, J1006 = 707, J1009 = 2028, J11 = 1209, J2 = 911, J3 = 811)
for (s in loo$removed) record(paste0("LOO_nDEG_without_", s), loo$n_DEG[loo$removed == s], expected[[s]])
record("CLDN2_DEG_in_n_of_6", sum(loo$CLDN2 == "DEG"), 6)
record("CXCL8_DEG_in_n_of_6", sum(loo$CXCL8 == "DEG"), 5)
record("NFKBIA_DEG_in_n_of_6", sum(loo$NFKBIA == "DEG"), 3)
record("IL1B_DEG_in_n_of_6", sum(loo$IL1B == "DEG"), 1)

