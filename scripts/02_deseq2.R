
# 02: filtering, DESeq2 (adult vs infant), VST, PCA, selected genes
source("R/config.R")
need(c("DESeq2", "SummarizedExperiment"))

inp <- readRDS("data/processed/input.rds")
counts <- inp$counts; coldata <- inp$coldata; annot <- inp$annot

keep <- rowSums(counts) >= MIN_TOTAL_COUNT
counts_f <- counts[keep, ]
record("genes_after_filter(total>=10)", nrow(counts_f), 20058)
record("genes_if_>=10_in_>=2_samples", sum(rowSums(counts >= 10) >= 2), NA)  # 15224: NOT the filter that reproduces 20,058

dds <- DESeq2::DESeqDataSetFromMatrix(counts_f, coldata, design = ~ group)
dds <- DESeq2::DESeq(dds)
res <- as.data.frame(DESeq2::results(dds, contrast = c("group", "adult", "infant")))  # + log2FC = higher in adult
res$ensembl <- rownames(res)
res$symbol  <- annot$symbol[match(res$ensembl, annot$ensembl)]
res$entrez  <- annot$entrez[match(res$ensembl, annot$ensembl)]
res$sig <- ifelse(!is.na(res$padj) & res$padj < PADJ_CUT & res$log2FoldChange >  LFC_CUT, "up_in_adult",
           ifelse(!is.na(res$padj) & res$padj < PADJ_CUT & res$log2FoldChange < -LFC_CUT, "up_in_infant", "NS"))
res$rank_padj <- match(res$ensembl, res$ensembl[order(res$padj)])

record("n_DEG", sum(res$sig != "NS"), 1202)
record("DEG_up_in_adult", sum(res$sig == "up_in_adult"), 640)
record("DEG_up_in_infant", sum(res$sig == "up_in_infant"), 562)
top10 <- res$symbol[order(res$padj)][1:10]
record("top10_by_padj", paste(top10, collapse = ";"))

# ---- VST and PCA (top 500 most variable genes; identical to DESeq2::plotPCA)
vsd  <- DESeq2::vst(dds, blind = TRUE)
vmat <- SummarizedExperiment::assay(vsd)
sel  <- order(apply(vmat, 1, stats::var), decreasing = TRUE)[1:500]
pc   <- stats::prcomp(t(vmat[sel, ]))
pct  <- 100 * pc$sdev^2 / sum(pc$sdev^2)
pca  <- data.frame(sample = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                   group = coldata[rownames(pc$x), "group"])
cen <- stats::aggregate(cbind(PC1, PC2) ~ group, pca, mean)
pca$dist_to_centroid <- sqrt((pca$PC1 - cen$PC1[match(pca$group, cen$group)])^2 +
                             (pca$PC2 - cen$PC2[match(pca$group, cen$group)])^2)
record("PC1_percent", round(pct[1], 1), 71.6, 0.1)
record("PC2_percent", round(pct[2], 1), 11.3, 0.1)
record("PC1+PC2_percent", round(pct[1] + pct[2], 1), 82.9, 0.2)
record("J1009_distance_to_centroid", round(pca$dist_to_centroid[pca$sample == "J1009"], 1), 32.4, 0.3)

# ---- selected genes (Table 3)
sel_genes <- c("PYGL", "CLDN2", "FABP6", "CXCL8", "NFKBIA", "IL1B", "TNF", "RELB")
t3 <- res[res$symbol %in% sel_genes, c("symbol", "ensembl", "baseMean", "log2FoldChange", "padj", "rank_padj", "sig")]
t3 <- t3[order(match(t3$symbol, sel_genes)), ]
print(t3)
utils::write.csv(t3, "results/tables/Table3_selected_genes.csv", row.names = FALSE)
record("CXCL8_rank_by_padj", res$rank_padj[which(res$symbol == "CXCL8")[1]], 509)
record("NFKBIA_log2FC", round(res$log2FoldChange[which(res$symbol == "NFKBIA")[1]], 2), 1.02, 0.01)

utils::write.csv(res, "results/tables/DESeq2_all_genes.csv", row.names = FALSE)
utils::write.csv(res[res$sig != "NS", ], "results/tables/DEGs_1202.csv", row.names = FALSE)
utils::write.csv(pca, "results/tables/PCA_coordinates.csv", row.names = FALSE)
saveRDS(dds,  "results/rds/dds.rds")
saveRDS(vmat, "results/rds/vst_matrix.rds")
saveRDS(list(res = res, pca = pca, pct = pct), "results/rds/deseq_results.rds")

