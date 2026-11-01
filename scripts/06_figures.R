
# 06: figures 1-9 (PNG 300 dpi + PDF) in results/figures
source("R/config.R")
need(c("ggplot2", "ggrepel", "pheatmap", "cowplot", "enrichplot", "WGCNA", "SummarizedExperiment"))
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel); library(WGCNA) })
theme_set(theme_minimal(base_size = 12))
cols_grp <- c(infant = "#E31A1C", adult = "#1F78B4")

inp <- readRDS("data/processed/input.rds"); coldata <- inp$coldata
dr  <- readRDS("results/rds/deseq_results.rds"); res <- dr$res; pca <- dr$pca; pct <- dr$pct
dds <- readRDS("results/rds/dds.rds"); vmat <- readRDS("results/rds/vst_matrix.rds")
ora <- readRDS("results/rds/ora.rds")
w <- readRDS("results/rds/wgcna_power12.rds")
net <- w$net; mt <- w$mt; hub <- w$hub; fi <- w$fi

# Fig 1: distributions before / after VST (all genes)
raw_counts <- SummarizedExperiment::assay(dds, "counts")
save_base("Fig1_QC_boxplots", function() {
  op <- graphics::par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))
  graphics::boxplot(log2(raw_counts + 1), las = 2, col = "grey80", main = "Raw counts", ylab = "log2(count + 1)")
  graphics::boxplot(vmat, las = 2, col = "grey80", main = "After VST", ylab = "VST value")
  graphics::par(op)
}, w = 9, h = 4.5)

# Fig 2: PCA
p2 <- ggplot(pca, aes(PC1, PC2, colour = group, label = sample)) + geom_point(size = 4) +
  geom_text_repel(show.legend = FALSE) + scale_colour_manual(values = cols_grp) +
  labs(x = sprintf("PC1: %.1f%% variance", pct[1]), y = sprintf("PC2: %.1f%% variance", pct[2]), colour = "Group",
       title = "PCA, 500 most variable genes")
save_gg(p2, "Fig2_PCA", 6.5, 5)

# Fig 3: volcano with gene symbols
rv <- res[!is.na(res$padj), ]
rv$label <- ifelse(is.na(rv$symbol), rv$ensembl, rv$symbol)
deg <- rv[rv$sig != "NS", ]; deg <- deg[order(deg$padj), ]
lab <- unique(rbind(utils::head(deg, 10), rv[rv$symbol %in% c("CLDN2", "CXCL8", "NFKBIA"), ]))
p3 <- ggplot(rv, aes(log2FoldChange, -log10(padj), colour = sig)) + geom_point(alpha = .5, size = 1) +
  geom_text_repel(data = lab, aes(label = label), colour = "black", size = 3, max.overlaps = Inf, fontface = "italic") +
  scale_colour_manual(values = c(up_in_adult = "#E31A1C", up_in_infant = "#1F78B4", NS = "grey75"),
                      labels = c(up_in_adult = "Higher in adult", up_in_infant = "Higher in infant", NS = "NS")) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = 2) + geom_hline(yintercept = -log10(PADJ_CUT), linetype = 2) +
  labs(x = "log2 fold change (adult / infant)", y = "-log10 adjusted P", colour = NULL)
save_gg(p3, "Fig3_volcano", 8, 6)

# Fig 4: GO-BP, both groups (significant terms only)
dp <- function(x, title) { x <- keep_sig(x); if (is.null(x)) NULL else enrichplot::dotplot(x, showCategory = 12, title = title) }
pl4 <- Filter(Negate(is.null), list(dp(ora$adult$go, "GO-BP: higher in adult"), dp(ora$infant$go, "GO-BP: higher in infant")))
save_gg(cowplot::plot_grid(plotlist = pl4, ncol = 1, labels = "AUTO"), "Fig4_GO_BP", 9, 11)

# Fig 5: KEGG and Hallmark
pl5 <- Filter(Negate(is.null), list(dp(ora$adult$kegg, "KEGG: higher in adult"), dp(ora$adult$hall, "Hallmark: higher in adult"),
                                    dp(ora$infant$kegg, "KEGG: higher in infant")))
save_gg(cowplot::plot_grid(plotlist = pl5, ncol = 1, labels = "AUTO"), "Fig5_KEGG_Hallmark", 9, 4 * length(pl5))

# Fig 6: soft-threshold scan (line at R2 = 0.8)
pa <- ggplot(fi, aes(Power, signed_R2, label = Power)) + geom_text(colour = "red") + geom_hline(yintercept = 0.8, colour = "blue") +
  labs(y = "Scale-free topology fit (signed R2)", title = "Scale independence")
pb <- ggplot(fi, aes(Power, mean.k., label = Power)) + geom_text(colour = "red") + labs(y = "Mean connectivity", title = "Mean connectivity")
save_gg(cowplot::plot_grid(pa, pb, nrow = 1), "Fig6_soft_threshold", 10, 4.5)

# Fig 7: dendrograms of the WGCNA blocks (one page / file per block)
cols <- WGCNA::labels2colors(net$colors)
draw_block <- function(b) WGCNA::plotDendroAndColors(net$dendrograms[[b]], cols[net$blockGenes[[b]]], "Module",
  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = paste("Gene dendrogram, block", b))
grDevices::pdf("results/figures/Fig7_dendrograms.pdf", width = 12, height = 5)
for (b in seq_along(net$dendrograms)) draw_block(b)
grDevices::dev.off()
for (b in seq_along(net$dendrograms)) {
  grDevices::png(sprintf("results/figures/Fig7_dendrogram_block%d.png", b), width = 12, height = 5, units = "in", res = 200)
  draw_block(b); grDevices::dev.off()
}

# Fig 8: 25 highest-kME genes of each main module (illustration only, not regulatory hubs)
top_of <- function(m, k = 25) { h <- hub[hub$module == m, ]; utils::head(h[order(-h$kME_own), ], k) }
hh <- rbind(top_of(w$adult_mod), top_of(w$infant_mod))
hh$label <- make.unique(ifelse(is.na(hh$symbol), hh$ensembl, hh$symbol))
mm <- vmat[hh$ensembl, ]; rownames(mm) <- hh$label
ann_col <- data.frame(Group = coldata[colnames(mm), "group"], row.names = colnames(mm))
ann_row <- data.frame(Module = rep(c("adult-associated", "infant-associated"), each = 25), row.names = hh$label)
for (ext in c("png", "pdf"))
  pheatmap::pheatmap(mm, scale = "row", annotation_col = ann_col, annotation_row = ann_row, fontsize_row = 7,
                     filename = paste0("results/figures/Fig8_heatmap_kME.", ext), width = 6, height = 9)

# Fig 9: module-trait relationships
m9 <- mt; m9$module <- factor(m9$module, levels = m9$module[order(m9$r)])
p9 <- ggplot(m9, aes("Adult vs infant", module, fill = r)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f (p=%.3f)", r, p_student)), size = 2.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1)) +
  labs(x = NULL, y = NULL, caption = "Student p, uncorrected. No module is significant after Benjamini-Hochberg correction.")
save_gg(p9, "Fig9_module_trait", 5, 11)
message("Figures written to results/figures")

