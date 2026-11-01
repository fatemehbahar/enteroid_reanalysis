
# 04: signed WGCNA on all genes (power 12), module-trait tests, hubs, module enrichment,
#     sensitivity (power 18; 5,000 most variable genes)
source("R/config.R")
need(c("WGCNA", "clusterProfiler", "org.Hs.eg.db", "msigdbr"))
suppressPackageStartupMessages(library(WGCNA))   # WGCNA masks cor(); we always call stats::cor explicitly
options(stringsAsFactors = FALSE)
if (WGCNA_THREADS > 1) WGCNA::enableWGCNAThreads(WGCNA_THREADS) else WGCNA::disableWGCNAThreads()

inp <- readRDS("data/processed/input.rds"); coldata <- inp$coldata; annot <- inp$annot
res <- readRDS("results/rds/deseq_results.rds")$res
datExpr <- t(readRDS("results/rds/vst_matrix.rds"))
stopifnot(WGCNA::goodSamplesGenes(datExpr, verbose = 0)$allOK)
trait <- as.numeric(coldata[rownames(datExpr), "group"] == "adult")
n <- nrow(datExpr)
record("WGCNA_genes", ncol(datExpr), 20058)

# ---- soft threshold
sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = 1:20, networkType = "signed", verbose = 0)
fi  <- sft$fitIndices
fi$signed_R2 <- -sign(fi$slope) * fi$SFT.R.sq
utils::write.csv(fi, "results/tables/soft_threshold_primary.csv", row.names = FALSE)
record("scale_free_R2_power12", round(fi$signed_R2[fi$Power == 12], 3), 0.505, 0.003)
record("scale_free_R2_power20", round(fi$signed_R2[fi$Power == 20], 2), 0.70, 0.01)
record("max_R2_up_to_power20", round(max(fi$signed_R2), 2), NA)
record("powerEstimate_R2>=0.85", ifelse(is.na(sft$powerEstimate), "none", sft$powerEstimate))

# ---- helpers
build_net <- function(power) {
  set.seed(SEED)
  WGCNA::blockwiseModules(datExpr, power = power, networkType = "signed", TOMType = "signed",
                          minModuleSize = MIN_MODULE, mergeCutHeight = MERGE_CUT, deepSplit = DEEP_SPLIT,
                          maxBlockSize = MAX_BLOCK, numericLabels = TRUE, saveTOMs = FALSE, verbose = 0)
}
idx_perm <- utils::combn(n, n / 2)          # 20 balanced label assignments for 3 vs 3
module_trait <- function(net) {
  MEs <- WGCNA::orderMEs(net$MEs)
  r   <- as.numeric(stats::cor(MEs, trait))
  p   <- as.numeric(WGCNA::corPvalueStudent(r, n))
  pp  <- vapply(seq_len(ncol(MEs)), function(j) {
    me <- MEs[[j]]; obs <- stats::cor(me, trait)
    rr <- apply(idx_perm, 2, function(i) { t <- rep(0, n); t[i] <- 1; stats::cor(me, t) })
    mean(abs(rr) >= abs(obs) - 1e-12) }, numeric(1))
  lab <- sub("ME", "", colnames(MEs))
  mt <- data.frame(module = colnames(MEs), label = lab, r = r, p_student = p, p_perm = pp,
                   n_genes = as.numeric(table(net$colors)[lab]))
  mt$p_bh <- stats::p.adjust(mt$p_student, "BH")
  mt[order(-abs(mt$r)), ]
}

# ---- primary network (power 12)
net <- build_net(SOFT_POWER)
mt  <- module_trait(net)
utils::write.csv(mt, "results/tables/Table2_module_trait_power12.csv", row.names = FALSE)
print(utils::head(mt, 6), digits = 3)
adult_mod  <- mt$label[which.max(mt$r)]     # most adult-associated module (positive r)
infant_mod <- mt$label[which.min(mt$r)]     # most infant-associated module (negative r)
mt_get <- function(l, col) mt[[col]][mt$label == l]
record("n_modules", length(unique(net$colors)), 36, 3)
record("adult_module_r", round(mt_get(adult_mod, "r"), 3), 0.869, 0.05)
record("adult_module_p_student", round(mt_get(adult_mod, "p_student"), 4), 0.0247, 0.02)
record("adult_module_genes", mt_get(adult_mod, "n_genes"), 3700, 400)
record("infant_module_r", round(mt_get(infant_mod, "r"), 3), -0.908, 0.05)
record("infant_module_p_student", round(mt_get(infant_mod, "p_student"), 4), 0.0124, 0.02)
record("infant_module_genes", mt_get(infant_mod, "n_genes"), 2617, 400)
record("n_modules_p_student<0.05", sum(mt$p_student < 0.05), 2)
record("min_BH_over_modules", round(min(mt$p_bh), 3), 0.445, 0.1)
record("min_permutation_p", min(mt$p_perm), 0.1)
rest <- mt[!mt$label %in% c(adult_mod, infant_mod), ]
record("next_best_module_r", round(rest$r[which.max(abs(rest$r))], 2), 0.74, 0.1)

# ---- overlap with DEGs
g_adult  <- colnames(datExpr)[net$colors == adult_mod]
g_infant <- colnames(datExpr)[net$colors == infant_mod]
deg_a <- res$ensembl[res$sig == "up_in_adult"]; deg_i <- res$ensembl[res$sig == "up_in_infant"]
record("adult_module_has_adult_up_DEGs", sum(g_adult %in% deg_a), 543, 40)
record("adult_module_has_infant_up_DEGs", sum(g_adult %in% deg_i), 0, 5)
record("infant_module_has_infant_up_DEGs", sum(g_infant %in% deg_i), 420, 40)
record("infant_module_has_adult_up_DEGs", sum(g_infant %in% deg_a), 0, 5)

# ---- hubs (kME in own module)
MEs <- WGCNA::orderMEs(net$MEs)
kME <- as.matrix(WGCNA::signedKME(datExpr, MEs, outputColumnName = "kME"))
own <- kME[cbind(seq_len(nrow(kME)), match(paste0("kME", net$colors), colnames(kME)))]
hub <- data.frame(ensembl = colnames(datExpr), symbol = annot$symbol[match(colnames(datExpr), annot$ensembl)],
                  module = as.character(net$colors), kME_own = own, stringsAsFactors = FALSE)
ha <- hub[hub$module == adult_mod, ]; ha$rank <- rank(-ha$kME_own)
ha <- ha[order(ha$rank), ]
utils::write.csv(ha, "results/tables/hub_ranks_adult_module.csv", row.names = FALSE)
utils::write.csv(utils::head(ha, 30), "results/tables/top30_hubs_adult_module.csv", row.names = FALSE)
for (g in c("NFKBIA", "CXCL1", "CXCL8", "IL1B")) {
  k <- which(ha$symbol == g)[1]
  record(paste0(g, "_rank_in_adult_module"), ha$rank[k], c(NFKBIA = 253, CXCL1 = 560, CXCL8 = 695, IL1B = 813)[[g]], 120)
}
record("adult_module_size", nrow(ha), 3700, 400)
record("fraction_adult_module_kME>0.85", round(mean(ha$kME_own > 0.85), 3), 0.435, 0.05)
record("fraction_all_genes_kME>0.85", round(mean(hub$kME_own > 0.85, na.rm = TRUE), 3), NA)
record("top15_hubs_adult_module", paste(utils::head(ha$symbol, 15), collapse = ";"))
relb <- hub$module[which(hub$symbol == "RELB")[1]]
record("RELB_module", relb); record("RELB_module_r", round(mt_get(relb, "r"), 2), 0.05, 0.15)

# ---- module enrichment (universe = all network genes; default GO size limits)
db <- org.Hs.eg.db::org.Hs.eg.db; t2g <- get_hallmark()
univ_ent <- unique(stats::na.omit(annot$entrez[match(colnames(datExpr), annot$ensembl)]))
univ_sym <- unique(stats::na.omit(annot$symbol[match(colnames(datExpr), annot$ensembl)]))
enrich_module <- function(genes, tag) {
  ent <- unique(stats::na.omit(annot$entrez[match(genes, annot$ensembl)]))
  sym <- unique(stats::na.omit(annot$symbol[match(genes, annot$ensembl)]))
  hm <- safe_df(clusterProfiler::enricher(sym, TERM2GENE = t2g, universe = univ_sym, pvalueCutoff = 1, qvalueCutoff = 1))
  go <- safe_df(clusterProfiler::enrichGO(ent, OrgDb = db, keyType = "ENTREZID", ont = "BP", universe = univ_ent,
                                          pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE))
  utils::write.csv(utils::head(hm[order(hm$p.adjust), ], 20), sprintf("results/tables/module_%s_hallmark.csv", tag), row.names = FALSE)
  utils::write.csv(utils::head(go[order(go$p.adjust), ], 30), sprintf("results/tables/module_%s_GO.csv", tag), row.names = FALSE)
  hm
}
hm_a <- enrich_module(g_adult, "adult"); hm_i <- enrich_module(g_infant, "infant")
record("adult_module_TNFA_count", hm_a$Count[hm_a$ID == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"], 80, 15)
record("infant_module_top_hallmark", hm_i$ID[which.min(hm_i$p.adjust)])

# ---- sensitivity 1: power 18
net18 <- build_net(SOFT_POWER_SENS); mt18 <- module_trait(net18)
utils::write.csv(mt18, "results/tables/module_trait_power18.csv", row.names = FALSE)
record("power18_n_modules", length(unique(net18$colors)), 43, 5)
record("power18_top_positive_r", round(max(mt18$r), 3), 0.878, 0.05)
record("power18_top_negative_r", round(min(mt18$r), 3), -0.904, 0.05)
record("power18_min_BH", round(min(mt18$p_bh), 3), 0.463, 0.1)
tab <- table(net$colors, net18$colors)
record("share_adult_module_kept_in_power18", round(max(tab[adult_mod, ]) / sum(tab[adult_mod, ]), 2), 0.85, 0.1)
record("share_infant_module_kept_in_power18", round(max(tab[infant_mod, ]) / sum(tab[infant_mod, ]), 2), 0.78, 0.1)

# ---- sensitivity 2: 5,000 most variable genes (selection is nearly circular)
top5k <- order(apply(datExpr, 2, stats::var), decreasing = TRUE)[1:5000]
deg_all <- res$ensembl[res$sig != "NS"]
record("DEGs_within_5000_most_variable", sum(deg_all %in% colnames(datExpr)[top5k]), 1199, 10)
sft5 <- WGCNA::pickSoftThreshold(datExpr[, top5k], powerVector = c(1:10, seq(12, 20, 2)), networkType = "signed", verbose = 0)
f5 <- sft5$fitIndices
record("5000genes_R2_power8", round(f5$SFT.R.sq[f5$Power == 8], 2), 0.15, 0.02)
record("5000genes_slope_power8", round(f5$slope[f5$Power == 8], 2), NA)

saveRDS(list(net = net, mt = mt, hub = hub, adult_mod = adult_mod, infant_mod = infant_mod, fi = fi),
        "results/rds/wgcna_power12.rds")
saveRDS(net18, "results/rds/wgcna_power18.rds")

