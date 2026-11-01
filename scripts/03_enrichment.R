
# 03: ORA (adult-up, infant-up) and GSEA. Universe = all tested genes.
source("R/config.R")
need(c("clusterProfiler", "org.Hs.eg.db", "msigdbr"))

res   <- readRDS("results/rds/deseq_results.rds")$res
db    <- org.Hs.eg.db::org.Hs.eg.db
t2g   <- get_hallmark()
tested <- !is.na(res$padj)
univ_ent <- unique(stats::na.omit(res$entrez[tested]))
univ_sym <- unique(stats::na.omit(res$symbol[tested]))

run_ora <- function(label) {
  ent <- unique(stats::na.omit(res$entrez[res$sig == label]))
  sym <- unique(stats::na.omit(res$symbol[res$sig == label]))
  list(
    go   = clusterProfiler::enrichGO(ent, OrgDb = db, keyType = "ENTREZID", ont = "BP", universe = univ_ent,
                                     pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                                     minGSSize = 10, maxGSSize = 5000, readable = TRUE),
    kegg = tryCatch(clusterProfiler::enrichKEGG(ent, organism = "hsa", universe = univ_ent,
                                                pvalueCutoff = 1, qvalueCutoff = 1),
                    error = function(e) { message("KEGG needs internet: ", e$message); NULL }),
    hall = clusterProfiler::enricher(sym, TERM2GENE = t2g, universe = univ_sym,
                                     pvalueCutoff = 1, qvalueCutoff = 1))
}
ora <- list(adult = run_ora("up_in_adult"), infant = run_ora("up_in_infant"))

for (g in names(ora)) for (d in names(ora[[g]])) {
  df <- safe_df(ora[[g]][[d]])
  if (nrow(df)) utils::write.csv(df[df$p.adjust < 0.05, ], sprintf("results/tables/ORA_%s_up_%s.csv", g, d), row.names = FALSE)
}

# ---- Table 1 (adult-up ORA)
pick <- function(df, id) {
  if (!nrow(df) || !(id %in% df$ID)) return(data.frame(ID = id, Count = NA, qvalue = NA))
  df[df$ID == id, c("ID", "Count", "qvalue")][1, ]
}
ids <- list(hall = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE"),
            kegg = c("hsa04668", "hsa04064"), go = c("GO:0001817", "GO:0002237"))
t1 <- do.call(rbind, c(lapply(ids$hall, function(i) pick(safe_df(ora$adult$hall), i)),
                       lapply(ids$kegg, function(i) pick(safe_df(ora$adult$kegg), i)),
                       lapply(ids$go,   function(i) pick(safe_df(ora$adult$go),   i))))
print(t1)
utils::write.csv(t1, "results/tables/Table1_ORA_adult_up.csv", row.names = FALSE)
cnt <- function(x, i) { v <- x$Count[x$ID == i]; if (length(v)) v else NA }
record("ORA_TNFA_via_NFKB_count", cnt(t1, "HALLMARK_TNFA_SIGNALING_VIA_NFKB"), 30)
record("ORA_inflammatory_response_count", cnt(t1, "HALLMARK_INFLAMMATORY_RESPONSE"), 28)
record("ORA_KEGG_TNF_count", cnt(t1, "hsa04668"), 13)
record("ORA_KEGG_NFkB_count", cnt(t1, "hsa04064"), 12)
record("ORA_GO0001817_count", cnt(t1, "GO:0001817"), 42)
record("ORA_GO0002237_count", cnt(t1, "GO:0002237"), 20)

# ---- GSEA on the Wald statistic (positive = adult)
set.seed(SEED)
r1 <- res[!is.na(res$stat) & !is.na(res$symbol) & !duplicated(res$symbol), ]
rk_sym <- sort(stats::setNames(r1$stat, r1$symbol), decreasing = TRUE)
r2 <- res[!is.na(res$stat) & !is.na(res$entrez) & !duplicated(res$entrez), ]
rk_ent <- sort(stats::setNames(r2$stat, r2$entrez), decreasing = TRUE)
gs_h <- clusterProfiler::GSEA(rk_sym, TERM2GENE = t2g, pvalueCutoff = 1, eps = 0, seed = TRUE, verbose = FALSE)
gs_k <- tryCatch(clusterProfiler::gseKEGG(rk_ent, organism = "hsa", pvalueCutoff = 1, eps = 0, seed = TRUE, verbose = FALSE),
                 error = function(e) { message("KEGG needs internet: ", e$message); NULL })
gh <- safe_df(gs_h); gk <- safe_df(gs_k)
gsea_tab <- rbind(gh[gh$ID %in% ids$hall, c("ID", "setSize", "NES", "qvalue")],
                  if (nrow(gk)) gk[gk$ID %in% ids$kegg, c("ID", "setSize", "NES", "qvalue")])
print(gsea_tab)
utils::write.csv(gsea_tab, "results/tables/Table1_GSEA.csv", row.names = FALSE)
nes <- function(i) if (i %in% gsea_tab$ID) round(gsea_tab$NES[gsea_tab$ID == i], 2) else NA
record("GSEA_NES_TNFA_via_NFKB", nes("HALLMARK_TNFA_SIGNALING_VIA_NFKB"), 2.46, 0.1)
record("GSEA_NES_inflammatory_response", nes("HALLMARK_INFLAMMATORY_RESPONSE"), 2.31, 0.1)
record("GSEA_NES_KEGG_TNF", nes("hsa04668"), 2.33, 0.1)
record("GSEA_NES_KEGG_NFkB", nes("hsa04064"), 2.05, 0.1)

saveRDS(ora, "results/rds/ora.rds")
utils::write.csv(gh, "results/tables/GSEA_hallmark_all.csv", row.names = FALSE)

