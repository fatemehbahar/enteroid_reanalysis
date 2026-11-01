
# 01: download GSE227205 metadata + gene counts; build sample table and gene annotation
source("R/config.R")
need(c("GEOquery", "Biobase", "AnnotationDbi", "org.Hs.eg.db"))
options(timeout = 600)

# ---- sample metadata (group is read from GEO "developmental stage", not typed by hand)
gse  <- GEOquery::getGEO(GEO_ID, GSEMatrix = TRUE, destdir = "data/raw")[[1]]
meta <- Biobase::pData(gse)
utils::write.csv(meta, "data/processed/geo_metadata.csv")

stage_col <- names(meta)[vapply(meta, function(x) any(grepl("developmental stage", x, ignore.case = TRUE)), logical(1))][1]
stopifnot(!is.na(stage_col))
stage <- tolower(trimws(sub(".*developmental stage:\\s*", "", as.character(meta[[stage_col]]), ignore.case = TRUE)))
coldata <- data.frame(sample = meta$title, gsm = rownames(meta),
                      group = factor(stage, levels = c("infant", "adult")),
                      row.names = meta$title, stringsAsFactors = FALSE)
stopifnot(!anyNA(coldata$group), all(table(coldata$group) == 3))
print(coldata)

# ---- counts (GEO supplementary file)
f <- list.files("data/raw", pattern = "featureCounts", recursive = TRUE, full.names = TRUE)
if (length(f) == 0) {
  GEOquery::getGEOSuppFiles(GEO_ID, baseDir = "data/raw", makeDirectory = TRUE)
  f <- list.files("data/raw", pattern = "featureCounts", recursive = TRUE, full.names = TRUE)
}
stopifnot(length(f) >= 1)
raw <- utils::read.csv(f[1], check.names = FALSE)
stopifnot(all(c("GeneID", coldata$sample) %in% colnames(raw)), !anyDuplicated(raw$GeneID))

counts <- as.matrix(raw[, coldata$sample]); rownames(counts) <- raw$GeneID
storage.mode(counts) <- "integer"

# The GeneSymbol column of the GEO file is corrupted by spreadsheet date conversion (1-Dec, 1-Mar ...),
# so symbols are re-derived from Ensembl IDs.
db <- org.Hs.eg.db::org.Hs.eg.db
annot <- data.frame(
  ensembl = rownames(counts),
  symbol  = unname(AnnotationDbi::mapIds(db, rownames(counts), "SYMBOL",   "ENSEMBL", multiVals = "first")),
  entrez  = unname(AnnotationDbi::mapIds(db, rownames(counts), "ENTREZID", "ENSEMBL", multiVals = "first")),
  biotype = raw$GeneBiotype, stringsAsFactors = FALSE)

record("genes_in_GEO_file", nrow(counts), 22445)
lib <- round(colSums(counts) / 1e6, 1)
for (s in names(lib)) record(paste0("library_size_million_", s), lib[[s]])

saveRDS(list(counts = counts, coldata = coldata, annot = annot), "data/processed/input.rds")

