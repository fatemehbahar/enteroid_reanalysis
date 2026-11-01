
need <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "),
                         ". Run scripts/00_setup.R first.", call. = FALSE)
}

# Print and append one number to results/key_numbers.csv.
# 'expected' = value reported in the earlier verified analysis; status OK/DIFFERS within tol.
record <- function(item, value, expected = NA, tol = 0) {
  value_chr <- as.character(value)
  status <- ""
  if (!is.na(expected)) {
    v <- suppressWarnings(as.numeric(value))
    status <- if (is.na(v)) "DIFFERS" else if (abs(v - expected) <= tol) "OK" else "DIFFERS"
  }
  row <- data.frame(item = item, value = value_chr, earlier_value = expected, status = status)
  f <- "results/key_numbers.csv"
  utils::write.table(row, f, sep = ",", row.names = FALSE,
                     col.names = !file.exists(f), append = file.exists(f))
  message(sprintf("%-46s %-14s %-14s %s", item, value_chr,
                  ifelse(is.na(expected), "", paste0("earlier:", expected)), status))
  invisible(row)
}

save_gg <- function(p, name, w = 7, h = 5) {
  ggplot2::ggsave(file.path("results/figures", paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path("results/figures", paste0(name, ".pdf")), p, width = w, height = h)
}

save_base <- function(name, fun, w = 7, h = 5) {
  grDevices::png(file.path("results/figures", paste0(name, ".png")), width = w, height = h, units = "in", res = 300)
  fun(); grDevices::dev.off()
  grDevices::pdf(file.path("results/figures", paste0(name, ".pdf")), width = w, height = h)
  fun(); grDevices::dev.off()
}

get_hallmark <- function() {
  h <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
                error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
  as.data.frame(h)[, c("gs_name", "gene_symbol")]
}

safe_df <- function(x) if (is.null(x)) data.frame() else as.data.frame(x)

keep_sig <- function(x, cut = 0.05) {
  if (is.null(x)) return(NULL)
  x@result <- x@result[x@result$p.adjust < cut, , drop = FALSE]
  if (nrow(x@result) == 0) NULL else x
}

