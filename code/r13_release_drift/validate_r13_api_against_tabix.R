#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(curl)
  library(parallel)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
api_path <- file.path(root, "outputs/R13_API_UNIQUE_VARIANT_EFFECTS.tsv")
query_path <- file.path(root, "outputs/R13_API_QUERY_INDEX.tsv")
stopifnot(file.exists(api_path), file.exists(query_path))

api <- fread(api_path)
queries <- fread(query_path)
setorder(queries, finngen_phenocode, query_chr, query_pos, query_ref, query_alt)
endpoint_codes <- unique(queries$finngen_phenocode)
endpoint_pick <- unique(round(seq(1, length(endpoint_codes), length.out = 12)))
selected_codes <- endpoint_codes[endpoint_pick]
selected_fractions <- seq(0.10, 0.90, length.out = length(selected_codes))

joined <- merge(
  queries,
  api[, .(query_id, api_status, beta_alt_api = beta_alt,
          pval_api = pval, se_api = sebeta_reconstructed)],
  by = "query_id",
  all.x = TRUE
)
sample_rows <- rbindlist(lapply(seq_along(selected_codes), function(i) {
  code <- selected_codes[i]
  x <- joined[finngen_phenocode == code & api_status == "ok"]
  if (!nrow(x)) return(NULL)
  setorder(x, query_chr, query_pos, query_ref, query_alt)
  x[max(1L, ceiling(.N * selected_fractions[i]))]
}))
if (nrow(sample_rows) < 10L) stop("Fewer than 10 deterministic precision-audit queries are usable")

cache_dir <- file.path(root, "metadata/r13_tbi_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
sample_rows[, summary_url := sprintf(
  "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_%s.gz",
  finngen_phenocode
)]
sample_rows[, index_url := paste0(summary_url, ".tbi")]
sample_rows[, index_path := file.path(cache_dir, sprintf("finngen_R13_%s.gz.tbi", finngen_phenocode))]

need <- sample_rows[!file.exists(index_path), .(index_url, index_path)]
if (nrow(need)) {
  multi_download(need$index_url, need$index_path, resume = TRUE, progress = FALSE)
}
if (!all(file.exists(sample_rows$index_path))) stop("Not all precision-audit Tabix indexes were downloaded")

tasks <- split(as.data.frame(sample_rows), seq_len(nrow(sample_rows)))
workers <- min(6L, nrow(sample_rows))
cl <- makeCluster(workers)
on.exit(stopCluster(cl), add = TRUE)
clusterEvalQ(cl, suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicRanges)
}))

extract_one <- function(x) {
  gr <- GRanges(as.character(x$query_chr), IRanges(as.integer(x$query_pos), as.integer(x$query_pos)))
  tf <- TabixFile(x$summary_url, index = x$index_path)
  records <- NULL
  err <- ""
  for (attempt in seq_len(3L)) {
    records <- tryCatch(scanTabix(tf, param = gr)[[1]], error = function(e) {
      err <<- conditionMessage(e)
      NULL
    })
    if (!is.null(records)) break
    Sys.sleep(attempt * 2)
  }
  if (is.null(records) || !length(records)) {
    return(data.frame(query_id = x$query_id, tabix_status = "no_record",
                      beta_tabix = NA_real_, pval_tabix = NA_real_, se_tabix = NA_real_,
                      tabix_error = err, stringsAsFactors = FALSE))
  }
  fields <- strsplit(records, "\t", fixed = TRUE)
  keep <- vapply(fields, function(z) {
    length(z) >= 11L && z[[1L]] == as.character(x$query_chr) &&
      suppressWarnings(as.integer(z[[2L]])) == as.integer(x$query_pos) &&
      toupper(z[[3L]]) == toupper(x$query_ref) && toupper(z[[4L]]) == toupper(x$query_alt)
  }, logical(1))
  hit <- fields[keep]
  if (length(hit) != 1L) {
    return(data.frame(query_id = x$query_id, tabix_status = paste0("exact_hits_", length(hit)),
                      beta_tabix = NA_real_, pval_tabix = NA_real_, se_tabix = NA_real_,
                      tabix_error = err, stringsAsFactors = FALSE))
  }
  z <- hit[[1L]]
  data.frame(
    query_id = x$query_id,
    tabix_status = "ok",
    beta_tabix = as.numeric(z[[9L]]),
    pval_tabix = as.numeric(z[[7L]]),
    se_tabix = as.numeric(z[[10L]]),
    tabix_error = "",
    stringsAsFactors = FALSE
  )
}

tabix_rows <- rbindlist(parLapplyLB(cl, tasks, extract_one), fill = TRUE)
audit <- merge(sample_rows, tabix_rows, by = "query_id", all.x = TRUE)
audit[, beta_abs_error := abs(beta_alt_api - beta_tabix)]
audit[, p_abs_error := abs(pval_api - pval_tabix)]
audit[, se_abs_error := abs(se_api - se_tabix)]
audit[, se_relative_error := se_abs_error / se_tabix]
audit[, row_pass := tabix_status == "ok" &
        is.finite(beta_abs_error) & beta_abs_error <= 1e-7 &
        is.finite(se_relative_error) & se_relative_error <= 1e-3]
setorder(audit, finngen_phenocode)
fwrite(audit, file.path(root, "qa/R13_API_TABIX_PRECISION_AUDIT.tsv"), sep = "\t", na = "NA")

gate <- data.table(
  check = c("deterministic_sample_size", "exact_tabix_records", "beta_precision", "se_precision"),
  observed = c(
    as.character(nrow(audit)),
    sprintf("%d/%d", sum(audit$tabix_status == "ok"), nrow(audit)),
    sprintf("max abs error %.6g", max(audit$beta_abs_error, na.rm = TRUE)),
    sprintf("max relative error %.6g", max(audit$se_relative_error, na.rm = TRUE))
  ),
  required = c(">=10", "all", "<=1e-7", "<=1e-3"),
  pass = c(
    nrow(audit) >= 10L,
    all(audit$tabix_status == "ok"),
    all(is.finite(audit$beta_abs_error) & audit$beta_abs_error <= 1e-7),
    all(is.finite(audit$se_relative_error) & audit$se_relative_error <= 1e-3)
  )
)
fwrite(gate, file.path(root, "qa/R13_API_TABIX_PRECISION_GATE.tsv"), sep = "\t", na = "NA")
writeLines(capture.output(sessionInfo()), file.path(root, "qa/sessionInfo_r13_precision_audit.txt"))
if (!all(gate$pass)) stop("API-to-Tabix precision gate failed")
message("API-to-Tabix precision gate passed for ", nrow(audit), " deterministic queries")
