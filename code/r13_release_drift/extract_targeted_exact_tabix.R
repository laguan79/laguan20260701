#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(parallel)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
lock_path <- file.path(root, "inputs/R13_TARGETED_EXACT_TABIX_QUERY_LOCK.tsv")
hash_path <- file.path(root, "inputs/R13_TARGETED_EXACT_TABIX_QUERY_LOCK.sha256")
stopifnot(file.exists(lock_path), file.exists(hash_path))
expected_hash <- strsplit(readLines(hash_path, warn = FALSE)[1], " ", fixed = TRUE)[[1]][1]
observed_hash <- digest(lock_path, algo = "sha256", file = TRUE)
if (!identical(tolower(expected_hash), tolower(observed_hash))) stop("Targeted Tabix query lock hash mismatch")

targets <- fread(lock_path)
cache_dir <- file.path(root, "metadata/r13_tbi_cache")
batch_dir <- file.path(root, "outputs/r13_targeted_exact_batches")
dir.create(batch_dir, recursive = TRUE, showWarnings = FALSE)
targets[, summary_url := sprintf(
  "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_%s.gz",
  finngen_phenocode
)]
targets[, index_path := file.path(cache_dir, sprintf("finngen_R13_%s.gz.tbi", finngen_phenocode))]
if (!all(file.exists(targets$index_path))) stop("Local R13 Tabix index cache is incomplete")

extract_one <- function(x) {
  suppressPackageStartupMessages({
    library(Rsamtools)
    library(GenomicRanges)
  })
  gr <- GRanges(as.character(x$query_chr), IRanges(as.integer(x$query_pos), as.integer(x$query_pos)))
  tf <- TabixFile(x$summary_url, index = x$index_path)
  records <- NULL
  err <- ""
  for (attempt in seq_len(4L)) {
    records <- tryCatch(scanTabix(tf, param = gr)[[1]], error = function(e) {
      err <<- conditionMessage(e)
      NULL
    })
    if (!is.null(records)) break
    Sys.sleep(2 * attempt)
  }
  base <- data.frame(query_id = x$query_id, stringsAsFactors = FALSE)
  if (is.null(records) || !length(records)) {
    return(cbind(base, data.frame(tabix_status = "no_record", beta_tabix = NA_real_,
                                  pval_tabix = NA_real_, se_tabix = NA_real_, error = err)))
  }
  fields <- strsplit(records, "\t", fixed = TRUE)
  keep <- vapply(fields, function(z) {
    length(z) >= 11L && z[[1L]] == as.character(x$query_chr) &&
      as.integer(z[[2L]]) == as.integer(x$query_pos) &&
      toupper(z[[3L]]) == toupper(x$query_ref) && toupper(z[[4L]]) == toupper(x$query_alt)
  }, logical(1))
  hit <- fields[keep]
  if (length(hit) != 1L) {
    return(cbind(base, data.frame(tabix_status = paste0("exact_hits_", length(hit)),
                                  beta_tabix = NA_real_, pval_tabix = NA_real_,
                                  se_tabix = NA_real_, error = err)))
  }
  z <- hit[[1L]]
  cbind(base, data.frame(
    tabix_status = "ok",
    beta_tabix = as.numeric(z[[9L]]),
    pval_tabix = as.numeric(z[[7L]]),
    se_tabix = as.numeric(z[[10L]]),
    error = "",
    stringsAsFactors = FALSE
  ))
}

batch_size <- 400L
targets[, batch_id := sprintf("batch_%03d", ceiling(seq_len(.N) / batch_size))]
batch_ids <- unique(targets$batch_id)
workers <- min(32L, nrow(targets))
cl <- makeCluster(workers)
on.exit(stopCluster(cl), add = TRUE)

for (bid in batch_ids) {
  path <- file.path(batch_dir, paste0(bid, ".tsv"))
  q <- targets[batch_id == bid]
  completed <- NULL
  if (file.exists(path)) {
    old <- tryCatch(fread(path), error = function(e) NULL)
    if (!is.null(old) && nrow(old) == nrow(q) && setequal(old$query_id, q$query_id) &&
        all(old$tabix_status == "ok")) {
      message(bid, " cached")
      next
    }
    if (!is.null(old) && all(c(
      "query_id", "tabix_status", "beta_tabix", "pval_tabix", "se_tabix", "error"
    ) %chin% names(old))) {
      completed <- old[
        tabix_status == "ok",
        .(query_id, tabix_status, beta_tabix, pval_tabix, se_tabix, error)
      ]
      completed <- unique(completed, by = "query_id")
    }
  }
  pending <- if (is.null(completed)) q else q[!query_id %chin% completed$query_id]
  message(
    bid, ": exact Tabix audit for ", nrow(pending), " unresolved of ",
    nrow(q), " queries using ", workers, " workers"
  )
  new_rows <- rbindlist(
    parLapplyLB(
      cl,
      split(as.data.frame(pending), seq_len(nrow(pending))),
      extract_one
    ),
    fill = TRUE
  )
  result_rows <- rbindlist(list(completed, new_rows), use.names = TRUE, fill = TRUE)
  result_rows <- unique(result_rows, by = "query_id")
  rows <- merge(q, result_rows, by = "query_id", all.x = TRUE)
  fwrite(rows, path, sep = "\t", na = "NA")
  if (!all(rows$tabix_status == "ok")) stop("Unresolved exact Tabix rows in ", bid)
}

files <- file.path(batch_dir, paste0(batch_ids, ".tsv"))
exact <- rbindlist(lapply(files, fread), use.names = TRUE, fill = TRUE)
setorder(exact, finngen_phenocode, query_chr, query_pos)
stopifnot(nrow(exact) == nrow(targets), uniqueN(exact$query_id) == nrow(exact))
exact[, `:=`(
  beta_abs_error = abs(beta_alt_api - beta_tabix),
  p_abs_error = abs(pval_api - pval_tabix),
  se_abs_error = abs(se_api - se_tabix),
  se_relative_error = abs(se_api - se_tabix) / se_tabix
)]
fwrite(exact, file.path(root, "outputs/R13_TARGETED_EXACT_TABIX.tsv"), sep = "\t", na = "NA")

gate <- data.table(
  check = c("locked_query_hash", "target_query_count", "exact_record_count",
            "beta_match", "p_match", "finite_exact_se"),
  observed = c(
    observed_hash,
    as.character(nrow(exact)),
    sprintf("%d/%d", sum(exact$tabix_status == "ok"), nrow(exact)),
    sprintf("max abs %.6g", max(exact$beta_abs_error)),
    sprintf("max abs %.6g", max(exact$p_abs_error)),
    sprintf("%d/%d", sum(is.finite(exact$se_tabix) & exact$se_tabix > 0), nrow(exact))
  ),
  required = c(expected_hash, as.character(nrow(targets)), "all", "<=1e-7", "<=1e-12", "all"),
  pass = c(
    identical(tolower(expected_hash), tolower(observed_hash)),
    nrow(exact) == nrow(targets),
    all(exact$tabix_status == "ok"),
    max(exact$beta_abs_error) <= 1e-7,
    max(exact$p_abs_error) <= 1e-12,
    all(is.finite(exact$se_tabix) & exact$se_tabix > 0)
  )
)
fwrite(gate, file.path(root, "qa/R13_TARGETED_EXACT_TABIX_GATE.tsv"), sep = "\t", na = "NA")
writeLines(capture.output(sessionInfo()), file.path(root, "qa/sessionInfo_targeted_exact_tabix.txt"))
if (!all(gate$pass)) stop("Targeted exact Tabix gate failed")
message("Targeted exact Tabix gate passed for ", nrow(exact), " queries")
