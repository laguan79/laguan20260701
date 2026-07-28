#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(httr2)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
required <- c(
  "inputs/LOCKED_COMMON_SNP_IDENTITIES.tsv",
  "inputs/LOCKED_COMMON_SNP_PAIR_COUNTS.tsv",
  "qa/PRE_EFFECT_CONTRACT_GATE.tsv"
)
stopifnot(all(file.exists(file.path(root, required))))

gate <- fread(file.path(root, "qa/PRE_EFFECT_CONTRACT_GATE.tsv"))
if (!all(gate$pass)) stop("Pre-effect contract gate is not fully passed")

out_dir <- file.path(root, "outputs")
batch_dir <- file.path(out_dir, "r13_api_batches")
qa_dir <- file.path(root, "qa")
dir.create(batch_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

locked <- fread(file.path(root, "inputs/LOCKED_COMMON_SNP_IDENTITIES.tsv"))
pair_counts <- fread(file.path(root, "inputs/LOCKED_COMMON_SNP_PAIR_COUNTS.tsv"))
stopifnot(nrow(locked) == 42200L, uniqueN(locked$trait_id) == 20L,
          uniqueN(locked$outcome_family) == 53L)

queries <- unique(locked[, .(
  finngen_phenocode,
  variant_key = r12_selected_variant_key,
  query_chr,
  query_pos,
  query_ref = toupper(query_ref),
  query_alt = toupper(query_alt)
)])
setorder(queries, finngen_phenocode, query_chr, query_pos, query_ref, query_alt)
queries[, query_id := .I]
queries[, api_url := sprintf(
  "https://r13.finngen.fi/api/variant/%s/%s",
  variant_key, finngen_phenocode
)]
stopifnot(nrow(queries) == 38244L, uniqueN(queries$query_id) == nrow(queries))
fwrite(queries, file.path(out_dir, "R13_API_QUERY_INDEX.tsv"), sep = "\t", na = "NA")

scalar_num <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_real_)
  suppressWarnings(as.numeric(x[[1L]]))
}

parse_response <- function(resp, qrow) {
  base <- qrow[, .(query_id, finngen_phenocode, variant_key, query_chr,
                   query_pos, query_ref, query_alt, api_url)]
  if (inherits(resp, "httr2_failure")) {
    return(base[, `:=`(http_status = NA_integer_, api_status = "transport_error",
                       beta_alt = NA_real_, pval = NA_real_, mlogp = NA_real_,
                       sebeta_reconstructed = NA_real_, se_source = NA_character_,
                       maf_case = NA_real_, maf_control = NA_real_, error = conditionMessage(resp))])
  }
  code <- resp_status(resp)
  if (code != 200L) {
    status <- if (code == 404L) "http_404" else paste0("http_", code)
    return(base[, `:=`(http_status = code, api_status = status,
                       beta_alt = NA_real_, pval = NA_real_, mlogp = NA_real_,
                       sebeta_reconstructed = NA_real_, se_source = NA_character_,
                       maf_case = NA_real_, maf_control = NA_real_, error = "")])
  }

  body <- tryCatch(resp_body_json(resp, simplifyVector = TRUE), error = identity)
  if (inherits(body, "error") || is.null(body$results)) {
    msg <- if (inherits(body, "error")) conditionMessage(body) else "missing results"
    return(base[, `:=`(http_status = code, api_status = "parse_error",
                       beta_alt = NA_real_, pval = NA_real_, mlogp = NA_real_,
                       sebeta_reconstructed = NA_real_, se_source = NA_character_,
                       maf_case = NA_real_, maf_control = NA_real_, error = msg)])
  }

  ans <- body$results
  beta <- scalar_num(ans$beta)
  pval <- scalar_num(ans$pval)
  mlogp <- scalar_num(ans$mlogp)
  maf_case <- scalar_num(ans$maf_case)
  maf_control <- scalar_num(ans$maf_control)

  z_abs <- NA_real_
  se_source <- NA_character_
  if (is.finite(pval) && pval > 0 && pval <= 1) {
    z_abs <- qnorm(pval / 2, lower.tail = FALSE)
    se_source <- "beta_two_sided_p"
  } else if (is.finite(mlogp) && mlogp > 0) {
    log_half_p <- -mlogp * log(10) - log(2)
    z_abs <- qnorm(log_half_p, lower.tail = FALSE, log.p = TRUE)
    se_source <- "beta_mlogp"
  }
  se <- if (is.finite(beta) && beta != 0 && is.finite(z_abs) && z_abs > 0) {
    abs(beta) / z_abs
  } else {
    NA_real_
  }
  status <- if (is.finite(beta) && is.finite(se) && se > 0) "ok" else "no_usable_effect"

  base[, `:=`(
    http_status = code,
    api_status = status,
    beta_alt = beta,
    pval = pval,
    mlogp = mlogp,
    sebeta_reconstructed = se,
    se_source = se_source,
    maf_case = maf_case,
    maf_control = maf_control,
    error = ""
  )]
}

make_request <- function(url) {
  request(url) |>
    req_user_agent("source-release-stability-r13/1.0") |>
    req_timeout(35) |>
    req_retry(max_tries = 3, backoff = function(i) min(2^(i - 1), 8),
              is_transient = function(resp) resp_status(resp) %in% c(408, 425, 429, 500, 502, 503, 504))
}

batch_size <- 400L
max_active <- 40L
queries[, batch_id := sprintf("batch_%03d", ceiling(query_id / batch_size))]
batch_ids <- unique(queries$batch_id)

for (bid in batch_ids) {
  batch_path <- file.path(batch_dir, paste0(bid, ".tsv"))
  q <- queries[batch_id == bid]
  if (file.exists(batch_path)) {
    old <- tryCatch(fread(batch_path), error = function(e) NULL)
    if (!is.null(old) && nrow(old) == nrow(q) &&
        setequal(old$query_id, q$query_id) &&
        !any(old$api_status %chin% c("transport_error", "parse_error") |
             grepl("^http_(408|425|429|500|502|503|504)$", old$api_status))) {
      message(bid, " cached")
      next
    }
  }

  message(bid, ": querying ", nrow(q), " locked endpoint-variant pairs")
  reqs <- lapply(q$api_url, make_request)
  resps <- req_perform_parallel(
    reqs,
    on_error = "continue",
    progress = FALSE,
    max_active = max_active
  )
  rows <- rbindlist(lapply(seq_len(nrow(q)), function(i) parse_response(resps[[i]], q[i])))
  stopifnot(nrow(rows) == nrow(q), identical(rows$query_id, q$query_id))
  fwrite(rows, batch_path, sep = "\t", na = "NA")
}

batch_files <- file.path(batch_dir, paste0(batch_ids, ".tsv"))
stopifnot(all(file.exists(batch_files)))
api <- rbindlist(lapply(batch_files, fread), use.names = TRUE, fill = TRUE)
setorder(api, query_id)
stopifnot(nrow(api) == nrow(queries), identical(api$query_id, queries$query_id))

retryable <- api$api_status %chin% c("transport_error", "parse_error") |
  grepl("^http_(408|425|429|500|502|503|504)$", api$api_status)
if (any(retryable)) {
  fwrite(api[retryable], file.path(qa_dir, "R13_API_UNRESOLVED_QUERIES.tsv"), sep = "\t", na = "NA")
  stop(sum(retryable), " retryable R13 API queries remain unresolved")
}

fwrite(api, file.path(out_dir, "R13_API_UNIQUE_VARIANT_EFFECTS.tsv"), sep = "\t", na = "NA")

api_join <- api[, .(
  finngen_phenocode, variant_key, api_status, http_status,
  beta_alt, pval_r13 = pval, mlogp_r13 = mlogp,
  se_r13 = sebeta_reconstructed, se_source,
  maf_case_r13 = maf_case, maf_control_r13 = maf_control,
  api_error = error
)]
extracted <- merge(
  locked,
  api_join,
  by.x = c("finngen_phenocode", "r12_selected_variant_key"),
  by.y = c("finngen_phenocode", "variant_key"),
  all.x = TRUE,
  sort = FALSE
)

extracted[, allele_orientation := fifelse(
  toupper(effect_allele_exposure) == toupper(query_alt) &
    toupper(other_allele_exposure) == toupper(query_ref),
  "exposure_effect_is_alt",
  fifelse(
    toupper(effect_allele_exposure) == toupper(query_ref) &
      toupper(other_allele_exposure) == toupper(query_alt),
    "exposure_effect_is_ref",
    "allele_mismatch"
  )
)]
extracted[, beta_r13_aligned := fifelse(
  allele_orientation == "exposure_effect_is_alt", beta_alt,
  fifelse(allele_orientation == "exposure_effect_is_ref", -beta_alt, NA_real_)
)]
extracted[, extraction_status := fcase(
  is.na(api_status), "join_failure",
  api_status != "ok", api_status,
  allele_orientation == "allele_mismatch", "allele_mismatch",
  !is.finite(beta_r13_aligned), "invalid_beta",
  !is.finite(se_r13) | se_r13 <= 0, "invalid_se",
  default = "retained"
)]
setorder(extracted, outcome_family, trait_id, common_snp_order)
stopifnot(nrow(extracted) == 42200L)
fwrite(extracted, file.path(out_dir, "R13_LOCKED_VARIANT_EXTRACTION.tsv"), sep = "\t", na = "NA")

pair_qa <- extracted[, .(
  n_common_snps = .N,
  n_r13_available = sum(extraction_status == "retained"),
  n_missing = sum(extraction_status != "retained"),
  allele_mismatches = sum(extraction_status == "allele_mismatch"),
  invalid_beta_se = sum(extraction_status %chin% c("invalid_beta", "invalid_se")),
  api_absent_or_unusable = sum(extraction_status %chin% c("http_404", "no_usable_effect")),
  transport_or_join_errors = sum(extraction_status %chin% c("transport_error", "parse_error", "join_failure"))
), by = .(trait_id, outcome_family, finngen_phenocode)]
pair_qa <- merge(pair_qa, pair_counts, by = c("trait_id", "outcome_family", "finngen_phenocode"),
                 suffixes = c("", "_frozen"), all.x = TRUE)
pair_qa[, retention_fraction := n_r13_available / n_common_snps]
pair_qa[, estimable := n_r13_available >= 3L & n_r13_available >= min_required_r13_snps]
setorder(pair_qa, outcome_family, trait_id)
stopifnot(nrow(pair_qa) == 1060L)
fwrite(pair_qa, file.path(out_dir, "THREE_WAY_INSTRUMENT_QA.tsv"), sep = "\t", na = "NA")

status_counts <- extracted[, .N, by = extraction_status][order(extraction_status)]
fwrite(status_counts, file.path(qa_dir, "R13_EXTRACTION_STATUS_COUNTS.tsv"), sep = "\t", na = "NA")

n_estimable <- sum(pair_qa$estimable)
gate_out <- data.table(
  check = c(
    "locked_query_contract",
    "api_transport_complete",
    "pair_cells_present",
    "estimable_cells",
    "retained_allele_qa",
    "retained_scale_qa"
  ),
  observed = c(
    sprintf("%d unique; %d expanded", nrow(api), nrow(extracted)),
    sprintf("%d unresolved", sum(retryable)),
    as.character(nrow(pair_qa)),
    sprintf("%d/1060 (%.2f%%)", n_estimable, 100 * n_estimable / 1060),
    sprintf("%d mismatches among retained", sum(extracted$extraction_status == "retained" & extracted$allele_orientation == "allele_mismatch")),
    sprintf("%d invalid among retained", sum(extracted$extraction_status == "retained" & (!is.finite(extracted$beta_r13_aligned) | !is.finite(extracted$se_r13) | extracted$se_r13 <= 0)))
  ),
  required = c(
    "38244 unique; 42200 expanded",
    "0 unresolved",
    "1060",
    ">=1007/1060",
    "0",
    "0"
  ),
  pass = c(
    nrow(api) == 38244L && nrow(extracted) == 42200L,
    sum(retryable) == 0L,
    nrow(pair_qa) == 1060L,
    n_estimable >= ceiling(0.95 * 1060),
    !any(extracted$extraction_status == "retained" & extracted$allele_orientation == "allele_mismatch"),
    !any(extracted$extraction_status == "retained" & (!is.finite(extracted$beta_r13_aligned) | !is.finite(extracted$se_r13) | extracted$se_r13 <= 0))
  )
)
fwrite(gate_out, file.path(qa_dir, "R13_EXTRACTION_GATE.tsv"), sep = "\t", na = "NA")

writeLines(capture.output(sessionInfo()), file.path(qa_dir, "sessionInfo_r13_extraction.txt"))
if (!all(gate_out$pass)) stop("R13 extraction gate failed; downstream effect estimation is blocked")
message("R13 extraction gate passed: ", n_estimable, "/1060 cells estimable")
