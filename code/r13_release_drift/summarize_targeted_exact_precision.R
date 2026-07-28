#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
path <- file.path(root, "outputs/R13_TARGETED_EXACT_TABIX.tsv")
stopifnot(file.exists(path))
d <- fread(path)
stopifnot(nrow(d) == 3517L, all(d$tabix_status == "ok"))

summarize_rows <- function(x, label) {
  data.table(
    stratum = label,
    rows = nrow(x),
    median_se_relative_error = median(x$se_relative_error, na.rm = TRUE),
    p95_se_relative_error = quantile(
      x$se_relative_error, 0.95, na.rm = TRUE, names = FALSE
    ),
    p99_se_relative_error = quantile(
      x$se_relative_error, 0.99, na.rm = TRUE, names = FALSE
    ),
    max_se_relative_error = max(x$se_relative_error, na.rm = TRUE),
    rows_over_1e3 = sum(x$se_relative_error > 1e-3, na.rm = TRUE),
    rows_over_1pct = sum(x$se_relative_error > 0.01, na.rm = TRUE),
    max_beta_abs_error = max(x$beta_abs_error, na.rm = TRUE),
    max_p_abs_error = max(x$p_abs_error, na.rm = TRUE)
  )
}

summary <- summarize_rows(d, "all_targeted_exact_rows")
strata <- rbindlist(list(
  summarize_rows(d[audit_historical35 == TRUE], "historical35"),
  summarize_rows(d[audit_axis_migration == TRUE], "axis_migration"),
  summarize_rows(d[audit_robust_drift == TRUE], "robust_drift"),
  summarize_rows(d[audit_new_unsupported == TRUE], "new_unsupported"),
  summarize_rows(
    d[audit_extreme_near_null == TRUE],
    "extreme_near_null"
  ),
  summarize_rows(
    d[audit_endpoint_p_extreme == TRUE],
    "endpoint_p_extreme"
  )
), use.names = TRUE)

top <- copy(d)
setorder(top, -se_relative_error)
top <- top[seq_len(min(100L, .N))]

fwrite(
  summary,
  file.path(root, "qa/R13_TARGETED_EXACT_PRECISION_SUMMARY.tsv"),
  sep = "\t",
  na = "NA"
)
fwrite(
  strata,
  file.path(root, "qa/R13_TARGETED_EXACT_PRECISION_BY_REASON.tsv"),
  sep = "\t",
  na = "NA"
)
fwrite(
  top,
  file.path(root, "qa/R13_TARGETED_EXACT_PRECISION_TOP_ERRORS.tsv"),
  sep = "\t",
  na = "NA"
)
writeLines(
  capture.output(sessionInfo()),
  file.path(root, "qa/sessionInfo_targeted_exact_precision.txt")
)
message("Targeted exact precision summary written")
