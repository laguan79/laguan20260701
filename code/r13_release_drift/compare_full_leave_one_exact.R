#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
primary <- fread(
  file.path(root, "outputs/FULL_LEAVE_ONE_DECISION_PIPELINE.tsv")
)
exact <- fread(
  file.path(
    root,
    "outputs/FULL_LEAVE_ONE_DECISION_PIPELINE_TARGETED_EXACT.tsv"
  )
)
keys <- c("leave_one_type", "leave_one_id")
setkeyv(primary, keys)
setkeyv(exact, keys)
stopifnot(nrow(primary) == 74L, nrow(exact) == 74L)

metrics <- setdiff(intersect(names(primary), names(exact)), keys)
comparison <- rbindlist(lapply(metrics, function(metric) {
  a <- primary[[metric]]
  b <- exact[[metric]]
  data.table(
    metric = metric,
    max_abs_change = max(abs(a - b), na.rm = TRUE),
    changed_omissions = sum(a != b, na.rm = TRUE)
  )
}))
count_metrics <- c(
  "cells",
  "original_calls_remaining",
  "original_call_retained_count",
  "r12_joint_fdr_count",
  "r13_joint_fdr_count"
)
gate <- data.table(
  check = c("row_contract", "count_decision_stability"),
  observed = c(
    sprintf("%d/%d", nrow(primary), nrow(exact)),
    sprintf(
      "%d changed metric-omission combinations",
      comparison[metric %chin% count_metrics, sum(changed_omissions)]
    )
  ),
  required = c("74/74", "0"),
  pass = c(
    nrow(primary) == 74L && nrow(exact) == 74L,
    comparison[
      metric %chin% count_metrics,
      sum(changed_omissions)
    ] == 0L
  )
)
fwrite(
  comparison,
  file.path(root, "qa/FULL_LEAVE_ONE_TARGETED_EXACT_COMPARISON.tsv"),
  sep = "\t",
  na = "NA"
)
fwrite(
  gate,
  file.path(root, "qa/FULL_LEAVE_ONE_TARGETED_EXACT_GATE.tsv"),
  sep = "\t",
  na = "NA"
)
writeLines(
  capture.output(sessionInfo()),
  file.path(root, "qa/sessionInfo_full_leave_one_exact_comparison.txt")
)
if (!all(gate$pass)) stop("Full leave-one exact-SE comparison gate failed")
message("Full leave-one exact-SE comparison passed")
