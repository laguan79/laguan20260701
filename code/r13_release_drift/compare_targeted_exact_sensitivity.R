#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
primary_cells <- fread(
  file.path(root, "outputs/RELEASE_STABILITY_CELL_RESULTS.tsv")
)
exact_cells <- fread(
  file.path(
    root,
    "outputs/RELEASE_STABILITY_CELL_RESULTS_TARGETED_EXACT.tsv"
  )
)
primary_metrics <- fread(file.path(root, "outputs/PATTERN_ABC_METRICS.tsv"))
exact_metrics <- fread(
  file.path(root, "outputs/PATTERN_ABC_METRICS_TARGETED_EXACT.tsv")
)
primary_patterns <- fread(
  file.path(root, "outputs/PATTERN_ABC_GATE_SUMMARY.tsv")
)
exact_patterns <- fread(
  file.path(root, "outputs/PATTERN_ABC_GATE_SUMMARY_TARGETED_EXACT.tsv")
)
locked <- fread(
  file.path(
    root,
    "outputs/R13_LOCKED_VARIANT_EXTRACTION_TARGETED_EXACT_SENSITIVITY.tsv"
  )
)
stopifnot(
  nrow(primary_cells) == 1060L,
  nrow(exact_cells) == 1060L,
  nrow(locked) == 42200L
)

setorder(primary_cells, outcome_family, trait_id)
setorder(exact_cells, outcome_family, trait_id)
stopifnot(
  identical(primary_cells$trait_id, exact_cells$trait_id),
  identical(primary_cells$outcome_family, exact_cells$outcome_family)
)

call_fields <- c(
  "source_fdr_r13",
  "source_equivalent_r13",
  "naive_asymmetry_r13",
  "class_r13",
  "strict_axis_migration",
  "full_class_migration",
  "robust_large_release_drift",
  "U_select",
  "unsupported_inflation_cell"
)
call_comparison <- rbindlist(lapply(call_fields, function(field) {
  a <- primary_cells[[field]]
  b <- exact_cells[[field]]
  data.table(
    field = field,
    primary_count = if (is.logical(a) || is.numeric(a)) {
      sum(a, na.rm = TRUE)
    } else {
      NA_real_
    },
    exact_sensitivity_count = if (is.logical(b) || is.numeric(b)) {
      sum(b, na.rm = TRUE)
    } else {
      NA_real_
    },
    changed_cells = sum(a != b, na.rm = TRUE)
  )
}))

class_transition <- exact_cells[, .(
  trait_id,
  outcome_family,
  primary_class_r13 = primary_cells$class_r13,
  exact_class_r13 = class_r13
)][primary_class_r13 != exact_class_r13]

metric_comparison <- merge(
  primary_metrics[, .(
    metric,
    primary_estimate = estimate,
    primary_ci_low = ci_low,
    primary_ci_high = ci_high
  )],
  exact_metrics[, .(
    metric,
    exact_estimate = estimate,
    exact_ci_low = ci_low,
    exact_ci_high = ci_high
  )],
  by = "metric",
  all = TRUE
)
metric_comparison[, estimate_change := exact_estimate - primary_estimate]

pattern_comparison <- merge(
  primary_patterns[, .(
    pattern,
    primary_conditions_passed = conditions_passed,
    primary_conditions_total = conditions_total,
    primary_pattern_pass = pattern_pass
  )],
  exact_patterns[, .(
    pattern,
    exact_conditions_passed = conditions_passed,
    exact_conditions_total = conditions_total,
    exact_pattern_pass = pattern_pass
  )],
  by = "pattern",
  all = TRUE
)
pattern_comparison[, decision_changed :=
  primary_pattern_pass != exact_pattern_pass]

critical_cells <- unique(rbindlist(list(
  primary_cells[original_common_heterogeneity == TRUE, .(
    trait_id, outcome_family, reason = "historical35"
  )],
  primary_cells[strict_axis_migration == TRUE, .(
    trait_id, outcome_family, reason = "axis_migration"
  )],
  primary_cells[robust_large_release_drift == TRUE, .(
    trait_id, outcome_family, reason = "robust_drift"
  )],
  primary_cells[unsupported_inflation_cell > 0, .(
    trait_id, outcome_family, reason = "new_unsupported"
  )]
)))
coverage <- merge(
  locked[, .(
    snps = .N,
    exact_snps = sum(exact_se_replaced)
  ), by = .(trait_id, outcome_family)],
  critical_cells,
  by = c("trait_id", "outcome_family"),
  all.y = TRUE,
  allow.cartesian = TRUE
)
coverage[, full_exact_coverage := exact_snps == snps]
coverage_summary <- coverage[, .(
  cells = .N,
  fully_exact_cells = sum(full_exact_coverage),
  all_snp_rows = sum(snps),
  exact_snp_rows = sum(exact_snps),
  full_coverage = all(full_exact_coverage)
), by = reason]

fwrite(
  call_comparison,
  file.path(root, "qa/TARGETED_EXACT_CALL_COMPARISON.tsv"),
  sep = "\t",
  na = "NA"
)
fwrite(
  class_transition,
  file.path(root, "qa/TARGETED_EXACT_CLASS_TRANSITIONS.tsv"),
  sep = "\t",
  na = "NA"
)
fwrite(
  metric_comparison,
  file.path(root, "qa/TARGETED_EXACT_METRIC_COMPARISON.tsv"),
  sep = "\t",
  na = "NA"
)
fwrite(
  pattern_comparison,
  file.path(root, "qa/TARGETED_EXACT_PATTERN_COMPARISON.tsv"),
  sep = "\t",
  na = "NA"
)
fwrite(
  coverage_summary,
  file.path(root, "qa/TARGETED_EXACT_CRITICAL_CELL_COVERAGE.tsv"),
  sep = "\t",
  na = "NA"
)

gate <- data.table(
  check = c(
    "cell_contract",
    "critical_cell_exact_coverage",
    "pattern_decision_stability",
    "no_pattern_passes_after_exact_se",
    "finite_metric_comparison"
  ),
  observed = c(
    as.character(nrow(exact_cells)),
    sprintf(
      "%d/%d reason groups",
      sum(coverage_summary$full_coverage),
      nrow(coverage_summary)
    ),
    sprintf(
      "%d/%d patterns unchanged",
      sum(!pattern_comparison$decision_changed),
      nrow(pattern_comparison)
    ),
    paste(
      pattern_comparison[exact_pattern_pass == TRUE, pattern],
      collapse = ","
    ),
    sprintf(
      "%d/%d",
      sum(is.finite(metric_comparison$exact_estimate)),
      nrow(metric_comparison)
    )
  ),
  required = c("1060", "all", "all", "none", "all"),
  pass = c(
    nrow(exact_cells) == 1060L,
    all(coverage_summary$full_coverage),
    all(!pattern_comparison$decision_changed),
    !any(pattern_comparison$exact_pattern_pass),
    all(is.finite(metric_comparison$exact_estimate))
  )
)
fwrite(
  gate,
  file.path(root, "qa/TARGETED_EXACT_SENSITIVITY_GATE.tsv"),
  sep = "\t",
  na = "NA"
)
writeLines(
  capture.output(sessionInfo()),
  file.path(root, "qa/sessionInfo_targeted_exact_comparison.txt")
)
if (!all(gate$pass)) stop("Targeted exact-SE sensitivity gate failed")
message("Targeted exact-SE sensitivity comparison passed")
