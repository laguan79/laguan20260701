#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
analysis_tag <- trimws(Sys.getenv("R13_ANALYSIS_TAG", unset = ""))
tagged_path <- function(path) {
  if (!nzchar(analysis_tag)) return(path)
  sub("(\\.[^.]+)$", paste0("_", analysis_tag, "\\1"), path)
}
estimate_path <- Sys.getenv(
  "R13_ESTIMATE_PATH",
  unset = file.path(root, "outputs/THREE_SOURCE_RELEASE_ESTIMATES.tsv")
)
estimate_gate <- Sys.getenv(
  "R13_ESTIMATE_GATE",
  unset = file.path(root, "qa/THREE_SOURCE_ESTIMATION_GATE.tsv")
)
stopifnot(file.exists(estimate_path), file.exists(estimate_gate))
if (!all(fread(estimate_gate)$pass)) stop("Three-source estimation gate is not fully passed")

d <- fread(estimate_path)
stopifnot(nrow(d) == 1060L, uniqueN(d$trait_id) == 20L, uniqueN(d$outcome_family) == 53L)
if (any(!d$estimable)) stop("Primary analysis requires all frozen cells to remain estimable")

rho_source_grid <- seq(-0.30, 0.50, by = 0.05)
rho_release_grid <- c(0, 0.50, 0.70, 0.80, 0.90, 0.95, 0.975)
margin_source <- log(1.20)
margin_release <- log(1.10)

contrast_stats <- function(b1, se1, b2, se2, rho) {
  delta <- b1 - b2
  se <- sqrt(pmax(se1^2 + se2^2 - 2 * rho * se1 * se2, 0))
  z <- delta / se
  p <- 2 * pnorm(-abs(z))
  data.table(delta = delta, se = se, z = z, p = p)
}

source_rho_rows <- vector("list", length(rho_source_grid))
for (j in seq_along(rho_source_grid)) {
  rho <- rho_source_grid[j]
  a <- contrast_stats(d$b_r12, d$se_r12, d$b_panukb, d$se_panukb, rho)
  b <- contrast_stats(d$b_r13, d$se_r13, d$b_panukb, d$se_panukb, rho)
  q_joint <- p.adjust(c(a$p, b$p), method = "BH")
  q_by_joint <- p.adjust(c(a$p, b$p), method = "BY")
  a[, `:=`(q_bh = q_joint[seq_len(nrow(d))], q_by = q_by_joint[seq_len(nrow(d))])]
  b[, `:=`(q_bh = q_joint[nrow(d) + seq_len(nrow(d))], q_by = q_by_joint[nrow(d) + seq_len(nrow(d))])]
  a[, `:=`(
    ci90_low = delta - qnorm(0.95) * se,
    ci90_high = delta + qnorm(0.95) * se
  )]
  b[, `:=`(
    ci90_low = delta - qnorm(0.95) * se,
    ci90_high = delta + qnorm(0.95) * se
  )]
  a[, equivalent := ci90_low >= -margin_source & ci90_high <= margin_source]
  b[, equivalent := ci90_low >= -margin_source & ci90_high <= margin_source]
  a[, fdr_difference := q_bh < 0.05]
  b[, fdr_difference := q_bh < 0.05]
  a[, release := "r12"]
  b[, release := "r13"]
  x <- rbindlist(list(a, b))
  x[, `:=`(
    trait_id = rep(d$trait_id, 2L),
    outcome_family = rep(d$outcome_family, 2L),
    field_domain = rep(d$field_domain, 2L),
    rho_source = rho
  )]
  source_rho_rows[[j]] <- x
}
source_rho <- rbindlist(source_rho_rows)
setcolorder(source_rho, c("trait_id", "outcome_family", "field_domain", "release", "rho_source"))

release_rho_rows <- vector("list", length(rho_release_grid))
for (j in seq_along(rho_release_grid)) {
  rho <- rho_release_grid[j]
  x <- contrast_stats(d$b_r13, d$se_r13, d$b_r12, d$se_r12, rho)
  x[, `:=`(
    q_bh = p.adjust(p, method = "BH"),
    q_by = p.adjust(p, method = "BY"),
    trait_id = d$trait_id,
    outcome_family = d$outcome_family,
    field_domain = d$field_domain,
    rho_release = rho,
    operationally_large = abs(delta) >= margin_release
  )]
  release_rho_rows[[j]] <- x
}
release_rho <- rbindlist(release_rho_rows)
setcolorder(release_rho, c("trait_id", "outcome_family", "field_domain", "rho_release"))

primary_source <- source_rho[abs(rho_source) < 1e-12]
r12s <- primary_source[release == "r12"]
r13s <- primary_source[release == "r13"]
setorder(r12s, outcome_family, trait_id)
setorder(r13s, outcome_family, trait_id)
setorder(d, outcome_family, trait_id)
stopifnot(identical(r12s$trait_id, d$trait_id), identical(r13s$trait_id, d$trait_id))

assign_class <- function(fdr_difference, equivalent, p_finn, p_pan, b_finn, b_pan) {
  concordant <- p_finn < 0.05 & p_pan < 0.05 & sign(b_finn) == sign(b_pan)
  asymmetric <- xor(p_finn < 0.05, p_pan < 0.05)
  fcase(
    fdr_difference, "evidence_of_source_heterogeneity",
    concordant, "concordant_cross_source_support",
    equivalent, "evidence_of_practical_equivalence",
    asymmetric, "asymmetric_support_without_demonstrated_heterogeneity",
    default = "inconclusive_or_no_evidence_in_either_source"
  )
}

primary <- d[, .(
  trait_id, outcome_family, finngen_phenocode, field_domain,
  n_three_way, b_r12, se_r12, pval_r12, b_r13, se_r13, pval_r13,
  b_panukb, se_panukb, pval_panukb,
  delta_release, delta_source_r12, delta_source_r13,
  original_common_heterogeneity = common_heterogeneity,
  original_common_naive_asymmetry = common_naive_asymmetry,
  original_common_class = common_class
)]
primary[, `:=`(
  source_p_r12 = r12s$p,
  source_q_r12 = r12s$q_bh,
  source_q_by_r12 = r12s$q_by,
  source_equivalent_r12 = r12s$equivalent,
  source_p_r13 = r13s$p,
  source_q_r13 = r13s$q_bh,
  source_q_by_r13 = r13s$q_by,
  source_equivalent_r13 = r13s$equivalent
)]
primary[, `:=`(
  source_fdr_r12 = source_q_r12 < 0.05,
  source_fdr_r13 = source_q_r13 < 0.05,
  naive_asymmetry_r12 = xor(pval_r12 < 0.05, pval_panukb < 0.05),
  naive_asymmetry_r13 = xor(pval_r13 < 0.05, pval_panukb < 0.05)
)]
primary[, class_r12 := assign_class(source_fdr_r12, source_equivalent_r12,
                                    pval_r12, pval_panukb, b_r12, b_panukb)]
primary[, class_r13 := assign_class(source_fdr_r13, source_equivalent_r13,
                                    pval_r13, pval_panukb, b_r13, b_panukb)]
primary[, `:=`(
  strict_axis_migration = (source_fdr_r12 != source_fdr_r13) |
    (source_equivalent_r12 != source_equivalent_r13),
  full_class_migration = class_r12 != class_r13,
  axis_agreement = (source_fdr_r12 == source_fdr_r13) &
    (source_equivalent_r12 == source_equivalent_r13),
  positive_fixed_any = (pval_r12 < 0.05) | (pval_panukb < 0.05),
  positive_select_any = (pval_r12 < 0.05) | (pval_r13 < 0.05) | (pval_panukb < 0.05),
  naive_asymmetry_select = naive_asymmetry_r12 | naive_asymmetry_r13,
  U_fixed = naive_asymmetry_r12 & !source_fdr_r12,
  U_select = (naive_asymmetry_r12 & !source_fdr_r12) |
    (naive_asymmetry_r13 & !source_fdr_r13)
)]
primary[, unsupported_inflation_cell := as.integer(U_select) - as.integer(U_fixed)]
primary[, `:=`(
  fixed_holm_min_p = pmin(1, 2 * pmin(pval_r12, pval_panukb)),
  selected_holm_min_p = pmin(1, 3 * pmin(pval_r12, pval_r13, pval_panukb))
)]
primary[, `:=`(
  fixed_holm_bh_q = p.adjust(fixed_holm_min_p, method = "BH"),
  selected_holm_bh_q = p.adjust(selected_holm_min_p, method = "BH")
)]

release_anchor <- release_rho[abs(rho_release - 0.50) < 1e-12]
setorder(release_anchor, outcome_family, trait_id)
primary[, `:=`(
  release_p_rho50 = release_anchor$p,
  release_q_rho50 = release_anchor$q_bh,
  release_q_by_rho50 = release_anchor$q_by,
  robust_large_release_drift = release_anchor$q_bh < 0.05 & release_anchor$operationally_large
)]

axis_grid <- source_rho[, .(
  fdr_difference = q_bh < 0.05,
  equivalent
), by = .(trait_id, outcome_family, release, rho_source)]
axis_wide <- dcast(axis_grid, trait_id + outcome_family + rho_source ~ release,
                   value.var = c("fdr_difference", "equivalent"))
axis_wide[, axis_agreement := fdr_difference_r12 == fdr_difference_r13 &
            equivalent_r12 == equivalent_r13]
axis_agreement_grid <- axis_wide[, .(
  agreement_n = sum(axis_agreement),
  agreement_rate = mean(axis_agreement)
), by = rho_source]

safe_spearman <- function(x, y, w = NULL) {
  ok <- is.finite(x) & is.finite(y)
  if (is.null(w)) return(cor(x[ok], y[ok], method = "spearman"))
  idx <- rep(which(ok), w[ok])
  if (length(idx) < 3L) return(NA_real_)
  suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
}

weighted_median <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  if (!length(x)) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  x[which(cumsum(w) >= sum(w) / 2)[1L]]
}

metric_observed <- data.table(
  metric = c(
    "source_delta_spearman",
    "axis_agreement_rho0",
    "minimum_axis_agreement_across_rho_source",
    "median_release_to_source_delta_ratio",
    "original_35_retained_in_r13",
    "original_35_retention_rate",
    "strict_axis_migration_rate",
    "robust_large_release_drift_rate_rho50",
    "unsupported_wording_inflation",
    "unsupported_wording_risk_ratio"
  ),
  estimate = c(
    safe_spearman(primary$delta_source_r12, primary$delta_source_r13),
    mean(primary$axis_agreement),
    min(axis_agreement_grid$agreement_rate),
    median(abs(primary$delta_release)) /
      median((abs(primary$delta_source_r12) + abs(primary$delta_source_r13)) / 2),
    sum(primary$original_common_heterogeneity & primary$source_fdr_r13),
    mean(primary$source_fdr_r13[primary$original_common_heterogeneity]),
    mean(primary$strict_axis_migration),
    mean(primary$robust_large_release_drift),
    mean(primary$unsupported_inflation_cell),
    mean(primary$U_select) / mean(primary$U_fixed)
  )
)

cluster_support <- data.table(
  pattern_component = c("B_migration", "B_robust_large_drift", "C_new_unsupported"),
  endpoints = c(
    uniqueN(primary[strict_axis_migration == TRUE]$outcome_family),
    uniqueN(primary[robust_large_release_drift == TRUE]$outcome_family),
    uniqueN(primary[unsupported_inflation_cell > 0]$outcome_family)
  ),
  exposures = c(
    uniqueN(primary[strict_axis_migration == TRUE]$trait_id),
    uniqueN(primary[robust_large_release_drift == TRUE]$trait_id),
    uniqueN(primary[unsupported_inflation_cell > 0]$trait_id)
  ),
  clinical_domains = c(
    uniqueN(primary[strict_axis_migration == TRUE]$field_domain),
    uniqueN(primary[robust_large_release_drift == TRUE]$field_domain),
    uniqueN(primary[unsupported_inflation_cell > 0]$field_domain)
  )
)

selection_metrics <- data.table(
  metric = c(
    "any_nominal_positive_fixed_r12_panukb",
    "any_nominal_positive_after_release_selection",
    "nominal_positive_inflation",
    "naive_asymmetry_fixed_r12_panukb",
    "naive_asymmetry_after_release_selection",
    "naive_asymmetry_inflation",
    "unsupported_wording_fixed",
    "unsupported_wording_after_release_selection",
    "unsupported_wording_inflation",
    "holm_within_cell_bh_positive_fixed",
    "holm_within_cell_bh_positive_selected"
  ),
  count = c(
    sum(primary$positive_fixed_any),
    sum(primary$positive_select_any),
    sum(primary$positive_select_any) - sum(primary$positive_fixed_any),
    sum(primary$naive_asymmetry_r12),
    sum(primary$naive_asymmetry_select),
    sum(primary$naive_asymmetry_select) - sum(primary$naive_asymmetry_r12),
    sum(primary$U_fixed),
    sum(primary$U_select),
    sum(primary$U_select) - sum(primary$U_fixed),
    sum(primary$fixed_holm_bh_q < 0.05),
    sum(primary$selected_holm_bh_q < 0.05)
  )
)
selection_metrics[, rate := count / nrow(primary)]

set.seed(20260726)
traits <- sort(unique(primary$trait_id))
ends <- sort(unique(primary$outcome_family))
primary[, trait_index := match(trait_id, traits)]
primary[, endpoint_index := match(outcome_family, ends)]
B <- 20000L
boot <- matrix(NA_real_, nrow = B, ncol = 8L)
colnames(boot) <- c(
  "source_delta_spearman", "axis_agreement_rho0",
  "median_release_to_source_delta_ratio", "original_35_retention_rate",
  "strict_axis_migration_rate", "robust_large_release_drift_rate_rho50",
  "unsupported_wording_inflation", "unsupported_wording_risk_ratio"
)

for (b in seq_len(B)) {
  tw <- tabulate(sample.int(length(traits), length(traits), replace = TRUE), nbins = length(traits))
  ew <- tabulate(sample.int(length(ends), length(ends), replace = TRUE), nbins = length(ends))
  w <- tw[primary$trait_index] * ew[primary$endpoint_index]
  boot[b, "source_delta_spearman"] <- safe_spearman(primary$delta_source_r12, primary$delta_source_r13, w)
  boot[b, "axis_agreement_rho0"] <- weighted.mean(primary$axis_agreement, w)
  boot[b, "median_release_to_source_delta_ratio"] <-
    weighted_median(abs(primary$delta_release), w) /
    weighted_median((abs(primary$delta_source_r12) + abs(primary$delta_source_r13)) / 2, w)
  h <- primary$original_common_heterogeneity
  boot[b, "original_35_retention_rate"] <- weighted.mean(primary$source_fdr_r13[h], w[h])
  boot[b, "strict_axis_migration_rate"] <- weighted.mean(primary$strict_axis_migration, w)
  boot[b, "robust_large_release_drift_rate_rho50"] <- weighted.mean(primary$robust_large_release_drift, w)
  boot[b, "unsupported_wording_inflation"] <- weighted.mean(primary$unsupported_inflation_cell, w)
  u_fixed <- weighted.mean(primary$U_fixed, w)
  u_select <- weighted.mean(primary$U_select, w)
  boot[b, "unsupported_wording_risk_ratio"] <- if (u_fixed > 0) u_select / u_fixed else NA_real_
  if (b %% 2000L == 0L) message("Bootstrap ", b, "/", B)
}

alpha_tail <- (1 - 0.9833) / 2
boot_ci <- rbindlist(lapply(colnames(boot), function(name) {
  x <- boot[, name]
  data.table(
    metric = name,
    bootstrap_replicates = sum(is.finite(x)),
    ci_level = 0.9833,
    ci_low = unname(quantile(x, alpha_tail, na.rm = TRUE, type = 7)),
    ci_high = unname(quantile(x, 1 - alpha_tail, na.rm = TRUE, type = 7))
  )
}))
metric_observed <- merge(metric_observed, boot_ci, by = "metric", all.x = TRUE)

get_est <- function(name) metric_observed[metric == name, estimate]
get_low <- function(name) metric_observed[metric == name, ci_low]
get_high <- function(name) metric_observed[metric == name, ci_high]
support <- function(name, col) cluster_support[pattern_component == name, get(col)]

A_conditions <- data.table(
  condition = c("spearman_point", "spearman_lower", "axis_agreement_point",
                "axis_agreement_lower", "axis_agreement_rho_grid", "drift_source_ratio_point",
                "drift_source_ratio_upper", "retained_count", "retention_point", "retention_lower"),
  pass = c(
    get_est("source_delta_spearman") >= 0.90,
    get_low("source_delta_spearman") >= 0.85,
    get_est("axis_agreement_rho0") >= 0.90,
    get_low("axis_agreement_rho0") >= 0.85,
    get_est("minimum_axis_agreement_across_rho_source") >= 0.80,
    get_est("median_release_to_source_delta_ratio") <= 0.50,
    get_high("median_release_to_source_delta_ratio") <= 0.75,
    get_est("original_35_retained_in_r13") >= 28,
    get_est("original_35_retention_rate") >= 0.80,
    get_low("original_35_retention_rate") >= 0.60
  )
)
B_conditions <- data.table(
  condition = c("migration_point", "migration_lower", "migration_endpoints",
                "migration_exposures", "migration_domains", "large_drift_point", "large_drift_lower"),
  pass = c(
    get_est("strict_axis_migration_rate") >= 0.10,
    get_low("strict_axis_migration_rate") >= 0.05,
    support("B_migration", "endpoints") >= 10,
    support("B_migration", "exposures") >= 5,
    support("B_migration", "clinical_domains") >= 3,
    get_est("robust_large_release_drift_rate_rho50") >= 0.05,
    get_low("robust_large_release_drift_rate_rho50") >= 0.02
  )
)
C_conditions <- data.table(
  condition = c("inflation_point", "inflation_lower", "risk_ratio",
                "new_endpoints", "new_exposures", "new_domains"),
  pass = c(
    get_est("unsupported_wording_inflation") >= 0.05,
    get_low("unsupported_wording_inflation") >= 0.025,
    get_est("unsupported_wording_risk_ratio") >= 1.20,
    support("C_new_unsupported", "endpoints") >= 10,
    support("C_new_unsupported", "exposures") >= 5,
    support("C_new_unsupported", "clinical_domains") >= 3
  )
)
pattern_gate <- rbindlist(list(
  cbind(pattern = "A", A_conditions),
  cbind(pattern = "B", B_conditions),
  cbind(pattern = "C", C_conditions)
))
pattern_summary <- pattern_gate[, .(conditions_passed = sum(pass), conditions_total = .N,
                                    pattern_pass = all(pass)), by = pattern]
integration_decision <- if (any(pattern_summary$pattern_pass)) "GO_MAIN_TEXT" else "REFRAME_SUPPLEMENT"

fwrite(primary, tagged_path(file.path(root, "outputs/RELEASE_STABILITY_CELL_RESULTS.tsv")), sep = "\t", na = "NA")
fwrite(source_rho, tagged_path(file.path(root, "outputs/RHO_SOURCE_SENSITIVITY.tsv")), sep = "\t", na = "NA")
fwrite(release_rho, tagged_path(file.path(root, "outputs/RHO_RELEASE_SENSITIVITY.tsv")), sep = "\t", na = "NA")
fwrite(axis_agreement_grid, tagged_path(file.path(root, "outputs/AXIS_AGREEMENT_RHO_GRID.tsv")), sep = "\t", na = "NA")
fwrite(metric_observed, tagged_path(file.path(root, "outputs/PATTERN_ABC_METRICS.tsv")), sep = "\t", na = "NA")
fwrite(cluster_support, tagged_path(file.path(root, "outputs/PATTERN_CLUSTER_SUPPORT.tsv")), sep = "\t", na = "NA")
fwrite(selection_metrics, tagged_path(file.path(root, "outputs/SOURCE_RELEASE_SELECTION_METRICS.tsv")), sep = "\t", na = "NA")
fwrite(boot_ci, tagged_path(file.path(root, "outputs/CLUSTER_BOOTSTRAP_CI.tsv")), sep = "\t", na = "NA")
fwrite(pattern_gate, tagged_path(file.path(root, "outputs/PATTERN_ABC_GATE_DETAILS.tsv")), sep = "\t", na = "NA")
fwrite(pattern_summary, tagged_path(file.path(root, "outputs/PATTERN_ABC_GATE_SUMMARY.tsv")), sep = "\t", na = "NA")
saveRDS(boot, tagged_path(file.path(root, "outputs/CLUSTER_BOOTSTRAP_20000.rds")), compress = "xz")
writeLines(c(
  paste0("decision\t", integration_decision),
  paste0("passed_patterns\t", paste(pattern_summary[pattern_pass == TRUE, pattern], collapse = ",")),
  "review_requirement\tIndependent scientific review recommended before interpretive use"
), tagged_path(file.path(root, "GO_NO_GO_RESULT_STAGE.tsv")))
writeLines(capture.output(sessionInfo()), tagged_path(file.path(root, "qa/sessionInfo_release_stability_analysis.txt")))
message("Analysis complete: ", integration_decision)
