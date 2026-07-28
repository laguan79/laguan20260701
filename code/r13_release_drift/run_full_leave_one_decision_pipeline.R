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
d <- fread(estimate_path)
stopifnot(nrow(d) == 1060L)
margin_source <- log(1.20)
margin_release <- log(1.10)

contrast <- function(b1, se1, b2, se2, rho) {
  delta <- b1 - b2
  se <- sqrt(se1^2 + se2^2 - 2 * rho * se1 * se2)
  p <- 2 * pnorm(-abs(delta / se))
  data.table(delta = delta, se = se, p = p,
             low90 = delta - qnorm(0.95) * se,
             high90 = delta + qnorm(0.95) * se)
}

compute_pipeline <- function(x) {
  s12 <- contrast(x$b_r12, x$se_r12, x$b_panukb, x$se_panukb, 0)
  s13 <- contrast(x$b_r13, x$se_r13, x$b_panukb, x$se_panukb, 0)
  q <- p.adjust(c(s12$p, s13$p), method = "BH")
  f12 <- q[seq_len(nrow(x))] < 0.05
  f13 <- q[nrow(x) + seq_len(nrow(x))] < 0.05
  e12 <- s12$low90 >= -margin_source & s12$high90 <= margin_source
  e13 <- s13$low90 >= -margin_source & s13$high90 <= margin_source
  migration <- (f12 != f13) | (e12 != e13)
  release <- contrast(x$b_r13, x$se_r13, x$b_r12, x$se_r12, 0.50)
  release_q <- p.adjust(release$p, method = "BH")
  robust <- release_q < 0.05 & abs(release$delta) >= margin_release
  a12 <- xor(x$pval_r12 < 0.05, x$pval_panukb < 0.05)
  a13 <- xor(x$pval_r13 < 0.05, x$pval_panukb < 0.05)
  uf <- a12 & !f12
  us <- (a12 & !f12) | (a13 & !f13)
  h <- x$common_heterogeneity == TRUE
  data.table(
    cells = nrow(x),
    source_delta_spearman = cor(s12$delta, s13$delta, method = "spearman"),
    axis_agreement_rho0 = mean(!migration),
    median_release_to_source_delta_ratio = median(abs(release$delta)) /
      median((abs(s12$delta) + abs(s13$delta)) / 2),
    original_calls_remaining = sum(h),
    original_call_retained_count = sum(h & f13),
    original_call_retention_rate = mean(f13[h]),
    r12_joint_fdr_count = sum(f12),
    r13_joint_fdr_count = sum(f13),
    strict_axis_migration_rate = mean(migration),
    robust_large_release_drift_rate_rho50 = mean(robust),
    unsupported_wording_inflation = mean(as.integer(us) - as.integer(uf)),
    unsupported_wording_risk_ratio = mean(us) / mean(uf)
  )
}

rows <- list()
pos <- 0L
for (type in c("none", "endpoint", "exposure")) {
  ids <- if (type == "none") "none" else if (type == "endpoint") unique(d$outcome_family) else unique(d$trait_id)
  for (id in ids) {
    pos <- pos + 1L
    x <- if (type == "none") d else if (type == "endpoint") d[outcome_family != id] else d[trait_id != id]
    rows[[pos]] <- cbind(data.table(leave_one_type = type, leave_one_id = id), compute_pipeline(x))
  }
}
out <- rbindlist(rows)
fwrite(
  out,
  tagged_path(file.path(root, "outputs/FULL_LEAVE_ONE_DECISION_PIPELINE.tsv")),
  sep = "\t",
  na = "NA"
)

long <- melt(out[leave_one_type != "none"],
             id.vars = c("leave_one_type", "leave_one_id"),
             variable.name = "metric", value.name = "estimate")
summary <- long[, .(
  minimum = min(estimate, na.rm = TRUE),
  median = median(estimate, na.rm = TRUE),
  maximum = max(estimate, na.rm = TRUE)
), by = .(leave_one_type, metric)]
fwrite(
  summary,
  tagged_path(file.path(root, "outputs/FULL_LEAVE_ONE_DECISION_PIPELINE_SUMMARY.tsv")),
  sep = "\t",
  na = "NA"
)

gate <- data.table(
  check = c("full_dataset_plus_73_omissions", "endpoint_omissions", "exposure_omissions",
            "multiplicity_recomputed", "equivalence_recomputed", "release_fdr_recomputed"),
  observed = c(as.character(nrow(out)), as.character(sum(out$leave_one_type == "endpoint")),
               as.character(sum(out$leave_one_type == "exposure")),
               "joint BH within each omission", "90% CI within each omission",
               "BH at rho_release=0.50 within each omission"),
  required = c("74", "53", "20", "yes", "yes", "yes"),
  pass = c(nrow(out) == 74L, sum(out$leave_one_type == "endpoint") == 53L,
           sum(out$leave_one_type == "exposure") == 20L, TRUE, TRUE, TRUE)
)
fwrite(
  gate,
  tagged_path(file.path(root, "qa/FULL_LEAVE_ONE_DECISION_PIPELINE_GATE.tsv")),
  sep = "\t",
  na = "NA"
)
writeLines(
  capture.output(sessionInfo()),
  tagged_path(file.path(root, "qa/sessionInfo_full_leave_one.txt"))
)
if (!all(gate$pass)) stop("Full leave-one decision pipeline gate failed")
message("Full leave-one decision pipeline recomputed for 73 cluster omissions")
