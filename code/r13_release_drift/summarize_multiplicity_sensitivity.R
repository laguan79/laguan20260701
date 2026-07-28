#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
cell <- fread(file.path(root, "outputs/RELEASE_STABILITY_CELL_RESULTS.tsv"))
source_rho <- fread(file.path(root, "outputs/RHO_SOURCE_SENSITIVITY.tsv"))
release_rho <- fread(file.path(root, "outputs/RHO_RELEASE_SENSITIVITY.tsv"))
stopifnot(nrow(cell) == 1060L)

p12 <- cell$source_p_r12
p13 <- cell$source_p_r13
h <- cell$original_common_heterogeneity == TRUE
q12_sep_bh <- p.adjust(p12, "BH")
q13_sep_bh <- p.adjust(p13, "BH")
q12_sep_by <- p.adjust(p12, "BY")
q13_sep_by <- p.adjust(p13, "BY")
q_joint_by <- p.adjust(c(p12, p13), "BY")
f12_joint_by <- q_joint_by[seq_len(nrow(cell))] < 0.05
f13_joint_by <- q_joint_by[nrow(cell) + seq_len(nrow(cell))] < 0.05

calc_secondary <- function(f12, f13, label) {
  migration <- (f12 != f13) | (cell$source_equivalent_r12 != cell$source_equivalent_r13)
  a12 <- cell$naive_asymmetry_r12
  a13 <- cell$naive_asymmetry_r13
  uf <- a12 & !f12
  us <- (a12 & !f12) | (a13 & !f13)
  data.table(
    family = label,
    r12_calls = sum(f12),
    r13_calls = sum(f13),
    shared_calls = sum(f12 & f13),
    historical35_retained_r12 = sum(h & f12),
    historical35_retained_r13 = sum(h & f13),
    strict_axis_migration_rate = mean(migration),
    unsupported_wording_inflation = mean(as.integer(us) - as.integer(uf)),
    unsupported_wording_risk_ratio = mean(us) / mean(uf)
  )
}

summary <- rbindlist(list(
  calc_secondary(cell$source_fdr_r12, cell$source_fdr_r13, "joint_2N_BH_primary"),
  calc_secondary(q12_sep_bh < 0.05, q13_sep_bh < 0.05, "separate_N_BH_sensitivity"),
  calc_secondary(f12_joint_by, f13_joint_by, "joint_2N_BY_sensitivity"),
  calc_secondary(q12_sep_by < 0.05, q13_sep_by < 0.05, "separate_N_BY_sensitivity")
))
fwrite(summary, file.path(root, "outputs/MULTIPLICITY_FAMILY_SENSITIVITY.tsv"), sep = "\t", na = "NA")

retention_rho <- source_rho[, .(
  fdr_bh = q_bh < 0.05,
  fdr_by = q_by < 0.05
), by = .(trait_id, outcome_family, release, rho_source)]
retention_rho <- merge(retention_rho,
                       cell[, .(trait_id, outcome_family, historical35 = original_common_heterogeneity)],
                       by = c("trait_id", "outcome_family"), all.x = TRUE)
retention_rho <- retention_rho[release == "r13", .(
  r13_bh_calls = sum(fdr_bh),
  r13_by_calls = sum(fdr_by),
  historical35_retained_bh = sum(historical35 & fdr_bh),
  historical35_retained_by = sum(historical35 & fdr_by)
), by = rho_source]
fwrite(retention_rho, file.path(root, "outputs/HISTORICAL35_RETENTION_RHO_SENSITIVITY.tsv"), sep = "\t", na = "NA")

release_summary <- release_rho[, .(
  bh_calls = sum(q_bh < 0.05),
  by_calls = sum(q_by < 0.05),
  bh_operationally_large = sum(q_bh < 0.05 & operationally_large),
  by_operationally_large = sum(q_by < 0.05 & operationally_large)
), by = rho_release]
fwrite(release_summary, file.path(root, "outputs/RELEASE_RHO_MULTIPLICITY_SENSITIVITY.tsv"), sep = "\t", na = "NA")
writeLines(capture.output(sessionInfo()), file.path(root, "qa/sessionInfo_multiplicity_sensitivity.txt"))
message("Multiplicity-family and rho sensitivity summaries written")
