#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
snp_path <- file.path(root, "outputs/THREE_SOURCE_SNP_EFFECTS.tsv")
cell_path <- file.path(root, "outputs/RELEASE_STABILITY_CELL_RESULTS.tsv")
stopifnot(file.exists(snp_path), file.exists(cell_path))
s <- fread(snp_path)
c <- fread(cell_path)
stopifnot(nrow(s) == 42200L, nrow(c) == 1060L)

endpoint <- s[, .(
  n_rows = .N,
  pearson_beta_r12_r13 = cor(beta_outcome_r12, beta_outcome_r13),
  spearman_beta_r12_r13 = cor(beta_outcome_r12, beta_outcome_r13, method = "spearman"),
  median_se_ratio_r13_r12 = median(se_outcome_r13 / se_outcome_r12),
  median_absolute_beta_drift = median(abs(beta_outcome_r13 - beta_outcome_r12)),
  beta_sign_agreement = mean(sign(beta_outcome_r13) == sign(beta_outcome_r12))
), by = .(outcome_family, finngen_phenocode)]
setorder(endpoint, pearson_beta_r12_r13)
fwrite(endpoint, file.path(root, "qa/SNP_LEVEL_RELEASE_ENDPOINT_AUDIT.tsv"), sep = "\t", na = "NA")

overall <- data.table(
  metric = c(
    "snp_rows", "pearson_beta_r12_r13", "spearman_beta_r12_r13",
    "median_se_ratio_r13_r12", "median_absolute_beta_drift", "beta_sign_agreement",
    "minimum_endpoint_pearson", "minimum_endpoint_spearman"
  ),
  value = c(
    nrow(s), cor(s$beta_outcome_r12, s$beta_outcome_r13),
    cor(s$beta_outcome_r12, s$beta_outcome_r13, method = "spearman"),
    median(s$se_outcome_r13 / s$se_outcome_r12),
    median(abs(s$beta_outcome_r13 - s$beta_outcome_r12)),
    mean(sign(s$beta_outcome_r13) == sign(s$beta_outcome_r12)),
    min(endpoint$pearson_beta_r12_r13), min(endpoint$spearman_beta_r12_r13)
  )
)
fwrite(overall, file.path(root, "qa/SNP_LEVEL_RELEASE_OVERALL_AUDIT.tsv"), sep = "\t", na = "NA")

original35 <- c[original_common_heterogeneity == TRUE, .(
  trait_id, outcome_family, field_domain,
  delta_source_r12, delta_source_r13, delta_release,
  source_q_r12, source_q_r13,
  retained_r12_joint_fdr = source_fdr_r12,
  retained_r13_joint_fdr = source_fdr_r13,
  class_r12, class_r13,
  strict_axis_migration, full_class_migration
)]
setorder(original35, source_q_r13)
stopifnot(nrow(original35) == 35L)
fwrite(original35, file.path(root, "qa/ORIGINAL_35_RELEASE_AUDIT.tsv"), sep = "\t", na = "NA")

gate <- data.table(
  check = c("snp_row_contract", "finite_snp_effects", "endpoint_coverage",
            "api_tabix_precision_gate", "interpretation_boundary"),
  observed = c(
    as.character(nrow(s)),
    sprintf("%d invalid", sum(!is.finite(s$beta_outcome_r12) | !is.finite(s$beta_outcome_r13) |
                                !is.finite(s$se_outcome_r12) | !is.finite(s$se_outcome_r13) |
                                s$se_outcome_r12 <= 0 | s$se_outcome_r13 <= 0)),
    as.character(nrow(endpoint)),
    if (all(fread(file.path(root, "qa/R13_API_TABIX_PRECISION_GATE.tsv"))$pass)) "PASS" else "FAIL",
    "descriptive QA; no independent-replication claim"
  ),
  required = c("42200", "0", "53", "PASS", "release-drift boundary"),
  pass = c(
    nrow(s) == 42200L,
    !any(!is.finite(s$beta_outcome_r12) | !is.finite(s$beta_outcome_r13) |
           !is.finite(s$se_outcome_r12) | !is.finite(s$se_outcome_r13) |
           s$se_outcome_r12 <= 0 | s$se_outcome_r13 <= 0),
    nrow(endpoint) == 53L,
    all(fread(file.path(root, "qa/R13_API_TABIX_PRECISION_GATE.tsv"))$pass),
    TRUE
  )
)
fwrite(gate, file.path(root, "qa/SNP_LEVEL_RELEASE_GATE.tsv"), sep = "\t", na = "NA")
writeLines(capture.output(sessionInfo()), file.path(root, "qa/sessionInfo_snp_release_audit.txt"))
if (!all(gate$pass)) stop("SNP-level release QA failed")
message("SNP-level release QA passed")
