#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
locked_path <- file.path(root, "outputs/R13_LOCKED_VARIANT_EXTRACTION.tsv")
exact_path <- file.path(root, "outputs/R13_TARGETED_EXACT_TABIX.tsv")
exact_gate_path <- file.path(root, "qa/R13_TARGETED_EXACT_TABIX_GATE.tsv")
rp_path <- file.path(root, "inputs/LOCKED_R12_PAN_COMMON_SNP_EFFECTS.tsv")
primary_path <- file.path(root, "outputs/THREE_SOURCE_RELEASE_ESTIMATES.tsv")
stopifnot(
  file.exists(locked_path), file.exists(exact_path), file.exists(exact_gate_path),
  file.exists(rp_path), file.exists(primary_path)
)
if (!all(fread(exact_gate_path)$pass)) stop("Targeted exact-Tabix gate is not fully passed")

locked <- fread(locked_path)
exact <- fread(exact_path)[tabix_status == "ok"]
rp <- fread(rp_path)
primary <- fread(primary_path)
stopifnot(nrow(locked) == 42200L, nrow(primary) == 1060L)

key_cols <- c(
  "finngen_phenocode", "query_chr", "query_pos", "query_ref", "query_alt"
)
if (exact[, uniqueN(query_id)] != nrow(exact)) stop("Exact query rows are not unique")
exact_key <- unique(
  exact[, c(key_cols, "se_tabix"), with = FALSE],
  by = key_cols
)

locked[, `:=`(
  se_r13_primary = se_r13,
  exact_se_replaced = FALSE
)]
locked[
  exact_key,
  on = c(
    "finngen_phenocode", "query_chr", "query_pos",
    "query_ref", "query_alt"
  ),
  `:=`(
    se_r13 = i.se_tabix,
    exact_se_replaced = TRUE
  )
]
if (any(!is.finite(locked$se_r13) | locked$se_r13 <= 0)) {
  stop("Non-finite or non-positive R13 SE after exact replacement")
}

r13 <- locked[, .(
  trait_id, outcome_family, finngen_phenocode, SNP,
  beta_outcome_r13 = beta_r13_aligned,
  se_outcome_r13 = se_r13
)]
x <- merge(
  rp,
  r13,
  by = c("trait_id", "outcome_family", "finngen_phenocode", "SNP"),
  all = FALSE
)
setorder(x, outcome_family, trait_id, SNP)
stopifnot(nrow(x) == 42200L)

run_ivw <- function(bx, by, sx, sy) {
  z <- mr_ivw(bx, by, sx, sy, parameters = default_parameters())
  data.table(
    nsnp = z$nsnp, b = z$b, se = z$se, pval = z$pval,
    Q = z$Q, Q_df = z$Q_df, Q_pval = z$Q_pval
  )
}

r13_est <- x[, run_ivw(
  beta_exposure, beta_outcome_r13, se_exposure, se_outcome_r13
), by = .(trait_id, outcome_family, finngen_phenocode)]
stopifnot(nrow(r13_est) == 1060L)

sensitivity <- copy(primary)
sensitivity[
  r13_est,
  on = .(trait_id, outcome_family, finngen_phenocode),
  `:=`(
    nsnp_r13 = i.nsnp,
    b_r13 = i.b,
    se_r13 = i.se,
    pval_r13 = i.pval,
    Q_r13 = i.Q,
    Q_df_r13 = i.Q_df,
    Q_pval_r13 = i.Q_pval
  )
]
sensitivity[, `:=`(
  delta_release = b_r13 - b_r12,
  delta_source_r12 = b_r12 - b_panukb,
  delta_source_r13 = b_r13 - b_panukb,
  algebra_identity_error =
    (b_r13 - b_panukb) - (b_r12 - b_panukb) - (b_r13 - b_r12)
)]
setorder(sensitivity, outcome_family, trait_id)

comparison <- merge(
  primary[, .(
    trait_id, outcome_family,
    b_r13_primary = b_r13,
    se_r13_primary_cell = se_r13,
    pval_r13_primary = pval_r13
  )],
  sensitivity[, .(
    trait_id, outcome_family,
    b_r13_exact_sensitivity = b_r13,
    se_r13_exact_sensitivity = se_r13,
    pval_r13_exact_sensitivity = pval_r13
  )],
  by = c("trait_id", "outcome_family")
)
comparison[, `:=`(
  b_abs_change = abs(b_r13_exact_sensitivity - b_r13_primary),
  se_abs_change = abs(se_r13_exact_sensitivity - se_r13_primary_cell),
  p_abs_change = abs(pval_r13_exact_sensitivity - pval_r13_primary)
)]

fwrite(
  locked,
  file.path(
    root,
    "outputs/R13_LOCKED_VARIANT_EXTRACTION_TARGETED_EXACT_SENSITIVITY.tsv"
  ),
  sep = "\t",
  na = "NA"
)
fwrite(
  sensitivity,
  file.path(
    root,
    "outputs/THREE_SOURCE_RELEASE_ESTIMATES_TARGETED_EXACT.tsv"
  ),
  sep = "\t",
  na = "NA"
)
fwrite(
  comparison,
  file.path(root, "qa/TARGETED_EXACT_ESTIMATE_COMPARISON.tsv"),
  sep = "\t",
  na = "NA"
)

gate <- data.table(
  check = c(
    "locked_variant_rows",
    "target_exact_query_rows",
    "expanded_rows_replaced",
    "estimate_cells",
    "finite_exact_sensitivity",
    "algebra_identity"
  ),
  observed = c(
    as.character(nrow(locked)),
    as.character(nrow(exact)),
    as.character(sum(locked$exact_se_replaced)),
    as.character(nrow(sensitivity)),
    sprintf(
      "%d/%d",
      sum(
        is.finite(sensitivity$b_r13) &
          is.finite(sensitivity$se_r13) &
          sensitivity$se_r13 > 0
      ),
      nrow(sensitivity)
    ),
    sprintf(
      "max %.3g",
      max(abs(sensitivity$algebra_identity_error), na.rm = TRUE)
    )
  ),
  required = c(
    "42200", "3517", ">0", "1060", "all", "<=1e-12"
  ),
  pass = c(
    nrow(locked) == 42200L,
    nrow(exact) == 3517L,
    sum(locked$exact_se_replaced) > 0,
    nrow(sensitivity) == 1060L,
    all(
      is.finite(sensitivity$b_r13) &
        is.finite(sensitivity$se_r13) &
        sensitivity$se_r13 > 0
    ),
    max(abs(sensitivity$algebra_identity_error), na.rm = TRUE) <= 1e-12
  )
)
fwrite(
  gate,
  file.path(root, "qa/TARGETED_EXACT_ESTIMATION_GATE.tsv"),
  sep = "\t",
  na = "NA"
)
writeLines(
  capture.output(sessionInfo()),
  file.path(root, "qa/sessionInfo_targeted_exact_estimation.txt")
)
if (!all(gate$pass)) stop("Targeted exact-SE estimation gate failed")
message(
  "Targeted exact-SE estimates written; expanded rows replaced=",
  sum(locked$exact_se_replaced)
)
