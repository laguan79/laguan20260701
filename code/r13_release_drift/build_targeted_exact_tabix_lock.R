#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
cells <- fread(file.path(root, "outputs/RELEASE_STABILITY_CELL_RESULTS.tsv"))
locked <- fread(file.path(root, "inputs/LOCKED_COMMON_SNP_IDENTITIES.tsv"))
api <- fread(file.path(root, "outputs/R13_API_UNIQUE_VARIANT_EFFECTS.tsv"))
queries <- fread(file.path(root, "outputs/R13_API_QUERY_INDEX.tsv"))

cells[, `:=`(
  audit_historical35 = original_common_heterogeneity == TRUE,
  audit_axis_migration = strict_axis_migration == TRUE,
  audit_robust_drift = robust_large_release_drift == TRUE,
  audit_new_unsupported = unsupported_inflation_cell > 0
)]
cell_targets <- cells[
  audit_historical35 | audit_axis_migration | audit_robust_drift | audit_new_unsupported,
  .(trait_id, outcome_family, audit_historical35, audit_axis_migration,
    audit_robust_drift, audit_new_unsupported)
]

row_targets <- merge(
  locked,
  cell_targets,
  by = c("trait_id", "outcome_family"),
  all = FALSE
)[, .(
  audit_historical35 = any(audit_historical35),
  audit_axis_migration = any(audit_axis_migration),
  audit_robust_drift = any(audit_robust_drift),
  audit_new_unsupported = any(audit_new_unsupported)
), by = .(
  finngen_phenocode,
  variant_key = r12_selected_variant_key,
  query_chr,
  query_pos,
  query_ref,
  query_alt
)]

api_meta <- merge(
  queries[, .(query_id, finngen_phenocode, variant_key, query_chr, query_pos, query_ref, query_alt)],
  api[, .(query_id, api_status, beta_alt_api = beta_alt, pval_api = pval,
          se_api = sebeta_reconstructed)],
  by = "query_id",
  all.x = TRUE
)
api_meta[, audit_extreme_near_null := is.finite(pval_api) & pval_api >= 0.999]
api_meta[, endpoint_p_rank := frank(pval_api, ties.method = "first"), by = finngen_phenocode]
api_meta[, endpoint_p_reverse_rank := frank(-pval_api, ties.method = "first"), by = finngen_phenocode]
api_meta[, audit_endpoint_p_extreme := endpoint_p_rank == 1L | endpoint_p_reverse_rank == 1L]

target <- merge(
  api_meta,
  row_targets,
  by = c("finngen_phenocode", "variant_key", "query_chr", "query_pos", "query_ref", "query_alt"),
  all.x = TRUE
)
for (nm in c("audit_historical35", "audit_axis_migration", "audit_robust_drift", "audit_new_unsupported")) {
  set(target, which(is.na(target[[nm]])), nm, FALSE)
}
target <- target[
  audit_historical35 | audit_axis_migration | audit_robust_drift |
    audit_new_unsupported | audit_extreme_near_null | audit_endpoint_p_extreme
]
target[, audit_reason := paste(
  c("historical35", "axis_migration", "robust_drift", "new_unsupported",
    "p_ge_0.999", "endpoint_p_extreme")[
      c(audit_historical35[1], audit_axis_migration[1], audit_robust_drift[1],
        audit_new_unsupported[1], audit_extreme_near_null[1], audit_endpoint_p_extreme[1])
    ],
  collapse = ";"
), by = query_id]
setorder(target, finngen_phenocode, query_chr, query_pos, query_ref, query_alt)
stopifnot(uniqueN(target$query_id) == nrow(target), all(target$api_status == "ok"))

path <- file.path(root, "inputs/R13_TARGETED_EXACT_TABIX_QUERY_LOCK.tsv")
fwrite(target, path, sep = "\t", na = "NA")
hash <- digest(path, algo = "sha256", file = TRUE)
writeLines(paste(hash, basename(path)), file.path(root, "inputs/R13_TARGETED_EXACT_TABIX_QUERY_LOCK.sha256"))

summary <- data.table(
  criterion = c("all_targets", "historical35", "axis_migration", "robust_drift",
                "new_unsupported", "p_ge_0.999", "endpoint_p_extreme"),
  n_queries = c(
    nrow(target), sum(target$audit_historical35), sum(target$audit_axis_migration),
    sum(target$audit_robust_drift), sum(target$audit_new_unsupported),
    sum(target$audit_extreme_near_null), sum(target$audit_endpoint_p_extreme)
  )
)
fwrite(summary, file.path(root, "qa/R13_TARGETED_EXACT_TABIX_LOCK_SUMMARY.tsv"), sep = "\t", na = "NA")
message("Locked ", nrow(target), " exact-Tabix forensic queries; sha256=", hash)

