#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(TwoSampleMR)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
rp_path <- file.path(root, "inputs/LOCKED_R12_PAN_COMMON_SNP_EFFECTS.tsv")
rp_hash_path <- file.path(root, "inputs/LOCKED_R12_PAN_COMMON_SNP_EFFECTS.sha256")
r13_path <- file.path(root, "outputs/R13_LOCKED_VARIANT_EXTRACTION.tsv")
reference_path <- file.path(root, "outputs/THREE_SOURCE_RELEASE_ESTIMATES_LONG.tsv")
stopifnot(file.exists(rp_path), file.exists(rp_hash_path), file.exists(r13_path), file.exists(reference_path))

expected <- strsplit(readLines(rp_hash_path, warn = FALSE)[1], " ", fixed = TRUE)[[1]][1]
observed <- digest(rp_path, algo = "sha256", file = TRUE)
if (!identical(tolower(expected), tolower(observed))) stop("Packaged R12/Pan input hash mismatch")

rp <- fread(rp_path)
r13 <- fread(r13_path)[extraction_status == "retained", .(
  trait_id, outcome_family, finngen_phenocode, SNP,
  beta_outcome_r13 = beta_r13_aligned, se_outcome_r13 = se_r13
)]
x <- merge(rp, r13, by = c("trait_id", "outcome_family", "finngen_phenocode", "SNP"), all = FALSE)
setorder(x, outcome_family, trait_id, SNP)
stopifnot(nrow(x) == 42200L)

run_ivw <- function(bx, by, sx, sy) {
  z <- mr_ivw(bx, by, sx, sy, parameters = default_parameters())
  data.table(nsnp = z$nsnp, b = z$b, se = z$se, pval = z$pval,
             Q = z$Q, Q_df = z$Q_df, Q_pval = z$Q_pval)
}

rows <- x[, {
  src <- list(
    r12 = list(by = beta_outcome_r12, sy = se_outcome_r12),
    r13 = list(by = beta_outcome_r13, sy = se_outcome_r13),
    panukb = list(by = beta_outcome_panukb, sy = se_outcome_panukb)
  )
  rbindlist(lapply(names(src), function(name) {
    cbind(data.table(source = name), run_ivw(beta_exposure, src[[name]]$by, se_exposure, src[[name]]$sy))
  }))
}, by = .(trait_id, outcome_family, finngen_phenocode)]
setorder(rows, outcome_family, trait_id, source)
ref <- fread(reference_path)
setorder(ref, outcome_family, trait_id, source)
cmp <- merge(
  rows,
  ref[, .(trait_id, outcome_family, source, b_ref = b, se_ref = se, pval_ref = pval)],
  by = c("trait_id", "outcome_family", "source"),
  all = TRUE
)
cmp[, `:=`(b_abs_error = abs(b - b_ref), se_abs_error = abs(se - se_ref),
           p_abs_error = abs(pval - pval_ref))]
fwrite(rows, file.path(root, "outputs/THREE_SOURCE_RELEASE_ESTIMATES_PACKAGED_REPRO.tsv"), sep = "\t", na = "NA")
fwrite(cmp, file.path(root, "qa/PACKAGED_INPUT_REPRODUCTION_AUDIT.tsv"), sep = "\t", na = "NA")

gate <- data.table(
  check = c("derived_input_hash", "derived_rows", "estimate_rows", "beta_reproduction",
            "se_reproduction", "p_reproduction"),
  observed = c(observed, as.character(nrow(rp)), as.character(nrow(rows)),
               sprintf("max %.3g", max(cmp$b_abs_error)),
               sprintf("max %.3g", max(cmp$se_abs_error)),
               sprintf("max %.3g", max(cmp$p_abs_error))),
  required = c(expected, "42200", "3180", "<=1e-10", "<=1e-10", "<=1e-10"),
  pass = c(identical(tolower(expected), tolower(observed)), nrow(rp) == 42200L,
           nrow(rows) == 3180L, max(cmp$b_abs_error) <= 1e-10,
           max(cmp$se_abs_error) <= 1e-10, max(cmp$p_abs_error) <= 1e-10)
)
fwrite(gate, file.path(root, "qa/PACKAGED_INPUT_REPRODUCTION_GATE.tsv"), sep = "\t", na = "NA")
writeLines(capture.output(sessionInfo()), file.path(root, "qa/sessionInfo_packaged_reproduction.txt"))
if (!all(gate$pass)) stop("Package-relative reproduction gate failed")
message("Package-relative three-way reproduction gate passed")
