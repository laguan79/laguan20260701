#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (is.na(script_path) || !nzchar(script_path)) {
  script_path <- sys.frame(1)$ofile
}
script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
release_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
project_parent <- dirname(release_root)
input_root <- Sys.getenv("CALIBRATION_INPUT_ROOT", unset = "")
if (!nzchar(input_root)) stop("Set CALIBRATION_INPUT_ROOT to the restored curated calibration-input directory before running this provenance script.")
input_root <- normalizePath(input_root, winslash = "/", mustWork = TRUE)
batch_root <- Sys.getenv("CALIBRATION_BATCH_ROOT", unset = file.path(input_root, "batches"))
original_path <- Sys.getenv("CALIBRATION_CLASSIFICATION_TABLE", unset = file.path(input_root, "paired_uncertainty_classification.tsv"))
out_dir <- file.path(release_root, "common_snp_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

original <- fread(original_path)
batch_dirs <- list.dirs(batch_root, recursive = FALSE, full.names = TRUE)
mr_files <- unlist(lapply(
  batch_dirs,
  function(d) list.files(
    file.path(d, "mr_results"),
    pattern = "_mr_results[.]tsv$",
    full.names = TRUE
  )
))

read_index_row <- function(path) {
  x <- fread(path, nrows = 1L)
  stem <- sub("_mr_results[.]tsv$", "", basename(path))
  harm_path <- file.path(dirname(dirname(path)), "harmonized", paste0(stem, "_harmonised.tsv"))
  data.table(
    batch_id = basename(dirname(dirname(path))),
    trait_id = x$trait_id[1],
    trait_label = x$trait_label[1],
    outcome_family = x$outcome_family[1],
    endpoint_key = x$endpoint_key[1],
    source = x$source[1],
    independent_source_group = x$independent_source_group[1],
    outcome_label = x$outcome_label[1],
    harmonized_path = harm_path,
    harmonized_exists = file.exists(harm_path)
  )
}

message("Indexing ", length(mr_files), " source-result files")
source_index <- rbindlist(lapply(mr_files, read_index_row), fill = TRUE)
source_index <- unique(
  source_index,
  by = c("trait_id", "outcome_family", "independent_source_group")
)

if (nrow(source_index) != 2120L || any(!source_index$harmonized_exists)) {
  stop(
    "Source index contract failed: rows=", nrow(source_index),
    "; missing harmonized files=", sum(!source_index$harmonized_exists)
  )
}

prepare_harmonized <- function(path) {
  d <- fread(path)
  d <- d[
    as.logical(mr_keep) &
      is.finite(beta.exposure) &
      is.finite(beta.outcome) &
      is.finite(se.exposure) &
      is.finite(se.outcome) &
      se.outcome > 0
  ]
  unique(d, by = "SNP")
}

run_ivw <- function(d) {
  if (nrow(d) < 3L) {
    return(data.table(
      nsnp_common = nrow(d),
      b_common = NA_real_,
      se_common = NA_real_,
      pval_common = NA_real_,
      Q_common = NA_real_,
      Q_df_common = NA_real_,
      Q_pval_common = NA_real_
    ))
  }
  z <- mr_ivw(
    d$beta.exposure,
    d$beta.outcome,
    d$se.exposure,
    d$se.outcome,
    parameters = default_parameters()
  )
  data.table(
    nsnp_common = z$nsnp,
    b_common = z$b,
    se_common = z$se,
    pval_common = z$pval,
    Q_common = z$Q,
    Q_df_common = z$Q_df,
    Q_pval_common = z$Q_pval
  )
}

pair_keys <- unique(source_index[, .(trait_id, outcome_family)])
source_rows <- vector("list", 2L * nrow(pair_keys))
qa_rows <- vector("list", nrow(pair_keys))
source_pos <- 0L

for (i in seq_len(nrow(pair_keys))) {
  key <- pair_keys[i]
  idx <- source_index[key, on = .(trait_id, outcome_family)]
  fg <- idx[independent_source_group == "finngen"]
  pu <- idx[independent_source_group == "panukb"]
  if (nrow(fg) != 1L || nrow(pu) != 1L) {
    stop("Expected one FinnGen and one Pan-UKB row for ", key$trait_id, " / ", key$outcome_family)
  }

  d_fg <- prepare_harmonized(fg$harmonized_path)
  d_pu <- prepare_harmonized(pu$harmonized_path)
  raw_common <- intersect(d_fg$SNP, d_pu$SNP)

  align <- merge(
    d_fg[SNP %chin% raw_common, .(
      SNP,
      effect_allele_fg = effect_allele.exposure,
      other_allele_fg = other_allele.exposure,
      beta_exposure_fg = beta.exposure,
      se_exposure_fg = se.exposure
    )],
    d_pu[SNP %chin% raw_common, .(
      SNP,
      effect_allele_pu = effect_allele.exposure,
      other_allele_pu = other_allele.exposure,
      beta_exposure_pu = beta.exposure,
      se_exposure_pu = se.exposure
    )],
    by = "SNP",
    all = FALSE
  )
  align[, allele_match := effect_allele_fg == effect_allele_pu &
    other_allele_fg == other_allele_pu]
  align[, exposure_match := abs(beta_exposure_fg - beta_exposure_pu) < 1e-12 &
    abs(se_exposure_fg - se_exposure_pu) < 1e-12]
  valid_common <- align[allele_match & exposure_match, SNP]

  d_fg_common <- d_fg[SNP %chin% valid_common][order(SNP)]
  d_pu_common <- d_pu[SNP %chin% valid_common][order(SNP)]

  for (src in list(
    list(meta = fg, dat = d_fg_common),
    list(meta = pu, dat = d_pu_common)
  )) {
    source_pos <- source_pos + 1L
    est <- run_ivw(src$dat)
    source_rows[[source_pos]] <- cbind(
      src$meta[, .(
        batch_id,
        trait_id,
        trait_label,
        outcome_family,
        endpoint_key,
        source,
        independent_source_group,
        outcome_label,
        harmonized_path
      )],
      est
    )
  }

  union_n <- length(union(d_fg$SNP, d_pu$SNP))
  qa_rows[[i]] <- data.table(
    trait_id = key$trait_id,
    outcome_family = key$outcome_family,
    nsnp_finngen_available = nrow(d_fg),
    nsnp_panukb_available = nrow(d_pu),
    nsnp_raw_intersection = length(raw_common),
    nsnp_common_analysis = length(valid_common),
    nsnp_union = union_n,
    finngen_only_snps = length(setdiff(d_fg$SNP, d_pu$SNP)),
    panukb_only_snps = length(setdiff(d_pu$SNP, d_fg$SNP)),
    exact_instrument_set_match = setequal(d_fg$SNP, d_pu$SNP),
    jaccard_instrument_overlap = if (union_n > 0L) length(raw_common) / union_n else NA_real_,
    allele_mismatch_in_intersection = sum(!align$allele_match),
    exposure_weight_mismatch_in_intersection = sum(!align$exposure_match)
  )

  if (i %% 50L == 0L || i == nrow(pair_keys)) {
    message("Processed ", i, "/", nrow(pair_keys), " paired comparisons")
  }
}

common_source <- rbindlist(source_rows, fill = TRUE)
pair_qa <- rbindlist(qa_rows, fill = TRUE)

wide <- dcast(
  common_source,
  trait_id + trait_label + outcome_family ~ independent_source_group,
  value.var = c(
    "b_common", "se_common", "pval_common", "nsnp_common",
    "Q_common", "Q_df_common", "Q_pval_common"
  )
)
wide <- merge(
  wide,
  pair_qa,
  by = c("trait_id", "outcome_family"),
  all.x = TRUE
)
wide <- merge(
  original[, .(
    trait_id,
    exposure_family,
    outcome_family,
    original_b_finngen = b_finngen,
    original_b_panukb = b_panukb,
    original_se_finngen = se_finngen,
    original_se_panukb = se_panukb,
    original_pval_finngen = pval_finngen,
    original_pval_panukb = pval_panukb,
    original_nsnp_finngen = nsnp_finngen,
    original_nsnp_panukb = nsnp_panukb,
    original_effect_difference = effect_difference,
    original_effect_difference_p = effect_difference_p,
    original_effect_difference_fdr = effect_difference_fdr,
    original_equivalent_or_margin_20 = equivalent_or_margin_20,
    original_both_nominal_concordant = both_nominal_concordant,
    original_asymmetric_nominal_support = asymmetric_nominal_support,
    original_uncertainty_class = uncertainty_class
  )],
  wide,
  by = c("trait_id", "outcome_family"),
  all.y = TRUE
)

wide[, complete_pair_common := is.finite(b_common_finngen) &
  is.finite(b_common_panukb) &
  is.finite(se_common_finngen) &
  is.finite(se_common_panukb)]
wide[, `:=`(
  effect_difference_common = b_common_finngen - b_common_panukb,
  effect_difference_se_common = sqrt(se_common_finngen^2 + se_common_panukb^2)
)]
wide[, effect_difference_z_common := effect_difference_common / effect_difference_se_common]
wide[, effect_difference_p_common := 2 * pnorm(-abs(effect_difference_z_common))]
wide[, effect_difference_fdr_common := NA_real_]
wide[
  is.finite(effect_difference_p_common),
  effect_difference_fdr_common := p.adjust(effect_difference_p_common, method = "BH")
]
wide[, `:=`(
  effect_difference_ci90_low_common =
    effect_difference_common - qnorm(0.95) * effect_difference_se_common,
  effect_difference_ci90_high_common =
    effect_difference_common + qnorm(0.95) * effect_difference_se_common
)]

margin20 <- log(1.20)
wide[, equivalent_or_margin_20_common := complete_pair_common &
  effect_difference_ci90_low_common >= -margin20 &
  effect_difference_ci90_high_common <= margin20]
wide[, both_nominal_concordant_common := complete_pair_common &
  pval_common_finngen < 0.05 &
  pval_common_panukb < 0.05 &
  sign(b_common_finngen) == sign(b_common_panukb)]
wide[, asymmetric_nominal_support_common := complete_pair_common &
  xor(pval_common_finngen < 0.05, pval_common_panukb < 0.05)]
wide[, uncertainty_class_common := fcase(
  !complete_pair_common, "non_estimable",
  effect_difference_fdr_common < 0.05, "evidence_of_source_heterogeneity",
  both_nominal_concordant_common, "concordant_cross_source_support",
  equivalent_or_margin_20_common, "evidence_of_practical_equivalence",
  asymmetric_nominal_support_common, "asymmetric_support_without_demonstrated_heterogeneity",
  default = "inconclusive_or_no_evidence_in_either_source"
)]
wide[, `:=`(
  class_changed_after_common_snp = original_uncertainty_class != uncertainty_class_common,
  heterogeneity_changed_after_common_snp =
    (original_uncertainty_class == "evidence_of_source_heterogeneity") !=
    (uncertainty_class_common == "evidence_of_source_heterogeneity"),
  effect_difference_shift_common_minus_original =
    effect_difference_common - original_effect_difference,
  finngen_beta_shift_common_minus_original =
    b_common_finngen - original_b_finngen,
  panukb_beta_shift_common_minus_original =
    b_common_panukb - original_b_panukb
)]

class_levels <- c(
  "evidence_of_source_heterogeneity",
  "concordant_cross_source_support",
  "evidence_of_practical_equivalence",
  "asymmetric_support_without_demonstrated_heterogeneity",
  "inconclusive_or_no_evidence_in_either_source",
  "non_estimable"
)

class_summary <- wide[, .(
  n = .N,
  percent = 100 * .N / nrow(wide)
), by = uncertainty_class_common][order(match(uncertainty_class_common, class_levels))]

class_transition <- dcast(
  wide[, .N, by = .(original_uncertainty_class, uncertainty_class_common)],
  original_uncertainty_class ~ uncertainty_class_common,
  value.var = "N",
  fill = 0
)

naive_reclassification <- wide[
  asymmetric_nominal_support_common == TRUE,
  .(n = .N),
  by = uncertainty_class_common
][order(match(uncertainty_class_common, class_levels))]

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  cor(x[ok], y[ok])
}

metrics <- data.table(
  metric = c(
    "paired_comparisons",
    "pairs_with_different_original_snp_counts",
    "pairs_with_exact_instrument_set_match",
    "common_snp_estimable_pairs",
    "median_common_snp_count",
    "minimum_common_snp_count",
    "maximum_common_snp_count",
    "original_source_heterogeneity_pairs",
    "common_snp_source_heterogeneity_pairs",
    "original_heterogeneity_retained_with_common_snps",
    "new_common_snp_heterogeneity_pairs",
    "classification_agreement_pairs",
    "classification_agreement_percent",
    "original_naive_source_restricted_pairs",
    "common_snp_naive_source_restricted_pairs",
    "effect_difference_correlation",
    "median_absolute_effect_difference_shift",
    "maximum_allele_mismatches_per_pair",
    "maximum_exposure_weight_mismatches_per_pair"
  ),
  value = c(
    nrow(wide),
    sum(wide$original_nsnp_finngen != wide$original_nsnp_panukb, na.rm = TRUE),
    sum(wide$exact_instrument_set_match, na.rm = TRUE),
    sum(wide$complete_pair_common, na.rm = TRUE),
    median(wide$nsnp_common_analysis, na.rm = TRUE),
    min(wide$nsnp_common_analysis, na.rm = TRUE),
    max(wide$nsnp_common_analysis, na.rm = TRUE),
    sum(wide$original_uncertainty_class == "evidence_of_source_heterogeneity", na.rm = TRUE),
    sum(wide$uncertainty_class_common == "evidence_of_source_heterogeneity", na.rm = TRUE),
    sum(
      wide$original_uncertainty_class == "evidence_of_source_heterogeneity" &
        wide$uncertainty_class_common == "evidence_of_source_heterogeneity",
      na.rm = TRUE
    ),
    sum(
      wide$original_uncertainty_class != "evidence_of_source_heterogeneity" &
        wide$uncertainty_class_common == "evidence_of_source_heterogeneity",
      na.rm = TRUE
    ),
    sum(!wide$class_changed_after_common_snp, na.rm = TRUE),
    100 * mean(!wide$class_changed_after_common_snp, na.rm = TRUE),
    sum(wide$original_asymmetric_nominal_support, na.rm = TRUE),
    sum(wide$asymmetric_nominal_support_common, na.rm = TRUE),
    safe_cor(wide$original_effect_difference, wide$effect_difference_common),
    median(abs(wide$effect_difference_shift_common_minus_original), na.rm = TRUE),
    max(wide$allele_mismatch_in_intersection, na.rm = TRUE),
    max(wide$exposure_weight_mismatch_in_intersection, na.rm = TRUE)
  )
)

metric_value <- function(name) metrics[metric == name, value]

fwrite(source_index, file.path(out_dir, "common_snp_source_file_index.tsv"), sep = "\t")
fwrite(pair_qa, file.path(out_dir, "common_snp_pair_instrument_overlap_qa.tsv"), sep = "\t")
fwrite(common_source, file.path(out_dir, "common_snp_source_specific_ivw.tsv"), sep = "\t")
fwrite(wide, file.path(out_dir, "common_snp_pair_classification.tsv"), sep = "\t")
fwrite(class_summary, file.path(out_dir, "common_snp_class_summary.tsv"), sep = "\t")
fwrite(class_transition, file.path(out_dir, "common_snp_class_transition.tsv"), sep = "\t")
fwrite(naive_reclassification, file.path(out_dir, "common_snp_naive_reclassification.tsv"), sep = "\t")
fwrite(metrics, file.path(out_dir, "common_snp_agreement_metrics.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_common_snp.txt"))

report <- c(
  "# Common-Instrument Calibration Reanalysis",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Design",
  "",
  "For each of the 1,060 exposure-endpoint pairs, FinnGen and Pan-UKB IVW estimates were recomputed using only SNPs that survived harmonization in both sources. Exposure alleles, exposure effects, and exposure standard errors were required to match exactly across the two harmonized files. The paired source-effect contrast, BH-FDR, 20% odds-ratio compatibility margin, and interpretation classes were then recalculated using the calibration-analysis rules.",
  "",
  "The common-SNP result isolates outcome-association differences more cleanly from instrument availability and composition. The original source-specific-instrument analysis remains a real-world portability analysis because it reflects the instruments actually estimable in each public source.",
  "",
  "## Instrument Overlap",
  "",
  paste0("- Paired comparisons: ", metric_value("paired_comparisons"), "."),
  paste0("- Pairs with different original SNP counts: ", metric_value("pairs_with_different_original_snp_counts"), "."),
  paste0("- Pairs with exactly matching instrument sets: ", metric_value("pairs_with_exact_instrument_set_match"), "."),
  paste0("- Common-SNP estimable pairs: ", metric_value("common_snp_estimable_pairs"), "."),
  paste0("- Median common-SNP count: ", sprintf("%.1f", metric_value("median_common_snp_count")), "."),
  paste0("- Minimum common-SNP count: ", metric_value("minimum_common_snp_count"), "."),
  paste0("- Maximum common-SNP count: ", metric_value("maximum_common_snp_count"), "."),
  paste0("- Maximum allele mismatches within an SNP intersection: ", metric_value("maximum_allele_mismatches_per_pair"), "."),
  paste0("- Maximum exposure-weight mismatches within an SNP intersection: ", metric_value("maximum_exposure_weight_mismatches_per_pair"), "."),
  "",
  "## Primary Comparison",
  "",
  paste0("- Original source-heterogeneity pairs: ", metric_value("original_source_heterogeneity_pairs"), "."),
  paste0("- Common-SNP source-heterogeneity pairs: ", metric_value("common_snp_source_heterogeneity_pairs"), "."),
  paste0("- Original heterogeneity pairs retained with common SNPs: ", metric_value("original_heterogeneity_retained_with_common_snps"), "."),
  paste0("- New common-SNP heterogeneity pairs: ", metric_value("new_common_snp_heterogeneity_pairs"), "."),
  paste0("- Exact interpretation-class agreement: ", metric_value("classification_agreement_pairs"), "/", metric_value("paired_comparisons"), " (", sprintf("%.2f", metric_value("classification_agreement_percent")), "%)."),
  paste0("- Correlation between original and common-SNP effect differences: ", sprintf("%.4f", metric_value("effect_difference_correlation")), "."),
  paste0("- Median absolute change in the effect difference: ", sprintf("%.5f", metric_value("median_absolute_effect_difference_shift")), "."),
  "",
  "## Significance-Asymmetry Comparison",
  "",
  paste0("- Original naive one-source-significance pairs: ", metric_value("original_naive_source_restricted_pairs"), "."),
  paste0("- Common-SNP naive one-source-significance pairs: ", metric_value("common_snp_naive_source_restricted_pairs"), "."),
  "",
  "Common-SNP naive-claim reclassification:",
  "",
  paste(capture.output(print(naive_reclassification)), collapse = "\n"),
  "",
  "## Common-SNP Class Summary",
  "",
  paste(capture.output(print(class_summary)), collapse = "\n"),
  "",
  "## Interpretation",
  "",
  "The reporting output distinguishes two estimands: a common-instrument source contrast, which holds the instrument set fixed, and a source-specific-instrument portability contrast, which includes public-source variant availability. The former is the cleaner analysis for source-effect language; the latter measures the practical consequence of analyzing each source with all available harmonized instruments.",
  "",
  "## Outputs",
  "",
  "- `common_snp_pair_classification.tsv`",
  "- `common_snp_source_specific_ivw.tsv`",
  "- `common_snp_pair_instrument_overlap_qa.tsv`",
  "- `common_snp_class_transition.tsv`",
  "- `common_snp_naive_reclassification.tsv`",
  "- `common_snp_agreement_metrics.tsv`",
  "- `sessionInfo_common_snp.txt`"
)
writeLines(report, file.path(out_dir, "COMMON_SNP_REANALYSIS_REPORT.md"))

message("Wrote common-SNP reanalysis outputs to ", out_dir)
