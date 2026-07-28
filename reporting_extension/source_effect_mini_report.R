#!/usr/bin/env Rscript

usage <- function() {
  cat(paste(
    "Usage:",
    "  Rscript source_effect_mini_report.R --input input.tsv --output output.tsv [--markdown report.md] [--or-margin 1.20] [--rho 0] [--alpha 0.05] [--fdr-source auto]",
    "",
    "Required input columns:",
    "  comparison_id, beta_a, se_a, p_a, beta_b, se_b, p_b",
    "",
    "Recommended descriptive columns:",
    "  source_a, source_b, exposure, outcome, endpoint_family, instrument_basis",
    "",
    "FDR handling:",
    "  --fdr-source auto     Use an input FDR column if present; otherwise compute BH FDR across input rows.",
    "  --fdr-source input    Require an input FDR column.",
    "  --fdr-source computed Compute BH FDR across input rows.",
    "",
    "Accepted input FDR columns:",
    "  source_difference_fdr, effect_difference_fdr, reported_effect_difference_fdr",
    sep = "\n"
  ))
}

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    key <- sub("^--", "", key)
    if (grepl("=", key, fixed = TRUE)) {
      parts <- strsplit(key, "=", fixed = TRUE)[[1]]
      out[[parts[[1]]]] <- paste(parts[-1], collapse = "=")
      i <- i + 1
    } else {
      if (i == length(args)) stop("Missing value for --", key)
      out[[key]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  out
}

first_present <- function(dat, choices) {
  hit <- choices[choices %in% names(dat)]
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

num_col <- function(dat, choices, required = TRUE) {
  nm <- first_present(dat, choices)
  if (is.null(nm)) {
    if (required) stop("Missing required column; accepted names: ", paste(choices, collapse = ", "))
    return(rep(NA_real_, nrow(dat)))
  }
  suppressWarnings(as.numeric(dat[[nm]]))
}

char_col <- function(dat, choices, default = "") {
  nm <- first_present(dat, choices)
  if (is.null(nm)) return(rep(default, nrow(dat)))
  as.character(dat[[nm]])
}

fmt_p <- function(x) {
  ifelse(
    is.na(x), "NA",
    ifelse(x < 0.001, formatC(x, format = "e", digits = 2), formatC(x, format = "f", digits = 3))
  )
}

classify_source_effect <- function(dat, or_margin = 1.20, rho = 0, alpha = 0.05, fdr_source = "auto") {
  beta_a <- num_col(dat, c("beta_a", "b_a", "b_source_a", "b_finngen"))
  beta_b <- num_col(dat, c("beta_b", "b_b", "b_source_b", "b_panukb"))
  se_a <- num_col(dat, c("se_a", "se_source_a", "se_finngen"))
  se_b <- num_col(dat, c("se_b", "se_source_b", "se_panukb"))
  p_a <- num_col(dat, c("p_a", "pval_a", "p_source_a", "pval_finngen"))
  p_b <- num_col(dat, c("p_b", "pval_b", "p_source_b", "pval_panukb"))

  comparison_id <- char_col(dat, c("comparison_id", "label"), "")
  missing_id <- !nzchar(comparison_id)
  if (any(missing_id)) {
    comparison_id[missing_id] <- paste0("comparison_", seq_len(nrow(dat))[missing_id])
  }

  source_a <- char_col(dat, c("source_a", "source_1"), "Source A")
  source_b <- char_col(dat, c("source_b", "source_2"), "Source B")
  exposure <- char_col(dat, c("exposure", "trait_label", "exposure_label"), "")
  outcome <- char_col(dat, c("outcome", "endpoint", "outcome_label"), "")
  endpoint_family <- char_col(dat, c("endpoint_family", "outcome_family"), "")
  instrument_basis <- char_col(dat, c("instrument_basis", "analysis_basis"), "")

  diff_var <- se_a^2 + se_b^2 - 2 * rho * se_a * se_b
  valid <- is.finite(beta_a) & is.finite(beta_b) &
    is.finite(se_a) & is.finite(se_b) & se_a > 0 & se_b > 0 &
    is.finite(diff_var) & diff_var > 0

  effect_difference <- beta_a - beta_b
  effect_difference_se <- sqrt(diff_var)
  effect_difference_z <- effect_difference / effect_difference_se
  effect_difference_p <- 2 * pnorm(-abs(effect_difference_z))
  effect_difference_p[!valid] <- NA_real_

  fdr_col <- first_present(dat, c("source_difference_fdr", "effect_difference_fdr", "reported_effect_difference_fdr"))
  fdr_source <- match.arg(fdr_source, c("auto", "input", "computed"))
  if (fdr_source == "input" && is.null(fdr_col)) {
    stop("--fdr-source input was requested, but no accepted FDR column was found")
  }
  if (fdr_source %in% c("auto", "input") && !is.null(fdr_col)) {
    effect_difference_fdr <- suppressWarnings(as.numeric(dat[[fdr_col]]))
    fdr_basis <- paste0("input column: ", fdr_col)
  } else {
    effect_difference_fdr <- rep(NA_real_, nrow(dat))
    effect_difference_fdr[is.finite(effect_difference_p)] <- p.adjust(effect_difference_p[is.finite(effect_difference_p)], method = "BH")
    fdr_basis <- "Benjamini-Hochberg FDR computed across input rows"
  }

  ci95_low <- effect_difference - qnorm(0.975) * effect_difference_se
  ci95_high <- effect_difference + qnorm(0.975) * effect_difference_se
  ci90_low <- effect_difference - qnorm(0.95) * effect_difference_se
  ci90_high <- effect_difference + qnorm(0.95) * effect_difference_se
  margin <- log(or_margin)

  practical_equivalence <- valid & ci90_low >= -margin & ci90_high <= margin
  same_direction <- sign(beta_a) == sign(beta_b) & sign(beta_a) != 0
  both_nominal_concordant <- valid & p_a < alpha & p_b < alpha & same_direction
  naive_significance_asymmetry <- valid & ((p_a < alpha & p_b >= alpha) | (p_a >= alpha & p_b < alpha))
  nominal_source_difference <- valid & effect_difference_p < alpha
  fdr_source_difference <- valid & is.finite(effect_difference_fdr) & effect_difference_fdr < alpha

  source_effect_class <- rep("inconclusive_or_no_evidence_in_either_source", nrow(dat))
  source_effect_class[!valid] <- "non_estimable"
  source_effect_class[fdr_source_difference] <- "evidence_of_source_heterogeneity"
  source_effect_class[valid & !fdr_source_difference & both_nominal_concordant] <- "concordant_cross_source_support"
  source_effect_class[valid & !fdr_source_difference & !both_nominal_concordant & practical_equivalence] <-
    "evidence_of_practical_equivalence"
  source_effect_class[valid & !fdr_source_difference & !both_nominal_concordant & !practical_equivalence & naive_significance_asymmetry] <-
    "asymmetric_support_without_demonstrated_heterogeneity"

  naive_claim <- ifelse(
    naive_significance_asymmetry,
    "naive one-source significance would invite source-restricted wording",
    "naive one-source significance is not present"
  )

  recommended_wording <- ifelse(
    source_effect_class == "evidence_of_source_heterogeneity",
    "Report evidence of source heterogeneity, followed by phenotype, ascertainment, instrument, and source-feature review.",
    ifelse(
      source_effect_class == "evidence_of_practical_equivalence",
      "Report practical equivalence under the declared margin; do not describe the result as source restriction or replication failure.",
      ifelse(
        source_effect_class == "concordant_cross_source_support",
        "Report concordant same-direction source support; avoid claiming equivalence unless the declared margin is met.",
        ifelse(
          source_effect_class == "asymmetric_support_without_demonstrated_heterogeneity",
          "Report asymmetric statistical support without demonstrated source heterogeneity; avoid replication-failure or source-specific language.",
          ifelse(
            source_effect_class == "non_estimable",
            "State that the source comparison is not estimable from the supplied source-specific estimates.",
            "Report the source comparison as inconclusive for heterogeneity, equivalence, or concordant support."
          )
        )
      )
    )
  )

  suggested_sentence <- paste0(
    "For ", ifelse(nzchar(exposure), exposure, comparison_id), " versus ",
    ifelse(nzchar(outcome), outcome, endpoint_family), ", ",
    source_a, " beta = ", formatC(beta_a, format = "f", digits = 3),
    " (SE ", formatC(se_a, format = "f", digits = 3), ", P = ", fmt_p(p_a), ") and ",
    source_b, " beta = ", formatC(beta_b, format = "f", digits = 3),
    " (SE ", formatC(se_b, format = "f", digits = 3), ", P = ", fmt_p(p_b), "); ",
    "the paired source-effect difference was ", formatC(effect_difference, format = "f", digits = 3),
    " (90% CI ", formatC(ci90_low, format = "f", digits = 3), " to ",
    formatC(ci90_high, format = "f", digits = 3), "; FDR = ", fmt_p(effect_difference_fdr),
    "), classified as ", source_effect_class, "."
  )

  data.frame(
    comparison_id = comparison_id,
    exposure = exposure,
    outcome = outcome,
    endpoint_family = endpoint_family,
    source_a = source_a,
    source_b = source_b,
    instrument_basis = instrument_basis,
    beta_a = beta_a,
    se_a = se_a,
    p_a = p_a,
    beta_b = beta_b,
    se_b = se_b,
    p_b = p_b,
    effect_difference = effect_difference,
    effect_difference_se = effect_difference_se,
    effect_difference_p = effect_difference_p,
    effect_difference_fdr = effect_difference_fdr,
    fdr_basis = fdr_basis,
    assumed_source_estimate_correlation = rho,
    effect_difference_ci95_low = ci95_low,
    effect_difference_ci95_high = ci95_high,
    effect_difference_ci90_low = ci90_low,
    effect_difference_ci90_high = ci90_high,
    equivalence_or_margin = or_margin,
    practical_equivalence = practical_equivalence,
    nominal_source_difference = nominal_source_difference,
    fdr_source_difference = fdr_source_difference,
    both_nominal_concordant = both_nominal_concordant,
    naive_significance_asymmetry = naive_significance_asymmetry,
    naive_claim = naive_claim,
    source_effect_class = source_effect_class,
    recommended_wording = recommended_wording,
    suggested_sentence = suggested_sentence,
    stringsAsFactors = FALSE
  )
}

write_markdown_report <- function(out, markdown_path, input_path, or_margin, rho, alpha) {
  class_counts <- as.data.frame(table(out$source_effect_class), stringsAsFactors = FALSE)
  names(class_counts) <- c("source_effect_class", "n")

  lines <- c(
    "# Source-Effect Mini-Report",
    "",
    paste0("- Input: `", basename(input_path), "`."),
    paste0("- Rows classified: ", nrow(out), "."),
    paste0("- Practical-equivalence odds-ratio margin: ", or_margin, "."),
    paste0("- Assumed source-estimate correlation rho: ", rho, "."),
    paste0("- Nominal alpha: ", alpha, "."),
    paste0("- FDR basis: ", unique(out$fdr_basis)[1], "."),
    "",
    "## Class Summary",
    "",
    "| Source-effect class | n |",
    "|---|---:|"
  )
  if (nrow(class_counts)) {
    lines <- c(lines, paste0("| ", class_counts$source_effect_class, " | ", class_counts$n, " |"))
  }
  lines <- c(
    lines,
    "",
    "## Recommended Wording",
    "",
    "| Comparison | Class | Recommended wording |",
    "|---|---|---|"
  )
  for (i in seq_len(nrow(out))) {
    label <- out$comparison_id[[i]]
    lines <- c(lines, paste0("| ", label, " | ", out$source_effect_class[[i]], " | ", out$recommended_wording[[i]], " |"))
  }
  dir.create(dirname(markdown_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, markdown_path, useBytes = TRUE)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (is.null(args$input) || is.null(args$output)) {
    usage()
    quit(status = 1)
  }
  or_margin <- if (!is.null(args[["or-margin"]])) as.numeric(args[["or-margin"]]) else 1.20
  rho <- if (!is.null(args$rho)) as.numeric(args$rho) else 0
  alpha <- if (!is.null(args$alpha)) as.numeric(args$alpha) else 0.05
  fdr_source <- if (!is.null(args[["fdr-source"]])) args[["fdr-source"]] else "auto"
  if (!is.finite(or_margin) || or_margin <= 1) stop("--or-margin must be greater than 1")
  if (!is.finite(rho) || rho <= -1 || rho >= 1) stop("--rho must be greater than -1 and less than 1")
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) stop("--alpha must be between 0 and 1")

  dat <- read.delim(args$input, check.names = FALSE, stringsAsFactors = FALSE)
  out <- classify_source_effect(dat, or_margin = or_margin, rho = rho, alpha = alpha, fdr_source = fdr_source)
  dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
  write.table(out, args$output, sep = "\t", quote = FALSE, row.names = FALSE)

  summary_path <- sub("(\\.tsv)?$", "_class_summary.tsv", args$output)
  class_counts <- as.data.frame(table(out$source_effect_class), stringsAsFactors = FALSE)
  names(class_counts) <- c("source_effect_class", "n")
  write.table(class_counts, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)

  if (!is.null(args$markdown)) {
    write_markdown_report(out, args$markdown, args$input, or_margin, rho, alpha)
  }

  session_path <- sub("(\\.tsv)?$", "_sessionInfo.txt", args$output)
  sink(session_path)
  print(sessionInfo())
  sink()

  message("Wrote mini-report table: ", args$output)
  message("Wrote class summary: ", summary_path)
  message("Wrote session info: ", session_path)
}

if (sys.nframe() == 0) {
  main()
}
