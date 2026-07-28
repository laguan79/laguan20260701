#!/usr/bin/env Rscript

usage <- function() {
  cat(paste(
    "Usage:",
    "  Rscript source_effect_uncertainty_classifier.R --input input.tsv --output output.tsv [--or-margin 1.20] [--rho 0]",
    "",
    "Required input columns, using either generic names or the FinnGen/Pan-UKB aliases:",
    "  generic: comparison_id, beta_a, se_a, p_a, beta_b, se_b, p_b",
    "  aliases: b_finngen, se_finngen, pval_finngen, b_panukb, se_panukb, pval_panukb",
    "",
    "--rho is the assumed correlation between paired source-specific log-effect estimates.",
    "Use 0 when estimates come from non-overlapping samples or overlap is unknown.",
    "",
    "Optional columns are carried through: exposure, outcome, trait_id, exposure_short,",
    "outcome_family, clinical_domain, source_a, source_b.",
    sep = "\n"
  ))
}

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    key <- sub("^--", "", key)
    if (grepl("=", key, fixed = TRUE)) {
      parts <- strsplit(key, "=", fixed = TRUE)[[1]]
      out[[parts[[1]]]] <- paste(parts[-1], collapse = "=")
      i <- i + 1
    } else {
      if (i == length(args)) {
        stop("Missing value for --", key)
      }
      out[[key]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  out
}

first_present <- function(x, choices) {
  hit <- choices[choices %in% names(x)]
  if (length(hit) == 0) return(NULL)
  x[[hit[[1]]]]
}

num_col <- function(x, choices, required = TRUE) {
  v <- first_present(x, choices)
  if (is.null(v)) {
    if (required) stop("Missing required column; accepted names: ", paste(choices, collapse = ", "))
    return(rep(NA_real_, nrow(x)))
  }
  suppressWarnings(as.numeric(v))
}

char_col <- function(x, choices, default) {
  v <- first_present(x, choices)
  if (is.null(v)) return(rep(default, nrow(x)))
  as.character(v)
}

safe_p <- function(z) {
  2 * pnorm(-abs(z))
}

classify_rows <- function(dat, or_margin = 1.20, alpha = 0.05, rho = 0) {
  beta_a <- num_col(dat, c("beta_a", "b_a", "b_source_a", "b_finngen"))
  beta_b <- num_col(dat, c("beta_b", "b_b", "b_source_b", "b_panukb"))
  se_a <- num_col(dat, c("se_a", "se_source_a", "se_finngen"))
  se_b <- num_col(dat, c("se_b", "se_source_b", "se_panukb"))
  p_a <- num_col(dat, c("p_a", "pval_a", "p_source_a", "pval_finngen"))
  p_b <- num_col(dat, c("p_b", "pval_b", "p_source_b", "pval_panukb"))

  diff_var <- se_a^2 + se_b^2 - 2 * rho * se_a * se_b
  valid <- is.finite(beta_a) & is.finite(beta_b) &
    is.finite(se_a) & is.finite(se_b) & se_a > 0 & se_b > 0
  valid <- valid & is.finite(diff_var) & diff_var > 0

  diff <- beta_a - beta_b
  diff_se <- sqrt(diff_var)
  diff_z <- diff / diff_se
  diff_p <- safe_p(diff_z)
  diff_p[!valid] <- NA_real_
  diff_fdr <- rep(NA_real_, length(diff_p))
  diff_fdr[is.finite(diff_p)] <- p.adjust(diff_p[is.finite(diff_p)], method = "BH")

  ci95_low <- diff - qnorm(0.975) * diff_se
  ci95_high <- diff + qnorm(0.975) * diff_se
  ci90_low <- diff - qnorm(0.95) * diff_se
  ci90_high <- diff + qnorm(0.95) * diff_se
  for (v in c("ci95_low", "ci95_high", "ci90_low", "ci90_high")) {
    masked_value <- get(v)
    masked_value[!valid] <- NA_real_
    assign(v, masked_value)
  }

  margin <- log(or_margin)
  practical_equivalence <- valid & ci90_low >= -margin & ci90_high <= margin
  same_direction <- sign(beta_a) == sign(beta_b) & sign(beta_a) != 0
  both_nominal_concordant <- valid & p_a < alpha & p_b < alpha & same_direction
  asymmetric_nominal_support <- valid & ((p_a < alpha & p_b >= alpha) | (p_a >= alpha & p_b < alpha))
  source_heterogeneity <- valid & is.finite(diff_fdr) & diff_fdr < alpha

  uncertainty_class <- rep("inconclusive_or_no_evidence_in_either_source", nrow(dat))
  uncertainty_class[!valid] <- "non_estimable"
  uncertainty_class[source_heterogeneity] <- "evidence_of_source_heterogeneity"
  uncertainty_class[valid & !source_heterogeneity & both_nominal_concordant] <- "concordant_cross_source_support"
  uncertainty_class[valid & !source_heterogeneity & !both_nominal_concordant & practical_equivalence] <- "evidence_of_practical_equivalence"
  uncertainty_class[valid & !source_heterogeneity & !both_nominal_concordant & !practical_equivalence & asymmetric_nominal_support] <-
    "asymmetric_support_without_demonstrated_heterogeneity"

  recommended_wording <- ifelse(
    uncertainty_class == "evidence_of_source_heterogeneity",
    "Direct source-effect difference is supported after multiplicity correction; describe as evidence of source heterogeneity.",
    ifelse(
      uncertainty_class == "evidence_of_practical_equivalence",
      "The source-effect difference is contained within the prespecified equivalence margin; describe as practical equivalence on the chosen scale.",
      ifelse(
        uncertainty_class == "concordant_cross_source_support",
        "Both sources show nominal same-direction support; describe as concordant cross-source statistical support.",
        ifelse(
          uncertainty_class == "asymmetric_support_without_demonstrated_heterogeneity",
          "Only one source shows nominal support, but the direct source-effect difference is not multiplicity-supported; describe as asymmetric statistical support, not replication failure.",
          ifelse(
            uncertainty_class == "non_estimable",
            "The pair is not estimable from available source-specific estimates.",
            "No concordant support, practical equivalence, or multiplicity-supported source heterogeneity is demonstrated; describe as inconclusive."
          )
        )
      )
    )
  )

  comparison_id <- char_col(dat, c("comparison_id", "label"), "")
  missing_id <- !nzchar(comparison_id)
  if (any(missing_id)) {
    lhs <- char_col(dat, c("trait_id", "exposure_short", "exposure"), "exposure")
    rhs <- char_col(dat, c("outcome_family", "outcome", "phenocode"), "outcome")
    comparison_id[missing_id] <- paste(lhs[missing_id], rhs[missing_id], sep = "__")
  }

  source_a <- char_col(dat, c("source_a", "source_1"), "source_a")
  source_b <- char_col(dat, c("source_b", "source_2"), "source_b")
  if (all(source_a == "source_a") && "b_finngen" %in% names(dat)) source_a <- rep("FinnGen R12", nrow(dat))
  if (all(source_b == "source_b") && "b_panukb" %in% names(dat)) source_b <- rep("Pan-UKB", nrow(dat))

  keep <- intersect(
    c("trait_id", "exposure_short", "exposure_unit", "outcome_family", "clinical_domain",
      "finngen_phenocode", "panukb_phenocode", "construct_equivalence_tier"),
    names(dat)
  )

  cbind(
    data.frame(
      comparison_id = comparison_id,
      source_a = source_a,
      source_b = source_b,
      beta_a = beta_a,
      se_a = se_a,
      p_a = p_a,
      beta_b = beta_b,
      se_b = se_b,
      p_b = p_b,
      effect_difference = diff,
      effect_difference_variance = diff_var,
      effect_difference_se = diff_se,
      assumed_source_estimate_correlation = rho,
      source_difference_model = "normal_approximation_for_paired_source_log_effect_difference",
      effect_difference_z = diff_z,
      effect_difference_p = diff_p,
      effect_difference_fdr = diff_fdr,
      effect_difference_ci95_low = ci95_low,
      effect_difference_ci95_high = ci95_high,
      effect_difference_ci90_low = ci90_low,
      effect_difference_ci90_high = ci90_high,
      equivalence_or_margin = or_margin,
      equivalent_or_margin = practical_equivalence,
      both_nominal_concordant = both_nominal_concordant,
      asymmetric_nominal_support = asymmetric_nominal_support,
      uncertainty_class = uncertainty_class,
      recommended_wording = recommended_wording,
      stringsAsFactors = FALSE
    ),
    dat[keep],
    stringsAsFactors = FALSE
  )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (is.null(args$input) || is.null(args$output)) {
    usage()
    quit(status = 1)
  }
  or_margin <- if (!is.null(args[["or-margin"]])) as.numeric(args[["or-margin"]]) else 1.20
  if (!is.finite(or_margin) || or_margin <= 1) {
    stop("--or-margin must be a finite value greater than 1")
  }
  rho <- if (!is.null(args[["rho"]])) as.numeric(args[["rho"]]) else 0
  if (!is.finite(rho) || rho <= -1 || rho >= 1) {
    stop("--rho must be a finite value greater than -1 and less than 1")
  }
  dat <- read.delim(args$input, check.names = FALSE, stringsAsFactors = FALSE)
  out <- classify_rows(dat, or_margin = or_margin, rho = rho)
  dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
  write.table(out, args$output, sep = "\t", quote = FALSE, row.names = FALSE)

  summary <- as.data.frame(table(out$uncertainty_class), stringsAsFactors = FALSE)
  names(summary) <- c("uncertainty_class", "n")
  summary$percent <- 100 * summary$n / nrow(out)
  summary_path <- sub("(\\.tsv)?$", "_summary.tsv", args$output)
  write.table(summary, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote classified rows: ", args$output)
  message("Wrote summary: ", summary_path)
}

if (sys.nframe() == 0) {
  main()
}
