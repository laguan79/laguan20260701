#!/usr/bin/env Rscript

usage <- function() {
  cat(paste(
    "Usage:",
    "  Rscript run_precision_asymmetry_simulation.R --input TableS2.tsv --outdir simulation",
    "       [--n-families 200] [--n-tests 1060]",
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

safe_p <- function(z) {
  2 * pnorm(-abs(z))
}

simulate_scenario <- function(se_a_obs, se_b_obs, common_effect, true_source_difference,
                              rho, precision_regime, n_families, n_tests,
                              alpha = 0.05, or_margin = 1.20) {
  n_obs <- length(se_a_obs)
  counts <- c(
    naive_asymmetric = 0,
    nominal_heterogeneity = 0,
    fdr_heterogeneity = 0,
    practical_equivalence = 0,
    both_nominal_concordant = 0,
    total = 0
  )

  for (family_id in seq_len(n_families)) {
    idx <- sample.int(n_obs, n_tests, replace = TRUE)
    se_a <- se_a_obs[idx]
    se_b <- se_b_obs[idx]

    if (precision_regime == "equalized_pair_mean") {
      pooled <- sqrt((se_a^2 + se_b^2) / 2)
      se_a <- pooled
      se_b <- pooled
    } else if (precision_regime == "source_b_50pct_less_precise") {
      se_b <- se_b * 1.5
    } else if (precision_regime == "source_b_50pct_more_precise") {
      se_b <- se_b * 0.67
    } else if (precision_regime != "observed_se_pairs") {
      stop("Unknown precision regime: ", precision_regime)
    }

    true_beta_a <- common_effect + true_source_difference / 2
    true_beta_b <- common_effect - true_source_difference / 2

    z_a <- rnorm(n_tests)
    z_ind <- rnorm(n_tests)
    z_b <- rho * z_a + sqrt(1 - rho^2) * z_ind

    beta_a <- true_beta_a + se_a * z_a
    beta_b <- true_beta_b + se_b * z_b
    p_a <- safe_p(beta_a / se_a)
    p_b <- safe_p(beta_b / se_b)

    diff <- beta_a - beta_b
    diff_se <- sqrt(se_a^2 + se_b^2 - 2 * rho * se_a * se_b)
    diff_p <- safe_p(diff / diff_se)
    diff_fdr <- p.adjust(diff_p, method = "BH")

    ci90_low <- diff - qnorm(0.95) * diff_se
    ci90_high <- diff + qnorm(0.95) * diff_se
    margin <- log(or_margin)
    practical_equivalence <- ci90_low >= -margin & ci90_high <= margin
    same_direction <- sign(beta_a) == sign(beta_b) & sign(beta_a) != 0
    both_nominal_concordant <- p_a < alpha & p_b < alpha & same_direction

    counts["naive_asymmetric"] <- counts["naive_asymmetric"] + sum(xor(p_a < alpha, p_b < alpha))
    counts["nominal_heterogeneity"] <- counts["nominal_heterogeneity"] + sum(diff_p < alpha)
    counts["fdr_heterogeneity"] <- counts["fdr_heterogeneity"] + sum(diff_fdr < alpha)
    counts["practical_equivalence"] <- counts["practical_equivalence"] + sum(practical_equivalence)
    counts["both_nominal_concordant"] <- counts["both_nominal_concordant"] + sum(both_nominal_concordant)
    counts["total"] <- counts["total"] + n_tests
  }

  total <- counts["total"]
  data.frame(
    common_effect = common_effect,
    common_or = exp(common_effect),
    true_source_difference = true_source_difference,
    true_source_difference_or_ratio = exp(true_source_difference),
    rho = rho,
    precision_regime = precision_regime,
    n_families = n_families,
    n_tests_per_family = n_tests,
    n_total_tests = total,
    naive_asymmetric_support_rate = counts["naive_asymmetric"] / total,
    nominal_heterogeneity_rate = counts["nominal_heterogeneity"] / total,
    fdr_heterogeneity_rate = counts["fdr_heterogeneity"] / total,
    practical_equivalence_rate = counts["practical_equivalence"] / total,
    both_nominal_concordant_rate = counts["both_nominal_concordant"] / total,
    false_source_restriction_rate = ifelse(true_source_difference == 0, counts["naive_asymmetric"] / total, NA_real_),
    stringsAsFactors = FALSE
  )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (is.null(args$input) || is.null(args$outdir)) {
    usage()
    quit(status = 1)
  }

  n_families <- if (!is.null(args[["n-families"]])) as.integer(args[["n-families"]]) else 200L
  n_tests <- if (!is.null(args[["n-tests"]])) as.integer(args[["n-tests"]]) else 1060L
  if (!is.finite(n_families) || n_families < 10) stop("--n-families must be at least 10")
  if (!is.finite(n_tests) || n_tests < 100) stop("--n-tests must be at least 100")

  set.seed(20260630)
  dat <- read.delim(args$input, check.names = FALSE, stringsAsFactors = FALSE)
  se_a <- suppressWarnings(as.numeric(dat$se_finngen))
  se_b <- suppressWarnings(as.numeric(dat$se_panukb))
  keep <- is.finite(se_a) & is.finite(se_b) & se_a > 0 & se_b > 0
  se_a <- se_a[keep]
  se_b <- se_b[keep]
  if (length(se_a) < 100) stop("Too few valid SE pairs")

  common_effects <- log(c(1.00, 1.10, 1.20))
  true_differences <- log(c(1.00, 1.10, 1.20, 1.30))
  rhos <- c(0, 0.25, 0.50)
  precision_regimes <- c(
    "equalized_pair_mean",
    "observed_se_pairs",
    "source_b_50pct_less_precise",
    "source_b_50pct_more_precise"
  )

  grid <- expand.grid(
    common_effect = common_effects,
    true_source_difference = true_differences,
    rho = rhos,
    precision_regime = precision_regimes,
    stringsAsFactors = FALSE
  )

  dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
  rows <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    rows[[i]] <- simulate_scenario(
      se_a_obs = se_a,
      se_b_obs = se_b,
      common_effect = grid$common_effect[[i]],
      true_source_difference = grid$true_source_difference[[i]],
      rho = grid$rho[[i]],
      precision_regime = grid$precision_regime[[i]],
      n_families = n_families,
      n_tests = n_tests
    )
  }
  out <- do.call(rbind, rows)
  out$scenario_id <- seq_len(nrow(out))
  out <- out[, c("scenario_id", setdiff(names(out), "scenario_id"))]

  full_path <- file.path(args$outdir, "precision_asymmetry_operating_characteristics.tsv")
  write.table(out, full_path, sep = "\t", quote = FALSE, row.names = FALSE)

  key <- out[out$true_source_difference == 0 & out$common_effect %in% log(c(1.00, 1.10, 1.20)), ]
  key <- key[, c(
    "common_or", "rho", "precision_regime",
    "naive_asymmetric_support_rate", "nominal_heterogeneity_rate",
    "fdr_heterogeneity_rate", "practical_equivalence_rate"
  )]
  key_path <- file.path(args$outdir, "precision_asymmetry_null_key_results.tsv")
  write.table(key, key_path, sep = "\t", quote = FALSE, row.names = FALSE)

  sink(file.path(args$outdir, "precision_asymmetry_simulation_sessionInfo.txt"))
  print(sessionInfo())
  sink()

  message("Wrote operating characteristics: ", full_path)
  message("Wrote null key results: ", key_path)
}

if (sys.nframe() == 0) {
  main()
}
