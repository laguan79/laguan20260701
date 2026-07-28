#!/usr/bin/env Rscript

sensitivity_usage <- function() {
  cat(paste(
    "Usage:",
    "  Rscript run_margin_rho_sensitivity.R --input input.tsv --outdir sensitivity_outputs",
    "       [--margins 1.10,1.15,1.20,1.30] [--rhos 0,0.25,0.50]",
    sep = "\n"
  ))
}

parse_sensitivity_args <- function(args) {
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

num_list <- function(x, default) {
  if (is.null(x)) x <- default
  vals <- suppressWarnings(as.numeric(strsplit(x, ",", fixed = TRUE)[[1]]))
  if (any(!is.finite(vals))) stop("Could not parse numeric list: ", x)
  vals
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

cmd <- commandArgs(FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[[1]]) else NA_character_
script_dir <- dirname(normalizePath(script_path %||% getwd(), mustWork = FALSE))
if (!dir.exists(script_dir)) script_dir <- getwd()
source(file.path(script_dir, "source_effect_uncertainty_classifier.R"))

main <- function() {
  args <- parse_sensitivity_args(commandArgs(trailingOnly = TRUE))
  if (is.null(args$input) || is.null(args$outdir)) {
    sensitivity_usage()
    quit(status = 1)
  }

  margins <- num_list(args$margins, "1.10,1.15,1.20,1.30")
  rhos <- num_list(args$rhos, "0,0.25,0.50")
  if (any(margins <= 1)) stop("All margins must be greater than 1")
  if (any(rhos <= -1 | rhos >= 1)) stop("All rho values must be greater than -1 and less than 1")

  dat <- read.delim(args$input, check.names = FALSE, stringsAsFactors = FALSE)
  dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

  reference <- classify_rows(dat, or_margin = 1.20, rho = 0)
  rows <- list()
  idx <- 1
  for (margin in margins) {
    for (rho in rhos) {
      out <- classify_rows(dat, or_margin = margin, rho = rho)
      counts <- as.data.frame(table(out$uncertainty_class), stringsAsFactors = FALSE)
      names(counts) <- c("uncertainty_class", "n")
      counts$or_margin <- margin
      counts$rho <- rho
      counts$n_total <- nrow(out)
      counts$n_changed_from_margin_1_20_rho_0 <- sum(out$uncertainty_class != reference$uncertainty_class)
      rows[[idx]] <- counts[, c("or_margin", "rho", "uncertainty_class", "n", "n_total", "n_changed_from_margin_1_20_rho_0")]
      idx <- idx + 1
    }
  }

  summary <- do.call(rbind, rows)
  write.table(summary, file.path(args$outdir, "margin_rho_sensitivity_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote sensitivity summary: ", file.path(args$outdir, "margin_rho_sensitivity_summary.tsv"))
}

if (sys.nframe() == 0) {
  main()
}
