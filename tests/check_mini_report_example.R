suppressPackageStartupMessages({
  library(data.table)
})

args0 <- commandArgs(FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else file.path("tests", "check_mini_report_example.R")
root <- normalizePath(file.path(dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE)), ".."), winslash = "/", mustWork = TRUE)

report_path <- file.path(root, "reporting_extension", "example_output_source_effect_report.tsv")
summary_path <- file.path(root, "reporting_extension", "example_output_source_effect_report_class_summary.tsv")

stopifnot(file.exists(report_path))
stopifnot(file.exists(summary_path))

report <- fread(report_path)
summary <- fread(summary_path)

stopifnot(nrow(report) == 3)
stopifnot(sum(summary$n) == 3)
stopifnot("source_effect_class" %in% names(report))
stopifnot("recommended_wording" %in% names(report))

cat("Mini-report example checks passed.\n")
