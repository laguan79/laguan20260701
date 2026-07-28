#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(cowplot)
  library(ggplot2)
})

script_arg <- commandArgs(FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
if (is.na(script_file) || !nzchar(script_file)) script_file <- file.path(getwd(), "render_method_hardening_figure3.R")
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
project_root <- dirname(script_dir)
figure_root <- file.path(project_root, "figures")
report_dir <- file.path(project_root, "reporting_extension")
main_dir <- figure_root
data_dir <- file.path(figure_root, "source_data")
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

input_path <- file.path(report_dir, "example_input_source_effect.tsv")
output_path <- file.path(report_dir, "example_output_source_effect_report.tsv")
if (any(!file.exists(c(input_path, output_path)))) stop("Reporting-extension inputs are missing")

input <- fread(input_path)
output <- fread(output_path)
example <- output[comparison_id == "cigarettes_per_day_gscan__tier1_023_constipation"][1]
if (nrow(example) != 1) stop("Worked example row not found")

input_schema <- data.table(
  group = c("Comparison", "Source A", "Source B", "Context", "Assumptions"),
  fields = c(
    "ID; exposure; outcome",
    "label; beta; SE; P",
    "label; beta; SE; P",
    "endpoint; instrument basis",
    "rho; OR margin; FDR family"
  ),
  returned = c(
    "Traceable unit",
    "Estimate + uncertainty",
    "Estimate + uncertainty",
    "Analysis labels",
    "Decision settings"
  )
)
fwrite(input_schema, file.path(data_dir, "figure3_panelA_input_schema.tsv"), sep = "\t")

col_hetero <- "#6B5CA5"
col_asym <- "#C36A5A"
col_equiv <- "#3B9B79"
col_concord <- "#4F87B7"
col_inconclusive <- "#9AA4AD"
col_dark <- "#20323F"
col_mid <- "#526875"

box <- function(plot, x, y, w, h, fill, color = "#C7D1D7", radius = 0.02) {
  plot + draw_grob(grid::roundrectGrob(
    x = x, y = y, width = w, height = h,
    r = grid::unit(radius, "snpc"),
    gp = grid::gpar(fill = fill, col = color, lwd = 1)
  ))
}

arrow_line <- function(plot, x, y, xend, yend, color = "#6A7E89") {
  plot + draw_grob(grid::segmentsGrob(
    x0 = x, y0 = y, x1 = xend, y1 = yend,
    arrow = grid::arrow(length = grid::unit(0.06, "inches"), type = "closed"),
    gp = grid::gpar(col = color, lwd = 1.5)
  ))
}

# Panel A: structured input schema.
panel_a <- ggdraw() +
  draw_label("Mini-report input", x = 0.12, y = 0.95, hjust = 0, fontface = "bold", size = 9, color = col_dark) +
  draw_label("Field", x = 0.06, y = 0.85, hjust = 0, fontface = "bold", size = 6.2, color = col_mid) +
  draw_label("Required values", x = 0.31, y = 0.85, hjust = 0, fontface = "bold", size = 6.2, color = col_mid) +
  draw_label("Role", x = 0.76, y = 0.85, hjust = 0, fontface = "bold", size = 6.2, color = col_mid)
ys <- seq(0.72, 0.18, length.out = nrow(input_schema))
for (i in seq_len(nrow(input_schema))) {
  panel_a <- box(panel_a, 0.50, ys[i], 0.91, 0.10, if (i %% 2 == 1) "#F5F8FA" else "#EDF3F6") +
    draw_label(input_schema$group[i], x = 0.06, y = ys[i], hjust = 0, size = 6.0, fontface = "bold", color = col_dark) +
    draw_label(input_schema$fields[i], x = 0.31, y = ys[i], hjust = 0, size = 5.5, color = col_dark) +
    draw_label(input_schema$returned[i], x = 0.76, y = ys[i], hjust = 0, size = 5.3, color = col_mid)
}

# Panel B: standard statistic and reporting decision process.
panel_b <- ggdraw() +
  draw_label("Paired decision engine", x = 0.12, y = 0.95, hjust = 0, fontface = "bold", size = 9, color = col_dark)
panel_b <- box(panel_b, 0.17, 0.66, 0.27, 0.33, "#E9F1F7", "#4F87B7") +
  draw_label("Paired estimates", x = 0.17, y = 0.76, fontface = "bold", size = 6.5, color = col_dark) +
  draw_label("beta1, SE1\nbeta2, SE2", x = 0.17, y = 0.62, size = 6.0, color = col_mid)
panel_b <- arrow_line(panel_b, 0.31, 0.66, 0.39, 0.66)
panel_b <- box(panel_b, 0.53, 0.66, 0.27, 0.42, "#F1F3F5", "#7B8B95") +
  draw_label("Direct contrast", x = 0.53, y = 0.80, fontface = "bold", size = 6.5, color = col_dark) +
  draw_label("zDelta = (beta1 - beta2)\n/ SEDelta", x = 0.53, y = 0.68, size = 5.5, color = col_dark) +
  draw_label("rho = 0: Q = zDelta^2", x = 0.53, y = 0.56, size = 5.4, color = col_mid)
panel_b <- arrow_line(panel_b, 0.67, 0.66, 0.75, 0.66)
panel_b <- box(panel_b, 0.86, 0.66, 0.22, 0.42, "#EEF3F6", "#7B8B95") +
  draw_label("Reporting layer", x = 0.86, y = 0.80, fontface = "bold", size = 6.5, color = col_dark) +
  draw_label("BH-FDR family\n90% CI + OR margin\nrho sensitivity\ninstrument basis", x = 0.86, y = 0.63, size = 5.2, color = col_mid, lineheight = 1.12)
panel_b <- arrow_line(panel_b, 0.86, 0.44, 0.86, 0.34)
panel_b <- box(panel_b, 0.62, 0.19, 0.70, 0.19, "#EEF3F6", "#7B8B95") +
  draw_label("Evidence class + wording + session information", x = 0.62, y = 0.19, size = 5.8, fontface = "bold", color = col_dark)

# Panel C: output classes.
classes <- data.table(
  y = c(0.76, 0.63, 0.50, 0.37, 0.24),
  label = c("FDR source heterogeneity", "Practical equivalence", "Asymmetric support", "Concordant support", "Inconclusive"),
  rule = c("paired contrast passes FDR", "90% CI within declared margin", "one-source nominal; contrast unsupported", "both sources nominally supported", "remaining paired evidence"),
  color = c(col_hetero, col_equiv, col_asym, col_concord, col_inconclusive)
)
panel_c <- ggdraw() +
  draw_label("Returned evidence classes", x = 0.12, y = 0.95, hjust = 0, fontface = "bold", size = 9, color = col_dark)
for (i in seq_len(nrow(classes))) {
  panel_c <- box(panel_c, 0.50, classes$y[i], 0.90, 0.105, paste0(classes$color[i], "20"), classes$color[i]) +
    draw_label(classes$label[i], x = 0.08, y = classes$y[i], hjust = 0, fontface = "bold", size = 5.8, color = classes$color[i]) +
    draw_label(classes$rule[i], x = 0.52, y = classes$y[i], hjust = 0, size = 5.2, color = col_dark)
}

# Panel D: measured worked example.
panel_d <- ggdraw() +
  draw_label("Worked example: smoking intensity and constipation", x = 0.12, y = 0.95, hjust = 0, fontface = "bold", size = 8.7, color = col_dark)
panel_d <- box(panel_d, 0.27, 0.69, 0.45, 0.32, "#F5F8FA", "#C8D4DA") +
  draw_label("Source estimates", x = 0.27, y = 0.80, fontface = "bold", size = 6.5, color = col_dark) +
  draw_label(sprintf("FinnGen: beta %+.3f\nSE %.3f; P = %.3f", example$beta_a, example$se_a, example$p_a), x = 0.08, y = 0.69, hjust = 0, size = 5.4, color = col_mid) +
  draw_label(sprintf("Pan-UKB: beta %+.3f\nSE %.3f; P = %.3g", example$beta_b, example$se_b, example$p_b), x = 0.30, y = 0.69, hjust = 0, size = 5.4, color = col_mid)
panel_d <- arrow_line(panel_d, 0.50, 0.69, 0.61, 0.69)
panel_d <- box(panel_d, 0.78, 0.69, 0.34, 0.32, "#F1F3F5", "#7B8B95") +
  draw_label("Paired contrast", x = 0.78, y = 0.80, fontface = "bold", size = 6.5, color = col_dark) +
  draw_label(sprintf("Difference %+.3f", example$effect_difference), x = 0.78, y = 0.70, size = 5.7, color = col_dark) +
  draw_label(sprintf("90%% CI %+.3f to %+.3f", example$effect_difference_ci90_low, example$effect_difference_ci90_high), x = 0.78, y = 0.61, size = 5.2, color = col_mid) +
  draw_label(sprintf("FDR = %.3f", example$effect_difference_fdr), x = 0.78, y = 0.56, size = 5.5, color = col_mid)
panel_d <- box(panel_d, 0.50, 0.30, 0.90, 0.25, "#FFF2EA", col_asym) +
  draw_label("Naive label", x = 0.08, y = 0.36, hjust = 0, fontface = "bold", size = 5.8, color = col_asym) +
  draw_label("Pan-UKB-only / source restricted", x = 0.29, y = 0.36, hjust = 0, size = 5.4, color = col_dark) +
  draw_label("Paired-evidence wording", x = 0.08, y = 0.24, hjust = 0, fontface = "bold", size = 5.8, color = col_asym) +
  draw_label("Asymmetric support; paired contrast FDR = 0.078", x = 0.38, y = 0.24, hjust = 0, size = 5.2, color = col_dark)

fig3 <- plot_grid(
  plot_grid(panel_a, panel_b, labels = c("A", "B"), label_size = 10.5, label_fontface = "bold", nrow = 1, rel_widths = c(1, 1.05)),
  plot_grid(panel_c, panel_d, labels = c("C", "D"), label_size = 10.5, label_fontface = "bold", nrow = 1, rel_widths = c(0.92, 1.08)),
  ncol = 1, rel_heights = c(1, 1)
)

stem <- "Figure3_Reusable_Mini_Report_Worked_Example"
ggsave(file.path(main_dir, paste0(stem, ".png")), fig3, width = 6.69, height = 5.25, dpi = 320, bg = "white")
ggsave(file.path(main_dir, paste0(stem, ".pdf")), fig3, width = 6.69, height = 5.25, bg = "white")
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(file.path(main_dir, paste0(stem, ".svg")), fig3, width = 6.69, height = 5.25, bg = "white")
}
message("Final Figure 3 rendered in ", main_dir)
