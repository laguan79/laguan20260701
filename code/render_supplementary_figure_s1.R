#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

script_arg <- commandArgs(FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
if (is.na(script_file) || !nzchar(script_file)) {
  script_file <- file.path(getwd(), "code", "render_supplementary_figure_s1.R")
}
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
root <- dirname(script_dir)
source_dir <- file.path(root, "figures", "source_data")
output_dir <- file.path(root, "supplementary_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

paths <- c(
  overlap = file.path(source_dir, "supplementary_figure_s1_panelA_common_snp_overlap.tsv"),
  transition = file.path(source_dir, "supplementary_figure_s1_panelB_class_transition_long.tsv"),
  overcall = file.path(source_dir, "supplementary_figure_s1_panelC_naive_overcall_stability.tsv"),
  effect = file.path(source_dir, "supplementary_figure_s1_panelD_effect_difference_stability.tsv")
)
if (any(!file.exists(paths))) {
  stop("Missing Supplementary Figure S1 inputs: ", paste(paths[!file.exists(paths)], collapse = "; "))
}

overlap <- fread(paths[["overlap"]])
transition <- fread(paths[["transition"]])
overcall <- fread(paths[["overcall"]])
effect <- fread(paths[["effect"]])

col_heterogeneity <- "#6B5CA5"
col_equivalence <- "#3B9B79"
col_asymmetry <- "#C36A5A"
col_concordant <- "#4F87B7"
col_inconclusive <- "#B7C1C8"
col_dark <- "#253744"

theme_supp <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", color = col_dark, size = base_size + 1),
      plot.subtitle = element_text(color = "#617380", size = base_size),
      axis.title = element_text(color = col_dark),
      panel.grid.minor = element_blank(),
      legend.title = element_blank()
    )
}

panel_a <- ggplot(overlap, aes(x = nsnp_common_analysis)) +
  geom_histogram(binwidth = 10, boundary = 0, fill = col_concordant, color = "white") +
  geom_vline(xintercept = median(overlap$nsnp_common_analysis), linetype = 2, color = col_dark) +
  annotate(
    "label", x = Inf, y = Inf, hjust = 1.02, vjust = 1.15,
    label = paste0(
      "Different original SNP counts: ", sum(overlap$original_nsnp_finngen != overlap$original_nsnp_panukb), "/1,060\n",
      "Different original SNP sets: ", sum(!overlap$exact_instrument_set_match), "/1,060\n",
      "Common-SNP estimable: ", nrow(overlap), "/1,060\n",
      "Exact original sets: ", sum(overlap$exact_instrument_set_match), "/1,060"
    ),
    size = 3.2, color = col_dark
  ) +
  labs(
    title = "Common-SNP reanalysis retained all paired comparisons",
    subtitle = "Distribution of shared instruments per FinnGen/Pan-UKB source pair",
    x = "Common SNPs retained", y = "Paired comparisons"
  ) +
  theme_supp()

class_order <- c(
  "Source\nheterogeneity", "Concordant\nsupport", "Practical\nequivalence",
  "Asymmetric\nsupport", "Inconclusive/\nno evidence"
)
transition[, original_label := factor(gsub('"', "", original_label), levels = class_order)]
transition[, common_label := factor(gsub('"', "", common_label), levels = class_order)]
panel_b <- ggplot(transition, aes(x = common_label, y = original_label, fill = n)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(n > 0, n, "")), color = col_dark, size = 3.4) +
  scale_fill_gradient(low = "#F2F5F7", high = "#2F7696") +
  labs(
    title = "Interpretation classes were stable after fixing instruments",
    subtitle = "Exact class agreement 986/1,060 (93.0%)",
    x = "Common-SNP class", y = "Source-specific-instrument class"
  ) +
  theme_supp(9) +
  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8))

class_colors <- c(
  "Source\nheterogeneity" = col_heterogeneity,
  "Practical\nequivalence" = col_equivalence,
  "Asymmetric\nsupport" = col_asymmetry
)
overcall[, analysis := factor(analysis, levels = unique(analysis))]
overcall[, class_label := factor(class_label, levels = names(class_colors))]
panel_c <- ggplot(overcall, aes(x = analysis, y = n, fill = class_label)) +
  geom_col(width = 0.56, color = "white") +
  geom_text(
    aes(label = n), position = position_stack(vjust = 0.5),
    color = "white", fontface = "bold", size = 3.3
  ) +
  geom_text(
    data = unique(overcall[, .(analysis, total)]),
    aes(x = analysis, y = total + 8, label = paste0("n=", total)),
    inherit.aes = FALSE, color = col_dark, size = 3.2
  ) +
  scale_fill_manual(values = class_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "One-source significance remained mostly overcalled",
    subtitle = "Naive source-restricted rows after direct paired contrast",
    x = NULL, y = "Naive one-source-significance rows"
  ) +
  theme_supp() +
  theme(legend.position = "bottom")

transition_colors <- c(
  "No FDR heterogeneity" = col_inconclusive,
  "Retained heterogeneity" = col_heterogeneity,
  "Lost after common SNPs" = "#D98C78",
  "New common-SNP heterogeneity" = col_equivalence
)
panel_d <- ggplot(
  effect,
  aes(x = original_effect_difference, y = effect_difference_common, color = transition_status)
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = col_dark) +
  geom_point(alpha = 0.75, size = 1.7) +
  annotate(
    "label", x = -Inf, y = Inf, hjust = -0.02, vjust = 1.1,
    label = "r = 0.9344\nMedian |shift| = 0.01911\nRetained: 33; new: 2",
    size = 3.1, color = col_dark
  ) +
  scale_color_manual(
    values = transition_colors,
    labels = c(
      "No FDR heterogeneity" = "No FDR heterogeneity",
      "Retained heterogeneity" = "Retained",
      "Lost after common SNPs" = "Lost",
      "New common-SNP heterogeneity" = "New"
    )
  ) +
  coord_equal() +
  labs(
    title = "Source-effect differences remained aligned",
    subtitle = "All-available versus common-SNP instruments",
    x = "Original source-effect difference", y = "Common-SNP source-effect difference"
  ) +
  theme_supp() +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    plot.margin = margin(5.5, 12, 5.5, 5.5)
  )

figure <- ((panel_a + panel_b) / (panel_c + panel_d)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.margin = margin(18, 20, 18, 20))
  ) &
  theme(plot.tag = element_text(face = "bold", size = 16, color = col_dark))

stem <- file.path(output_dir, "Supplementary_Figure_S1_Common_SNP_Source_Effect_Stability")
ggsave(paste0(stem, ".png"), figure, width = 13.2, height = 9.4, dpi = 320, bg = "white")
ggsave(paste0(stem, ".pdf"), figure, width = 13.2, height = 9.4, bg = "white")
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(paste0(stem, ".svg"), figure, width = 13.2, height = 9.4, bg = "white")
}
message("Supplementary Figure S1 rendered in ", output_dir)
