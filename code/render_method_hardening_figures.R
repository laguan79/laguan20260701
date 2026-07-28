#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
mode <- if ("--final" %in% args) "final" else "pilot"
script_arg <- commandArgs(FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
if (is.na(script_file) || !nzchar(script_file)) script_file <- file.path(getwd(), "render_method_hardening_figures.R")
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
package_root <- dirname(script_dir)
figure_root <- file.path(package_root, "figures")
data_dir <- file.path(figure_root, "source_data")
main_dir <- figure_root
supp_dir <- file.path(package_root, "supplementary_figures")
pilot_dir <- file.path(figure_root, "pilot")
portal_dir <- package_root
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pilot_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(
  file.path(data_dir, "figure1_panelB_claim_overcall_decision_impact.tsv"),
  file.path(data_dir, "figure1_panelC_common_snp_stability.tsv"),
  file.path(data_dir, "figure2_source_interpretation_impact_index.tsv"),
  file.path(data_dir, "figure2_supporting_third_source_stress_test.tsv"),
  file.path(data_dir, "figure2_panelD_non_mr_2d_matrix.tsv"),
  file.path(data_dir, "figure2_panelD_non_mr_2d_cell_summary.tsv"),
  file.path(data_dir, "figure2_panelB_axis_counts_rho.tsv"),
  file.path(data_dir, "figure3_simulation_operating_characteristics.tsv"),
  file.path(data_dir, "figure2_panelA_rho_sensitivity.tsv")
)
if (any(!file.exists(required))) stop("Missing figure inputs: ", paste(required[!file.exists(required)], collapse = "; "))

decision <- fread(required[1])
common <- fread(required[2])
impact <- fread(required[3])
third <- fread(required[4])
nonmr_matrix <- fread(required[5])
nonmr_summary <- fread(required[6])
rho_axes <- fread(required[7])
sim <- fread(required[8])
rho_summary <- fread(required[9])

col_hetero <- "#6B5CA5"
col_asym <- "#C36A5A"
col_equiv <- "#3B9B79"
col_concord <- "#4F87B7"
col_inconclusive <- "#9AA4AD"
col_dark <- "#263746"
col_light <- "#E8EEF2"

theme_method <- function(base_size = 10) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1.5, color = "#152533", margin = margin(b = 5, l = 13)),
      plot.subtitle = element_text(size = base_size - 0.2, color = "#586873", margin = margin(b = 5)),
      axis.title = element_text(color = "#263746"),
      axis.text = element_text(color = "#364954"),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(8, 10, 8, 10),
      strip.text = element_text(face = "bold", color = "#263746")
    )
}

tag_theme <- theme(
  plot.tag = element_text(face = "bold", size = 11, color = "#152533"),
  plot.tag.position = c(0.005, 0.995)
)

# Figure 1A: measured study map.
study_nodes <- data.table(
  x = rep(1.8, 4),
  y = c(3.45, 2.55, 1.65, 0.75),
  label = c(
    "20 exposures across 6 families",
    "53 curated endpoint pairs",
    "FinnGen and Pan-UKB estimates",
    "1,060 paired comparisons"
  )
)
p1a <- ggplot() +
  annotate("rect", xmin = study_nodes$x - 1.35, xmax = study_nodes$x + 1.35,
           ymin = study_nodes$y - 0.30, ymax = study_nodes$y + 0.30, fill = c("#E9F1F7", "#EDF4EE", "#F3EEF7", "#FFF2EA"),
           color = c("#4F87B7", "#3B9B79", "#6B5CA5", "#C36A5A"), linewidth = 0.8) +
  annotate("segment", x = rep(1.8, 3), xend = rep(1.8, 3),
           y = c(3.10, 2.20, 1.30), yend = c(2.92, 2.02, 1.12),
           arrow = arrow(length = unit(0.08, "inches")), color = "#607481", linewidth = 0.65) +
  geom_text(data = study_nodes, aes(x, y, label = label), size = 2.65, fontface = "bold", color = col_dark) +
  annotate("text", x = 1.8, y = 0.20, label = "Evidence unit = exposure × endpoint × source pair",
           color = "#465D6B", size = 2.45) +
  coord_cartesian(xlim = c(0.20, 3.40), ylim = c(0.02, 3.90), clip = "off") +
  labs(title = "Calibration design") +
  theme_void(base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 9.5, color = "#152533", margin = margin(b = 4, l = 13)),
        plot.margin = margin(8, 8, 6, 18))

# Figure 1B: all-instrument adjudication.
decision[, display := factor(comparison, levels = rev(c(
  "Naive one-source significance", "Retained as heterogeneity", "Practical equivalence", "Asymmetric support"
)))]
decision[, color_key := fifelse(evidence_class == "evidence_of_source_heterogeneity", "heterogeneity",
                         fifelse(evidence_class == "evidence_of_practical_equivalence", "equivalence",
                         fifelse(evidence_class == "asymmetric_support_without_demonstrated_heterogeneity", "asymmetry", "naive")))]
p1b <- ggplot(decision, aes(x = n, y = display, fill = color_key)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = n), hjust = -0.18, fontface = "bold", size = 3.4, color = col_dark) +
  scale_fill_manual(values = c(naive = "#526875", heterogeneity = col_hetero, equivalence = col_equiv, asymmetry = col_asym)) +
  scale_x_continuous(limits = c(0, 300), expand = expansion(mult = c(0, 0))) +
  labs(title = "Claim adjudication", subtitle = "All instruments; rho = 0; n = 1,060", x = "Comparisons", y = NULL) +
  theme_method(8.0) + theme(legend.position = "none", panel.grid.major.y = element_blank())

# Figure 1C: common-SNP adjudication plus stability metrics.
common_bar <- common[metric %in% c("Common-SNP naive rows", "FDR heterogeneity", "Practical equivalence", "Asymmetric support")]
common_bar[, display := factor(metric, levels = rev(c("Common-SNP naive rows", "FDR heterogeneity", "Practical equivalence", "Asymmetric support")))]
common_bar[, color_key := c("naive", "heterogeneity", "equivalence", "asymmetry")]
p1c <- ggplot(common_bar, aes(x = value, y = display, fill = color_key)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = format(value, trim = TRUE)), hjust = -0.18, fontface = "bold", size = 3.4, color = col_dark) +
  scale_fill_manual(values = c(naive = "#526875", heterogeneity = col_hetero, equivalence = col_equiv, asymmetry = col_asym)) +
  scale_x_continuous(limits = c(0, 265), expand = expansion(mult = c(0, 0))) +
  labs(title = "Common-instrument analysis", subtitle = "986/1,060 classes agree; r = 0.9344; 33/37 calls retained", x = "Comparisons", y = NULL) +
  theme_method(8.0) + theme(legend.position = "none", panel.grid.major.y = element_blank())

# Figure 1D: subgroup context.
impact_sub <- impact[group_type %in% c("exposure_family", "field_domain")]
impact_long <- melt(
  impact_sub,
  id.vars = c("group_type", "group_id", "n_pairs"),
  measure.vars = c("naive_source_restricted_heterogeneity_rate_pct", "naive_source_restricted_equivalence_rescue_pct", "naive_source_restricted_downgrade_pct"),
  variable.name = "metric", value.name = "pct"
)
impact_long[, metric := factor(metric,
  levels = c("naive_source_restricted_heterogeneity_rate_pct", "naive_source_restricted_equivalence_rescue_pct", "naive_source_restricted_downgrade_pct"),
  labels = c("Retained heterogeneity", "Practical equivalence", "Asymmetric support")
)]
impact_long[, group_label := gsub("_", " ", group_id)]
impact_long[, group_label := tools::toTitleCase(group_label)]
impact_long[, group_label := fifelse(group_label == "Immune Inflammatory Musculoskeletal", "Immune/inflammatory/MSK",
                              fifelse(group_label == "Digestive Respiratory", "Digestive/respiratory",
                              fifelse(group_label == "Glycemic Diabetes", "Glycemic/diabetes", group_label)))]
impact_long[, group_type := factor(group_type, levels = c("exposure_family", "field_domain"), labels = c("Exposure\nfamilies", "Endpoint\ndomains"))]
p1d <- ggplot(impact_long, aes(x = pct, y = reorder(group_label, pct), color = metric)) +
  geom_point(size = 2.15, alpha = 0.95) +
  facet_grid(group_type ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = c("Retained heterogeneity" = col_hetero, "Practical equivalence" = col_equiv, "Asymmetric support" = col_asym)) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), labels = label_percent(scale = 1)) +
  labs(title = "Stratified decision classes", x = "Share of naive-asymmetry rows", y = NULL) +
  theme_method(7.2) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 6.2),
    strip.text.y.right = element_text(angle = 0, size = 6.0, lineheight = 0.9)
  )

fig1 <- ((p1a + p1b) / (p1c + p1d)) + plot_annotation(tag_levels = "A") & tag_theme

# Figure 2A: rho sensitivity.
p2a <- ggplot(rho_summary, aes(x = rho, y = heterogeneity_n)) +
  geom_hline(yintercept = 243, linetype = "dashed", color = col_asym, linewidth = 0.8) +
  geom_line(color = col_hetero, linewidth = 1.1) +
  geom_point(color = col_hetero, size = 2) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "#526875") +
  annotate("text", x = -0.28, y = 236, label = "243 naive asymmetry rows", hjust = 0, vjust = 1, color = col_asym, size = 3.1) +
  annotate("label", x = 0, y = 37, label = "rho = 0\n37", vjust = -0.45, size = 2.8, fill = "white", linewidth = 0.2, color = col_dark) +
  scale_x_continuous(breaks = seq(-0.3, 0.5, 0.1), labels = label_number(accuracy = 0.1)) +
  scale_y_continuous(limits = c(0, 260), breaks = c(0, 50, 100, 150, 200, 250)) +
  labs(title = "Correlation sensitivity", subtitle = "FDR heterogeneity: 32-70 of 1,060", x = "Assumed correlation (rho)", y = "Comparisons") +
  theme_method(7.5) + theme(legend.position = "none")

# Figure 2B: independent evidence-axis trajectories.
axis_colors <- c("FDR-supported source difference" = col_hetero, "Within declared 20% margin" = col_equiv, "One-source nominal support" = col_asym)
rho_axes[, metric := factor(metric, levels = names(axis_colors))]
p2b <- ggplot(rho_axes, aes(x = rho, y = n, color = metric)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.3) +
  scale_color_manual(values = axis_colors) +
  scale_x_continuous(breaks = seq(-0.3, 0.5, 0.1), labels = label_number(accuracy = 0.1)) +
  labs(title = "Evidence axes across rho", x = "Assumed correlation (rho)", y = "Comparisons") +
  theme_method(7.5) + guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.text = element_text(size = 6.2))

# Figure 2C: calibrated simulation.
sim[, common_or := factor(sprintf("%.2f", exp(common_effect)), levels = sprintf("%.2f", sort(unique(exp(common_effect)))))]
sim[, diff_or := factor(sprintf("%.2f", exp(true_source_difference)), levels = sprintf("%.2f", sort(unique(exp(true_source_difference)))))]
sim_long <- melt(sim, id.vars = c("scenario_id", "common_or", "diff_or"),
                 measure.vars = c("old_asymmetric_support_rate", "fdr_heterogeneity_rate"),
                 variable.name = "rule", value.name = "rate")
sim_long[, rule := factor(rule, levels = c("old_asymmetric_support_rate", "fdr_heterogeneity_rate"),
                          labels = c("Naive asymmetry", "FDR heterogeneity"))]
p2c <- ggplot(sim_long, aes(x = diff_or, y = common_or, fill = rate * 100)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.1f", rate * 100)), size = 2.45, color = "#17252F") +
  facet_wrap(~rule, ncol = 2) +
  scale_fill_gradientn(colors = c("#F3F6F8", "#B7CEE0", "#6B5CA5"), limits = c(0, 100), oob = squish) +
  labs(title = "Calibration simulation", subtitle = "Operating rates (%)", x = "True source-difference OR ratio", y = "Common-effect OR", fill = "%") +
  theme_method(7.1) + theme(legend.position = "right", legend.title = element_text(size = 6.5), legend.text = element_text(size = 6.2))

# Figure 2D: prespecified direct-association transfer check.
nonmr_matrix[, statistical_source_difference := factor(statistical_source_difference, levels = c("unsupported", "supported"))]
nonmr_matrix[, practical_compatibility := factor(practical_compatibility, levels = c("not_demonstrated", "within_declared_margin"))]
nonmr_matrix[, nominal_asymmetry := nominal_support_pattern == "one"]
nonmr_summary[, statistical_source_difference := factor(statistical_source_difference, levels = c("unsupported", "supported"))]
nonmr_summary[, practical_compatibility := factor(practical_compatibility, levels = c("not_demonstrated", "within_declared_margin"))]
cell_bg <- data.table(
  statistical_source_difference = factor(c("unsupported", "unsupported", "supported", "supported"), levels = c("unsupported", "supported")),
  practical_compatibility = factor(c("not_demonstrated", "within_declared_margin", "not_demonstrated", "within_declared_margin"), levels = c("not_demonstrated", "within_declared_margin")),
  fill = c("#F1F3F5", "#E4F2EB", "#EEE8F6", "#E7ECEB")
)
p2d <- ggplot() +
  geom_tile(data = cell_bg, aes(x = statistical_source_difference, y = practical_compatibility, fill = fill), width = 0.98, height = 0.98, color = "white", linewidth = 1) +
  scale_fill_identity() +
  geom_point(data = nonmr_matrix, aes(x = statistical_source_difference, y = practical_compatibility, shape = nominal_asymmetry),
             position = position_jitter(width = 0.22, height = 0.19, seed = 20260714), size = 2.35, stroke = 0.85, color = col_dark, fill = "white") +
  geom_label(data = nonmr_summary, aes(x = statistical_source_difference, y = practical_compatibility, label = paste0("n = ", n)),
             vjust = 1.95, size = 3.0, fontface = "bold", fill = "white", color = col_dark, linewidth = 0.18) +
  annotate("text", x = 1.5, y = 2.42, label = "two nominal-asymmetry units", size = 2.55, color = col_asym) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 24), labels = c(`FALSE` = "Other units", `TRUE` = "One-source nominal support")) +
  scale_x_discrete(labels = c(unsupported = "FDR unsupported", supported = "FDR supported")) +
  scale_y_discrete(labels = c(not_demonstrated = "Compatibility\nnot demonstrated", within_declared_margin = "Within declared\n10% margin")) +
  labs(title = "Direct-association transfer", subtitle = "25 paired contrasts", x = "Statistical source difference", y = "Practical compatibility") +
  theme_method(7.1) + theme(legend.text = element_text(size = 6.2), legend.position = "bottom", panel.grid = element_blank())

# Supplementary Figure S2: previous third-source stress test retained for completeness.
third_total <- data.table(
  metric = factor(c("Planned", "Complete", "Nominal asymmetry", "FDR heterogeneity"), levels = rev(c("Planned", "Complete", "Nominal asymmetry", "FDR heterogeneity"))),
  n = c(sum(third$planned_rows), sum(third$complete_rows), sum(third$asymmetry_rows), sum(third$heterogeneity_rows)),
  color_key = c("planned", "complete", "asymmetry", "heterogeneity")
)
pS2 <- ggplot(third_total, aes(x = n, y = metric, fill = color_key)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = n), hjust = -0.25, fontface = "bold", size = 3.6, color = col_dark) +
  scale_fill_manual(values = c(planned = "#526875", complete = col_concord, asymmetry = col_asym, heterogeneity = col_hetero)) +
  scale_x_continuous(limits = c(0, 82), expand = expansion(mult = c(0, 0))) +
  labs(title = "Supporting third-source stress test", subtitle = "12 endpoint constructs; BioBank Japan or consortium outcomes", x = "Number of comparisons", y = NULL) +
  theme_method(9.2) + theme(legend.position = "none", panel.grid.major.y = element_blank())

fig2 <- ((p2a + p2b) / (p2c + p2d)) + plot_annotation(tag_levels = "A") & tag_theme

save_pair <- function(plot, stem, width = 6.69, height = 5.25) {
  png_path <- file.path(main_dir, paste0(stem, ".png"))
  pdf_path <- file.path(main_dir, paste0(stem, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = 320, bg = "white")
  ggsave(pdf_path, plot, width = width, height = height, bg = "white")
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(file.path(main_dir, paste0(stem, ".svg")), plot, width = width, height = height, bg = "white")
  }
  c(png_path, pdf_path)
}

if (mode == "pilot") {
  ggsave(file.path(pilot_dir, "Figure1_Claim_Impact_Pilot.png"), fig1, width = 13.2, height = 9.4, dpi = 200, bg = "white")
  ggsave(file.path(pilot_dir, "Figure1_Claim_Impact_Pilot.pdf"), fig1, width = 13.2, height = 9.4, bg = "white")
  ggsave(file.path(pilot_dir, "Figure2_Robustness_and_Transfer_Pilot.png"), fig2, width = 13.2, height = 9.4, dpi = 200, bg = "white")
  ggsave(file.path(pilot_dir, "Figure2_Robustness_and_Transfer_Pilot.pdf"), fig2, width = 13.2, height = 9.4, bg = "white")
  message("Pilots rendered: ", file.path(pilot_dir, "Figure1_Claim_Impact_Pilot.png"), " and ", file.path(pilot_dir, "Figure2_Robustness_and_Transfer_Pilot.png"))
} else {
  save_pair(fig1, "Figure1_Claim_Impact_Paired_Source_Reporting")
  save_pair(fig2, "Figure2_Robustness_and_Stress_Tests")
  s2 <- pS2 + plot_annotation(tag_levels = "A") & tag_theme
  s2_stem <- file.path(supp_dir, "Supplementary_Figure_S2_Third_Source_Stress_Test")
  ggsave(paste0(s2_stem, ".png"), s2, width = 6.69, height = 4.7, dpi = 320, bg = "white")
  ggsave(paste0(s2_stem, ".pdf"), s2, width = 6.69, height = 4.7, bg = "white")
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(paste0(s2_stem, ".svg"), s2, width = 6.69, height = 4.7, bg = "white")
  }
  message("Final Figures 1-2 rendered in ", main_dir)
}
