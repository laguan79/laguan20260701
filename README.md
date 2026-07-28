# Paired-Source Statistical Reporting for Cross-Biobank Genomic Evidence

This reproducibility package implements the paired source-effect reporting procedure used in the manuscript. The direct comparison is the standard one-degree-of-freedom Wald contrast (equivalent to a two-study Cochran Q under independence); the contribution is the executable reporting process around that contrast.

## Package status

The expanded code tree and `v1.0.0` tag are public. The scientific content, Supplementary Tables S1-S14, figure scripts, and reproducibility checks are frozen for archive deposit. Journal submission remains blocked until the permanent derived-data archive has been registered and its URL and DOI have been inserted into the submission files.

## Evidence outputs

Each estimable comparison returns three non-exclusive axes: FDR-supported versus unsupported statistical source difference, practical compatibility within versus outside a declared margin, and nominal support in both, one, or neither source. Legacy priority-ordered labels are retained only for audit mapping.

## Directly reproducible checks

```powershell
Rscript reporting_extension/source_effect_mini_report.R --input reporting_extension/example_input_source_effect.tsv --output reporting_extension/example_output_source_effect_report.tsv --markdown reporting_extension/example_output_source_effect_report.md --or-margin 1.20 --rho 0 --fdr-source auto
Rscript tests/check_mini_report_example.R
```

The same release also rebuilds Figures 1-3 and Supplementary Figures S1-S2:

```powershell
Rscript code/render_method_hardening_figures.R --final
Rscript code/render_method_hardening_figure3.R
Rscript code/render_supplementary_figure_s1.R
```

`RUN_ORDER.md` documents the analysis scripts, expected inputs, and R13 supplementary addon. The expanded code tree and `v1.0.0` tag are available at https://github.com/laguan79/laguan20260701. The permanent archive will include Supplementary Tables S1-S14, derived data, figure source data, code, tests, expected outputs, session metadata, and SHA-256 manifests. Supplementary Table S14 is a frozen release-drift stress test.

## Archive identifiers

- Code repository: https://github.com/laguan79/laguan20260701
- Derived-data archive: pending permanent archive registration before journal submission
