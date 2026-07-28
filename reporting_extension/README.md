# Source-Effect Reporting Mini-Report

This folder contains a lightweight reporting script for public-source comparisons in the MR calibration setting. It is not an R package. It is a script-level companion to the proposed minimum source-effect reporting items.

## Files

- `source_effect_mini_report.R`: classifies paired source-specific estimates and writes a report table.
- `example_input_source_effect.tsv`: three worked examples selected from the common-SNP calibration reanalysis.
- `example_output_source_effect_report.tsv`: mini-report output for the worked examples.
- `example_output_source_effect_report.md`: reader-facing mini-report summary.
- `example_output_source_effect_report_class_summary.tsv`: class-count summary for the worked examples.

## Example command

```powershell
Rscript .\source_effect_mini_report.R `
  --input .\example_input_source_effect.tsv `
  --output .\example_output_source_effect_report.tsv `
  --markdown .\example_output_source_effect_report.md `
  --or-margin 1.20 `
  --rho 0 `
  --fdr-source auto
```

The example input includes source-difference FDR values from the full 1,060-row common-SNP calibration family. For a new study, omit the FDR column to compute Benjamini-Hochberg FDR across the supplied input rows, or supply a family-level FDR column from the full comparison set.
