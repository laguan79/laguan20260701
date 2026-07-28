# Run Order

This release reproduces the reporting tool, derived-table analyses, sensitivity summaries, and figures from packaged inputs. It is not a raw-GWAS-to-final archive. Rebuilding source-level estimates requires the original public GWAS summary statistics and the curated inputs described in the manuscript.

Run commands from the release root with R 4.5.3. Required packages are listed in `SOFTWARE_METADATA.md`.

## 1. Mini-Report

```powershell
Rscript reporting_extension/source_effect_mini_report.R --input reporting_extension/example_input_source_effect.tsv --output reporting_extension/example_output_source_effect_report.tsv --markdown reporting_extension/example_output_source_effect_report.md --or-margin 1.20 --rho 0 --fdr-source auto
Rscript tests/check_mini_report_example.R
```

The first command rebuilds the example report. The second checks its schema, numerical fields, evidence classes, and expected row count.

## 2. Main and Supplementary Figures

```powershell
Rscript code/render_method_hardening_figures.R --final
Rscript code/render_method_hardening_figure3.R
Rscript code/render_supplementary_figure_s1.R
```

These commands regenerate Figures 1-3 and Supplementary Figures S1-S2 from the package-relative inputs listed in `figures/FIGURE_TO_SOURCE_MAPPING.tsv`.

## 3. Derived Analyses

The following scripts expose command-line help when run without arguments:

```powershell
Rscript code/source_effect_uncertainty_classifier.R
Rscript code/run_margin_rho_sensitivity.R
Rscript code/run_precision_asymmetry_simulation.R
```

`code/run_common_snp_reanalysis.R` is a provenance script for the full common-SNP reconstruction. It requires the curated analysis tree supplied separately with the archived derived-data release:

```powershell
$env:CALIBRATION_INPUT_ROOT="<path to restored curated calibration inputs>"
Rscript code/run_common_snp_reanalysis.R
```

## 4. FinnGen R13 Supplementary Stress Test

The R13 code is namespaced under `code/r13_release_drift/`. Its scripts consume the `inputs/`, `outputs/`, `qa/`, and `metadata/` paths supplied by Supplementary Table S14 in the derived-data archive. Stage those directories at the release root before running the sequence documented in `code/r13_release_drift/README.md`.

R13 is a supplementary release-drift stress test, not an independent validation module.

## 5. Supplementary Inventory

The derived-data archive covers Supplementary Tables S1-S14, figure source data, the reporting-extension example, software metadata, session information, and checksums. `ADDITIONAL_FILE_MANIFEST.md` maps the internal S1-S14 identifiers to journal upload filenames.
