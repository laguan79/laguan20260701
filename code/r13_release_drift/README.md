# FinnGen R13 Release-Drift Analysis Code

This directory contains calculation and quality-control code for the
within-source temporal release-drift stress test.

## Scope

The scripts operate on a frozen 20-exposure by 53-endpoint derived-data
benchmark. They calculate:

- FinnGen R12, FinnGen R13, and Pan-UKB fixed-instrument estimates;
- paired source and release contrasts;
- BH/BY and covariance-scenario sensitivity;
- practical-equivalence and evidence-class assignments;
- crossed exposure-by-endpoint bootstrap intervals;
- full leave-one-endpoint and leave-one-exposure decision reruns;
- targeted exact-Tabix standard-error sensitivity.

## Required Derived Inputs

The archived derived-data release supplies the expected `inputs/`, `outputs/`,
and `qa/` objects. A raw-GWAS-to-final rebuild requires the original public
GWAS summary statistics and is outside this code-only directory.

## Software

- R 4.5.3
- `data.table`
- `TwoSampleMR`
- `digest`
- `httr2`
- `Rsamtools`
- `GenomicRanges`

## Primary Commands

```bash
Rscript code/r13_release_drift/analyse_release_stability.R
Rscript code/r13_release_drift/summarize_multiplicity_sensitivity.R
Rscript code/r13_release_drift/run_full_leave_one_decision_pipeline.R
```

The targeted exact-SE sensitivity uses the exact-query lock and release
indexes in the derived-data archive:

```bash
Rscript code/r13_release_drift/extract_targeted_exact_tabix.R
Rscript code/r13_release_drift/summarize_targeted_exact_precision.R
Rscript code/r13_release_drift/estimate_targeted_exact_se_sensitivity.R
```

The scripts use package-relative paths and write results under `outputs/` and
`qa/`.

## License

MIT License. Derived data remain subject to the terms of the original public
GWAS sources.
