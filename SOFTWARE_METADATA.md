# Software Metadata

- Project name: Paired-Source Statistical Reporting for Cross-Biobank Genomic Evidence
- Release version: 1.0.0
- Article type: Method
- Code repository: https://github.com/laguan79/laguan20260701
- Derived-data archive: pending permanent archive registration before journal submission
- Operating system: tested on Windows 11; package-relative R scripts are platform independent.
- Programming language: R 4.5.3
- Core requirements: data.table, ggplot2, patchwork, grid, scales, cowplot, svglite, and TwoSampleMR; exact versions are recorded in the release session information.
- License: MIT for release code and documentation.
- Use restrictions: derived data are provided for reproducibility of this manuscript; reuse of original public GWAS summary statistics remains governed by the source repositories.
- Directly reproducible components: mini-report example, classifier and sensitivity analyses using packaged inputs, Figures 1-3, Supplementary Figures S1-S2, selected published-case numerical subset, derived tables, and expected-output tests.
- Supplementary addon: `code/r13_release_drift/` reproduces the FinnGen R13 release-drift stress test after the archived S14 inputs are staged at the release root.
- Not directly reproduced from this release: full raw-GWAS-to-final rebuilding, which requires original public GWAS summary statistics and curated source-specific inputs.
