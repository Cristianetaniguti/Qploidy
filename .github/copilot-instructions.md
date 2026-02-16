# Copilot Instructions for Qploidy

## Project Overview
Qploidy is an R package for ploidy and aneuploidy estimation using genotyping array and targeted sequencing data. It supports Axiom, Illumina, and similar platforms, but not RADseq/GBS. The package includes both R functions and a Shiny app interface.

## Key Components
- **R/**: Core R source files. Main logic for ploidy estimation, HMM, BAF calculations, and plotting.
- **inst/app/**: Shiny app UI and server logic for interactive analysis.
- **man/**: R documentation files for all exported functions.
- **example_data_run/**: Example datasets and scripts for testing and demonstration.
- **tests/**: Testthat-based tests for package validation.

## Developer Workflows
- **Install dependencies**: Use standard R package installation (`install.packages()` or `devtools::install_github()`).
- **Run the Shiny app**: `library(Qploidy); run_app()`
- **Testing**: Run `devtools::test()` or use the testthat workflow in RStudio.
- **Build/check**: Use `devtools::check()` for R CMD check compliance.
- **Documentation**: Update `.Rd` files in `man/` using roxygen2 comments in `R/` and run `devtools::document()`.

## Project Conventions
- All exported functions are documented in `man/` and use roxygen2-style comments.
- Data input/output is typically via CSV or RData files; see `example_data_run/` for formats.
- Shiny app configuration is in `R/app_config.R`.
- HMM and BAF logic is in `R/hmm_main.R`, `R/hmm_utils.R`, and `R/BAF_distributions.R`.
- Use `run_app()` as the main entry point for interactive analysis.

## Integration & Dependencies
- Relies on standard R packages (see `DESCRIPTION` for dependencies).
- Shiny app uses the `shiny` package and related UI/server modules.
- No external services or APIs are required; all computation is local.

## Examples
- To estimate ploidy from a VCF: see `R/qploidy_read_vcf.R` and `R/ploidy_est.R`.
- For plotting: see `R/plot_baf.R`, `R/plot_cn_track.R`, and related files.
- Example workflow: `example_data_run/plot_dev.R`.

## Special Notes
- Do not combine datasets from different sequencing batches in a single analysis (see README for rationale).
- Always ensure a subset of samples with known or predominant ploidy is present for accurate estimation.

---
For more details, see the [README.md](../README.md) and function documentation in `man/`.
