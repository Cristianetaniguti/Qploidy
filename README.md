<!-- badges: start -->
[![Development](https://img.shields.io/badge/development-active-blue.svg)](https://img.shields.io/badge/development-active-blue.svg)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![codecov](https://codecov.io/github/cristianetaniguti/qploidy/branch/main/graphs/badge.svg)](https://codecov.io/github/cristianetaniguti/qploidy)

<!-- badges: end -->

# Qploidy <img src="https://github.com/Cristianetaniguti/Qploidy/assets/7572527/88ef9fad-7f86-4a84-9e1a-5dd4625dd1c8" align="right" width="230"/>

`Qploidy` is and R package and Shiny app to perform ploidy and aneuploid estimation using genotyping platforms data such as intensities (normalized X and Y) from Axiom and Illumina array genotyping platforms and read counts for each allele from target sequencing platforms. 

## Installation

``` r
#install.packages("devtools")
devtools::install_github("cristianetaniguti/Qploidy")
```

## How to use

You can use Qploidy using its graphical interface by:

``` r
library(Qploidy)
run_app()
```

Or you can also use its functions individually.

* [Interpolation and ploidy estimation - code version]()


