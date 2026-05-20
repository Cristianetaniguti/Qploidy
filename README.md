<!-- badges: start -->
[![Development](https://img.shields.io/badge/development-active-blue.svg)](https://img.shields.io/badge/development-active-blue.svg)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/Breeding-Insight/Qploidy2/workflows/R-CMD-check/badge.svg)](https://github.com/Breeding-Insight/Qploidy2/actions)
[![codecov](https://codecov.io/gh/Breeding-Insight/Qploidy2/graph/badge.svg?token=DQBM227JSY)](https://codecov.io/gh/Breeding-Insight/Qploidy2)

<!-- badges: end -->

# Qploidy2

<img src="inst/app/www/Qploidy2_logo.png" align="right" width="230"/>


**`Qploidy2`** is an R package designed for large-scale copy number variation, ploidy and aneuploidy estimation using genotyping platforms data. 

`Qploidy2` builds upon the original [`Qploidy`](https://github.com/Cristianetaniguti/Qploidy) package. While `Qploidy` focused on signal standardization to reduce noise for ploidy estimation and relied heavily on visual inspection of diagnostic plots, `Qploidy2` introduces a multipoint approach using a Hidden Markov Model (HMM). This advanced method combines information from standardized allele ratios and total depth/intensities to provide automated and robust copy number estimation.

The package offers flexible implementation: you can access Qploidy's original standardization methods within `Qploidy2`, or apply the HMM directly to your data without standardization if noise reduction isn't needed.

## Installation

``` r
# Install the development version from GitHub
#install.packages("devtools")
devtools::install_github("Breeding-Insight/Qploidy2")
```

## Documentation

* [Tutorial](https://Breeding-Insight.github.io/Qploidy2/Qploidy_alfalfa_tutorial.html) - large-scale copy number estimation in alfalfa mapping population using DArTag sequencing data

## Contributing

Contributions are welcome! If you'd like to contribute, please fork the repository and submit a pull request to the development branch. For major changes, open an issue first to discuss your ideas.

## Bug Reports

If you find a bug or want an enhancement, please submit an issue [here](https://github.com/Breeding-Insight/Qploidy2/issues).

## How to cite

### Qploidy - Standardization 

Taniguti, C. H., Lau, J., Hochhaus, T., Arias, D. C. L., Hokanson, S. C., Zlesak, D. C., Byrne, D. H.,
Klein, P. E., & Riera-Lizarazu, O. (2025). Exploring chromosomal variations in garden roses: Insights 
from high-density SNP array data and a new tool, Qploidy. The Plant Genome, e70044. 
https://doi.org/10.1002/tpg2.70044

### Qploidy2 - HMM and grid approach BAF model selection

Manuscript in preparation. Please contact the author for more information.

### nQuack - EM approach BAF model selection

Gaynor, M., Landis, J., O'Connor, T., Laport, R., Doyle, J., Soltis, D., Ponciano, J., & Soltis, P. (2024). 
"nQuack: An R package for predicting ploidal level from sequence data using site-based heterozygosity." 
Applications in Plant Sciences, 12(4), e11606. doi:10.1002/aps3.11606

## Acknowledgments

### Qploidy - version 1.0.0

This work is funded in part by the Robert E. Basye Endowment in Rose Genetics, Dept. of Horticultural Sciences, Texas A&M University, and USDA’s National Institute of Food and Agriculture (NIFA), Specialty Crop Research Initiative (SCRI) projects: ‘‘Tools for Genomics-Assisted Breeding in Polyploids: Development of a Community Resource’’ (Award No. 2020-51181-32156); and ‘‘Developing Sustainable Rose Landscapes via Rose Rosette Disease Education, Socioeconomic Assessments, and Breeding RRD-Resistant Roses with Stable Black Spot Resistance’’ (Award No. 2022-51181-38330).

### Qploidy versions > 1.0.0 and Qploidy2

Supported by [Breeding Insight](https://www.breedinginsight.org/).
