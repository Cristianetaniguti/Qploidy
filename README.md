<!-- badges: start -->
[![Development](https://img.shields.io/badge/development-active-blue.svg)](https://img.shields.io/badge/development-active-blue.svg)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/Cristianetaniguti/Qploidy/workflows/R-CMD-check/badge.svg)](https://github.com/Cristianetaniguti/Qploidy/actions)
![](https://img.shields.io/badge/RRID-SCR_026724-yellow.svg)
[![codecov](https://codecov.io/gh/Cristianetaniguti/Qploidy/graph/badge.svg?token=DQBM227JSY)](https://codecov.io/gh/Cristianetaniguti/Qploidy)
[![CRAN version](https://www.r-pkg.org/badges/version/Qploidy)](https://CRAN.R-project.org/package=Qploidy)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/Qploidy)](https://CRAN.R-project.org/package=Qploidy)

<!-- badges: end -->

# Qploidy 

<img src="https://github.com/Cristianetaniguti/Qploidy/assets/7572527/88ef9fad-7f86-4a84-9e1a-5dd4625dd1c8" align="right" width="230"/>


**`Qploidy`** is an R package designed for ploidy and aneuploidy estimation using genotyping platform data. 

### When does the `Qploidy` methodology work?

The `Qploidy` approach is effective under the following conditions:  
- Your marker data originates from **Axiom** or **Illumina genotyping arrays**.  
- Your marker data is derived from **targeted sequencing platforms** (e.g., DArTag, GTseq, AgriSeq).  
- All DNA samples were prepared following the **same library preparation protocol**.  
- You know the **ploidy** of at least a subset of 60 samples *or* you know the **most common ploidy** in the dataset.  
- Your dataset includes **heterozygous samples**.  

### When does the `Qploidy` methodology NOT work?

The methodology will not be effective under the following circumstances:  
- Your marker data comes from **RADseq** or **GBS** (Genotyping-by-Sequencing) platforms.  
- You intend to **combine datasets from different sequencing batches**.  
   - For example: If you extracted DNA and sequenced two plates (192 samples) as one batch, and later sequenced an additional three plates (288 samples) as a second batch, you would need to analyze the two batches **separately** in `Qploidy`. Combining all 480 samples into a single analysis will lead to incorrect results.  
- You **do not have a subset of samples with known ploidy** or **lack a predominant ploidy** in your dataset.  
- Your samples consist of **inbred lines** (homozygous individuals).  


## Installation

``` r
# Install the latest stable version from CRAN
install.packages("Qploidy")

# Install the development version from GitHub
#install.packages("devtools")
devtools::install_github("cristianetaniguti/Qploidy")
```

## Access Shiny interface

``` r
library(Qploidy)
run_app()
```

## Documentation

* [`Qploidy` tutorial](https://cristianetaniguti.github.io/Qploidy/Qploidy.html) for directions on how to run

## Contributing

Contributions are welcome! If you'd like to contribute, please fork the repository and submit a pull request. For major changes, open an issue first to discuss your ideas.

## Bug Reports

If you find a bug or want an enhancement, please submit an issue [here](https://github.com/Cristianetaniguti/Qploidy/issues).

## How to cite

Taniguti, C. H., Lau, J., Hochhaus, T., Arias, D. C. L., Hokanson, S. C., Zlesak, D. C., Byrne, D. H., Klein, P. E., & Riera-Lizarazu, O. (2025). Exploring chromosomal variations in garden roses: Insights from high-density SNP array data and a new tool, Qploidy. The Plant Genome, e70044. https://doi.org/10.1002/tpg2.70044

## Acknowledgments

### Version 1.0.0

This work is funded in part by the Robert E. Basye Endowment in Rose Genetics, Dept. of Horticultural Sciences, Texas A&M University, and USDA’s National Institute of Food and Agriculture (NIFA), Specialty Crop Research Initiative (SCRI) projects: ‘‘Tools for Genomics-Assisted Breeding in Polyploids: Development of a Community Resource’’ (Award No. 2020-51181-32156); and ‘‘Developing Sustainable Rose Landscapes via Rose Rosette Disease Education, Socioeconomic Assessments, and Breeding RRD-Resistant Roses with Stable Black Spot Resistance’’ (Award No. 2022-51181-38330).

### Versions > 1.0.0

Supported by [Breeding Insight](https://www.breedinginsight.org/).
