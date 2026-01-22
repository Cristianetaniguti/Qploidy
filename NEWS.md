# Qploidy 1.6.0 

* Major refactor and improvements to BAF likelihood/model selection workflow and plotting 
* Added support for testing different emission distributions for BAF (nQuack idea), including negative binomial, with automated model selection via BIC
* Refactoring of BAF template generation and likelihood computation into dedicated functions
* Improved plotting functions:
    - plot_cn_track has now a fallback to hmm_CN$updated_data if qploidy_standarize_result is NULL
    - function compare_cn_track added
    - Improved color palette for CN plots
    - Chromosome sorting and x-axis labeling improved for clarity
    - Legends and color mapping clarified and improved
    - Subsetting chromosomes on plot - not only on hmm
* HMM segmentation and CN calling:
    - Implementation of changepoint-based z-score segmentation as an alternative to fixed-window approaches
    - Single HMM across all chromosomes (not per-chromosome)
    - Robust outlier handling for z-scores (rm_outlier supports z argument and outlier column)
    - Exception handling for single-window case (assigns CN by BAF likelihood only)
    - exp_ploidy argument added and integrated throughout (including Shiny UI)
    - Improved verbose output and summary messages
    - By default, z_range is now automatically calculated based on the inter-quantile range of z-scores and the number of CN states tested, for more robust and adaptive HMM initialization
* User experience and error handling:
    - Improved error messages and fallback logic in plotting and HMM functions
    - Documentation expanded and clarified for all major functions
* Shiny UI:
    - Sample-level BAF distributions plot added 
    - exp_ploidy is now an advanced option for z_only=TRUE
    - Argument passing and advanced options printing improved

# Qploidy 1.5.2

* Z-score calculation is now only applied to markers not filtered by genotype probabilities and missing data
* Markers filtered out during standardization are not considered to define the window boundaries for HMM CNV estimation 
* Last window of each chromosome are never smaller than the defined window size
* CN HMM plots follow same scale for easier comparison between metrics by window

# Qploidy 1.5.1

* Fix HMM CNV estimation when some of the windows don't have heterozygous genotypes (use only z-score for them)
* Improve Shiny interface for HMM CNV estimation

# Qploidy 1.5.0

* Add HMM to estimate CN (beta version)
* Shiny interface improvements
* Add CSV or TSV as possible input files

# Qploidy 1.1.0

* Beta version of Qploidy + BIGapp interface
* Fix Bioconductor dependency call
* Fix %>% call

# Qploidy 1.0.0

* First release
* Functions improvements including comprehensive documentation
* Vignette with functional simulated example
* Testthat tests included
* CRAN submission

# Qploidy 0.0.0.9000

* This is a beta version
* Initial addition of functions
* First version of the vignette



