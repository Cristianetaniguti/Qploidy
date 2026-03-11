# Qploidy 1.8.1

* Added cross-parameter validation in `standardize()`: `threshold.n.clusters` is now checked to not exceed `ploidy.standardization + 1`, with an informative error message.
* `select_best_baf_model()` now supports dual-input mode: BAF values can be supplied directly via `baf_vec` (with an optional `chr_vec`) or extracted from a `qploidy_standardization` object and a sample name. Both paths are mutually exclusive and validated.
* Added diagnostic warning in `hmm_estimate_CN()` when all heterozygous markers in one or more windows are discarded by the `dosage_threshold` filter, reporting the affected windows by chromosome and window ID.
* Added `parallel_type` argument to `hmm_estimate_CN_multi()` (default `"auto"`): automatically selects `"FORK"` on Unix/macOS (faster, no symbol re-export needed) and `"PSOCK"` on Windows. Can be set explicitly to `"FORK"` or `"PSOCK"`.
* Warnings emitted inside parallel workers are now captured and re-issued on the main R session after `parLapply`, preventing diagnostic messages from being silently lost in PSOCK clusters.
* `select_best_baf_model()` now emits a warning when the top three BIC-ranked models disagree on the best CN estimate, alerting users to ambiguous model selection.
* `hmm_estimate_CN()` now excludes CN = 1 as a candidate state for windows that contain BAF values within `het_range` (default `c(0.2, 0.8)`), preventing spurious CN = 1 calls driven by a low heterozygous-to-homozygous ratio.

# Qploidy 1.8.0

* Enhanced console messages for standardization and HMM steps, improving clarity and user feedback.
* Improved CN grid selection: unlike CN values are now removed after the EM loop by checking z-score means. The CN grid is reordered from lower to higher ploidy, and EM is rerun if disruptive CN values are detected.
* The min_snp_per_window argument now defaults to 10% of the smallest chromosome size, with a minimum of 5 SNPs, for more adaptive windowing.
* Added an argument to define z-score intervals by subtracting z_range (instead of adding), resulting in smaller intervals.
* Refactored EM loop into a dedicated function for better modularity and maintainability.
* Updated dosage call plot colors for improved visual distinction and interpretability.
* Exception created in HMM when expected ploidy = 1. In this case, heterozygous are not expected so the BAF weight is set to 0.5 regardless if heterozygous are present or not. Parameter can still be controled by baf_weight argument (default 1).

# Qploidy 1.7.1

* New function depth_pca_plot to check for batch effects using a pca based on total depth
* Allow input data for function chrom_ttest, window_ttest and depth_pca_plot to have only MarkerName, SampleName, and R columns if geno.pos data.frame is provided with Chromosome and Position information.
* Bugfix on DESCRIPTION authors 
* Improve documentation for new functions

# Qploidy 1.7.0

* New functions chrom_ttest and window_ttest
* Adding to plot_geno_by_marker:
  - transparency argument alpha
  - identify decreasing ratios by dosage for expected dot plot

# Qploidy 1.6.9

* Add filter_R argument on standardize. If FALSE (default) filters defined are not applied for z score calculation
* Add depth_zero_as_x to plot_cn_track. Total counts (R) of missing genotypes are counted as 0 and converted to z score. depth_zero_as_x = TRUE mark on the graphic which z values refer to the 0 total counts

# Qploidy 1.6.8

* Increase speed to export VCF
* Add re_standardize function to run standardization based on HMM results
* Reduce text on plot_cn_track

# Qploidy 1.6.7

* Make initial probability 0.95 for the expected ploidy
* Bugfix plot_cn_track

# Qploidy 1.6.6

* Add correction factor for Z and BAF likelihood to have same weight if BAF_weight = 0.5 (new default)
* Filter low probability heterozygous for counts on determining the BAF_weight
* If correct_scale = TRUE (default) BAF likelihood is corrected by the number of markers (with BAF values) used for the distribution
* plot_cn_track now contains the sample-level BAF histogram and parameters descriptions written on the figure

# Qploidy 1.6.5

* Number of heterozygous to define the BAF weight now is counted using dosage call based only on BAF
* Dosage call functions added based only on BAF and based on BAF according to HMM CN call result
* Export called dosages in a VCF format
* Add co-pilot-instructions to the repository

# Qploidy 1.6.2

* Avoid markers without chromosome information in the standardized dataset
* Add plot to compare ratios and standardazed ratios (BAF) in a geno vs value format
* Add new type of plot for plot_qploidy_standardized for raw total depth (R)
* Add checks for standardize input files. It now requires specific column names and ordering.


# Qploidy 1.6.1

* write and read function for hmm_CN object
* print function summarizing estimated paramters in hmm_CN object
* hmm now return two data.frames, one with information by window and another by marker
* hmm returns list of parameters defined by user and estimated by function for reproducibility

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



