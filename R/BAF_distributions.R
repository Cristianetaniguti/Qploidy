# This script contains functions for B-Allele Frequency (BAF) modeling and likelihood estimation:
#
# - compute_baf_likelihoods: Computes BAF likelihoods for a vector of BAF values across multiple copy number (CN) states and distributions. Optionally generates a plot comparing observed and template densities.
# - select_best_baf_model: Performs grid search over distributions, bandwidths, and uniform noise parameters to select the best BAF model by BIC. Returns best model info and (optionally) the best plot.
# - generate_baf_template: Builds a polyploid BAF comb template for a given copy number and distribution, with optional uniform noise component.
# - baf_log_likelihood: Computes the multinomial-style log-likelihood of observed BAF histogram counts against a template.

# ------------------------------------------------------------------------------
# Context
# -------
# The aim to choose the best copy-number / ploidy model by comparing
# how well candidate models explain allele-balance information.
#
# IMPORTANT: The idea of testing multiple candidate models (e.g., normal/beta/
# beta-binomial) and selecting the best one is inspired by the nQuack paper.
#
# Citation (nQuack)
# -----------------
# Gaynor, M., Landis, J., O'Connor, T., Laport, R., Doyle, J., Soltis, D.,
# Ponciano, J., & Soltis, P. (2024).
# "nQuack: An R package for predicting ploidal level from sequence data using
# site-based heterozygosity."
# Applications in Plant Sciences, 12(4), e11606. doi:10.1002/aps3.11606
#
# Notes:
# - The paper describes model selection and expectation-maximization (EM)
#   algorithms for normal, beta, and beta-binomial distributions.
#
# ------------------------------------------------------------------------------
# 1) What Qploidy current approach does
# ------------------------------------------------------------------------------
#
# A) Fixed "templates" rather than fitted mixture parameters
#    - We pre-construct expected BAF templates for each CN state:
#        * peaks at genotype fractions d/c for d = 0..c
#        * peak heights often from Binomial(c, 0.5) weights
#        * peak width controlled by bw (provided by user or chosen by grid)
#        * optional uniform component (if add_uniform = TRUE)
#    - These are "shape priors" / comb templates: the distribution is largely
#      determined a priori for each CN.
#
# B) Likelihood evaluated against binned BAF data (histogram)
#    - Convert baf values -> histogram counts across M bins on [0, 1].
#    - Compute log-likelihood under each CN template (baf_log_likelihood).
#    - Softmax across CN states -> relative probabilities.
#
# C) Model selection = best matching template configuration
#    - CN chosen by max log-likelihood (or max softmax probability).
#    - If expanded to multiple distributions and/or grid over bw/uniform_weight,
#      select best model using BIC/AIC.
#
# ------------------------------------------------------------------------------
# 2) What nQuack does (high-level)
# ------------------------------------------------------------------------------
#
# A) Explicit mixture modeling + EM estimation
#    - nQuack fits mixture models using EM for candidate distributions
#      (normal, beta, beta-binomial).
#    - EM iteratively estimates model parameters (mixture weights and/or
#      dispersion/shape parameters) to maximize likelihood.
#
# B) Uses count/depth information (site-based heterozygosity)
#    - nQuack is designed around site-based sequencing information (allele
#      counts + total depth), which is naturally modeled by (beta-)binomial
#      families and supports depth-aware dispersion.
#
# C) Model selection integrated
#    - nQuack compares candidate models and selects the best supported one,
#      improving ploidy predictions (as described in the paper).
#
# ------------------------------------------------------------------------------
# 3) Key differences (summary bullets)
# ------------------------------------------------------------------------------
#
# (1) Parameter fitting
#     - Qploidy: no EM fitting; templates are predefined (or chosen from a small grid).
#     - nQuack: EM is used to estimate mixture parameters for each candidate model.
#
# (2) Data representation
#     - Qploidy: BAF values -> histogram -> likelihood vs discrete template.
#     - nQuack: per-site allele counts + depth -> likelihood under fitted model.
#
# (3) Flexibility vs speed
#     - Qploidy: fast, deterministic per configuration; fewer parameters.
#     - nQuack: flexible but iterative; potentially slower and init-sensitive.
#
# ------------------------------------------------------------------------------
# 4) Why Qploidy DOES NOT fully reproduce the nQuack EM approach, despite admiring it
# ------------------------------------------------------------------------------
#
# We admire nQuack’s EM-based model fitting and its principled model-selection
# framework (Gaynor et al. 2024). However, our pipeline design is intentionally
# different at this stage for two main reasons:
#
# (A) Efficiency: this is an early stage in a heavier downstream model
#     - Our current step is meant to be a lightweight, efficient scoring stage.
#     - The next stage of our pipeline will already run an EM + HMM framework
#       that *jointly* integrates BAF information with z-score / coverage signals.
#     - Running full EM mixture fitting here as well would duplicate expensive
#       optimization twice (EM now + EM/HMM later), increasing runtime and
#       potentially complicating convergence/debugging without much gain.
#
#     In other words:
#       *Stage 1:* fast template-based scoring to get strong candidates/prior info
#       *Stage 2:* EM + HMM joint inference (BAF + z-score) for final segmentation
#
# (B) Avoiding over-flexibility: dispersion is already controlled upstream
#     - We already have a prior step (standardization) designed to control the
#       dispersion / scale of the BAF-derived signal and/or related statistics.
#     - If we allow very flexible mixture fitting (e.g., freely adapting peak
#       widths, overdispersion, and weights) at this stage, the model can start
#       explaining away artifacts rather than reflecting true CN structure.
#     - By keeping templates more constrained (and only allowing limited tuning
#       via small grids for bw and uniform noise), we preserve identifiability
#       and keep this step aligned with the standardized data assumptions.
#
# Practical outcome:
#   - We still adopt the nQuack *principle* of comparing multiple candidate models,
#     but we implement it as a computationally cheap template-selection problem.
#   - We defer full parameter estimation (EM) to the joint EM+HMM stage, where it
#     is most valuable because it integrates multiple signals (BAF + z-score).
#
# ------------------------------------------------------------------------------
# 5) How our approach borrows the nQuack idea while remaining template-based
# ------------------------------------------------------------------------------
#
# Inspired by nQuack’s model selection across distributions, we extend our
# template scoring to test:
#   - multiple distribution families (gaussian/beta/beta-binomial/negative-binomial)
#   - a small grid of bw values (smoothing/dispersion proxy)
#   - optional uniform noise weights (off-comb contamination)
#
# Then we choose the best (dist, CN, bw, uniform_weight) combination by BIC.
#
# This captures the nQuack spirit (compare models; pick best) without reproducing
# full EM fitting at this stage.
#
# ------------------------------------------------------------------------------


if(getRversion() >= "2.15.1") utils::globalVariables(c(
  "BAF", "Density", "Type", "label", "CN"
))

#' Compute BAF Likelihoods
#'
#' Given a vector of BAF values and a vector of tested copy numbers (cn_grid),
#' this function creates BAF templates for each CN, computes the likelihood of the observed BAF
#' histogram under each template, and returns a vector of likelihoods (one per CN state).
#'
#' @param baf_vec Numeric vector of BAF values for a single window.
#' @param cn_grid Integer vector of copy-number states to test (e.g., 2:6).
#' @param M Integer. Number of BAF histogram bins on [0,1]. Default 100.
#' @param bw Numeric. Bandwidth (SD) of the Gaussian kernels used to generate the BAF comb templates. Default 0.03.
#' @param dist Character. Distribution type for template peaks: "gaussian" (default), "beta", "beta_binomial", or "negative_binomial".
#' @param reflect Logical. If TRUE (default), apply reflection for continuous
#'   kernels to keep mass within [0,1] (useful for gaussian).
#' @param add_uniform Logical. If TRUE, add a uniform noise component to the
#'   mixture before renormalization. Default FALSE.
#' @param uniform_weight Numeric in [0,1]. Mixture weight of the uniform
#'   component when \code{add_uniform = TRUE}. Default 0.05.
#' @param min_het_frac Numeric in [0,1]. Threshold for the fraction of BAF values in \code{het_range} considered heterozygous. If the observed heterozygous fraction exceeds this value, CN=1 is excluded from \code{cn_grid}, as a meaningful proportion of heterozygous loci makes haploid (CN=1) implausible. Default 0.05.
#' @param het_range Numeric vector of length 2. BAF interval used to define heterozygous loci (default \code{c(0.2, 0.8)}). Values outside this range are treated as homozygous for the purpose of the \code{min_het_frac} filter.
#' @param plot Logical. If TRUE, return a ggplot object showing observed and template distributions with likelihood and probability text for each CN. Default FALSE.
#'
#' @return A list containing:
#'   \item{ll_vec}{Vector of log-likelihoods (n_cn), one per CN template.}
#'   \item{prob_vec}{Vector of probabilities (n_cn), softmax-transformed from log-likelihoods.}
#'   \item{plot}{Optional ggplot object if plot=TRUE.}
#' @details Uses the generate_baf_template and baf_log_likelihood helpers from Qploidy.
#'
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#'
#' @export
#' @author Cristiane Taniguti
compute_baf_likelihoods <- function(baf_vec, cn_grid, M = 100, bw = 0.03,
                                   plot = FALSE, dist="gaussian", reflect = TRUE,
                                   add_uniform = FALSE, uniform_weight = 0.05,
                                   min_het_frac = 0.05,
                                   het_range = c(0.2, 0.8)) {

  # Exclude CN=1 from the grid when the data has sufficient heterozygosity (het_frac > min_het_frac),
  # as a meaningful proportion of heterozygous loci makes haploid (CN=1) implausible.
  usable_baf <- baf_vec[!is.na(baf_vec) & baf_vec >= 0 & baf_vec <= 1]
  no_one <- FALSE
  het_frac <- NA_real_
  if (length(usable_baf) > 0) {
    het_frac <- mean(usable_baf >= het_range[1] & usable_baf <= het_range[2])
    if (het_frac > min_het_frac) {
      if(any(cn_grid %in% 1L)) no_one <- TRUE
      cn_grid <- cn_grid[!cn_grid %in% 1L]
      if (length(cn_grid) == 0) stop(
        sprintf("All CN states removed after heterozygous fraction filter (het_frac=%.3f > min_het_frac=%.3f). Consider expanding cn_grid or lowering min_het_frac.", het_frac, min_het_frac)
      )
    }
  }

  K <- length(cn_grid)
  breaks <- seq(0, 1, length.out = M + 1)
  if (length(baf_vec) == 0 || all(is.na(baf_vec))) {
    hist_counts <- rep(0, M)
  } else {
    hist_counts <- hist(baf_vec, breaks = breaks, plot = FALSE)$counts
  }
  # Use correct dist spelling for generate_baf_template
  templates <- lapply(cn_grid, function(c) generate_baf_template(c, M = M, bw = bw, dist = dist[1], reflect = reflect,
                                                         add_uniform = add_uniform, uniform_weight = uniform_weight))
  names(templates) <- as.character(cn_grid)
  templates <- lapply(templates, function(t) { t <- pmax(t, 1e-8); t / sum(t) })
  ll_vec <- vapply(seq_len(K), function(k) {
    templ <- templates[[as.character(cn_grid[k])]]
    baf_log_likelihood(hist_counts, templ)
  }, numeric(1))
  prob_vec <- {
    exp_x <- exp(ll_vec - max(ll_vec, na.rm = TRUE))
    exp_x / sum(exp_x, na.rm = TRUE)
  }
  plot_obj <- NULL
  if (plot) {
    obs_density <- hist_counts / sum(hist_counts)
    df <- data.frame(
      BAF = seq(0, 1, length.out = M),
      Observed = obs_density
    )
    for (k in seq_len(K)) {
      df[[paste0("CN", cn_grid[k])]] <- templates[[as.character(cn_grid[k])]]
    }
    df_long <- pivot_longer(df, cols = -BAF, names_to = "Type", values_to = "Density")
    # Ensure Observed is a factor with levels matching the fill values
    df_long$Type <- factor(df_long$Type, levels = c("Observed", paste0("CN", cn_grid)))
    legend_labels <- c("Observed (black line)", paste0("CN", cn_grid, " (logLik=", sprintf("%.2f", ll_vec), ", P=", sprintf("%.2f", prob_vec), ")"))
    names(legend_labels) <- c("Observed", paste0("CN", cn_grid))
    plot_obj <- ggplot(df_long, aes(x = BAF, y = Density)) +
      geom_area(data = subset(df_long, Type != "Observed"), aes(fill = Type), position = "identity", alpha = 0.4) +
      geom_line(data = subset(df_long, Type == "Observed"), aes(color = "Observed", linetype = "Observed"), linewidth = 1) +
      scale_fill_manual(
        values = setNames(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")[seq_len(min(8, K))], paste0("CN", cn_grid)),
        breaks = paste0("CN", cn_grid),
        labels = legend_labels[paste0("CN", cn_grid)],
        drop = FALSE
      ) +
      scale_color_manual(
        values = c("Observed" = "black"),
        labels = c("Observed (black line)"),
        breaks = c("Observed")
      ) +
      scale_linetype_manual(
        values = c("Observed" = "solid"),
        labels = c("Observed (black line)"),
        breaks = c("Observed")
      ) +
      theme_bw() +
      labs(y = "Density", fill = "", color = "", linetype = "", title = paste(dist,"distribution")) +
      theme(legend.position = "top", legend.box = "vertical") +
      guides(
        fill = guide_legend(order = 2, ncol = 3, byrow = TRUE),
        color = guide_legend(order = 1, override.aes = list(linetype = "solid", linewidth = 1)),
        linetype = guide_legend(order = 1)
      )
  }
  if(no_one){
    out <- list(
      ll_vec = c(min(ll_vec),ll_vec),
      prob_vec = c(0, prob_vec)
    )
  } else {
    out <- list(
      ll_vec = ll_vec,
      prob_vec = prob_vec
    )
  }
  if (!is.null(plot_obj)) out$plot <- plot_obj
  return(out)
}

#' Select the Best BAF Template Model via Grid Search and BIC
#'
#' This function performs a grid search over candidate BAF template models, varying distribution
#' family, kernel bandwidth, and optional uniform noise component, to select the best model for a
#' given BAF vector. For each configuration, it computes the BAF likelihoods across candidate copy
#' number (CN) states, selects the best CN, and evaluates model fit using the Bayesian Information
#' Criterion (BIC). The function returns the best model configuration, a summary of all grid results,
#' and (optionally) a plot for the best model.
#'
#' @param qploidy_standardization An object of class \code{qploidy_standardization} containing standardized BAF values and related info. The function will extract the relevant BAF vector for the specified sample from this object. Either this (together with \code{sample}) or \code{baf_vec} must be provided.
#' @param sample Character. Sample name. Required when \code{qploidy_standardization} is provided; may be omitted when \code{baf_vec} is supplied directly.
#' @param baf_vec Optional numeric vector of BAF values. If provided, \code{qploidy_standardization} and \code{sample} may be omitted. Mutually exclusive with the \code{qploidy_standardization}/\code{sample} pair: if \code{baf_vec} is \code{NULL}, the BAF vector is extracted from the standardization object.
#' @param cn_grid Integer vector. Copy-number states to test (e.g., 2:6).
#' @param dists Character vector. Distribution families to test (e.g., c("gaussian", "beta", ...)).
#' @param reflect Logical. If TRUE (default), applies reflection for continuous kernels to keep mass within [0,1]. Passed to `generate_baf_template()`.
#' @param bw_grid Numeric vector. Bandwidth (SD or concentration) values to try for kernel smoothing. Default: c(0.02, 0.03, 0.04).
#' @param add_uniform_grid Logical vector. Whether to include a uniform noise component in the template. Default: c(FALSE, TRUE).
#' @param uniform_weight_grid Numeric vector. Mixture weights for the uniform component (when enabled). Default: c(0.01, 0.03, 0.05, 0.10, 0.15).
#' @param M Integer. Number of histogram bins for BAF (default 100).
#' @param plot Logical. If TRUE, returns a ggplot object for the best model only (default FALSE).
#' @param param_count Optional named integer vector. Number of free parameters per distribution (for BIC penalty). If NULL, defaults to 0 for all.
#' @param count_grid_as_params Logical. If TRUE (default), adds +1 to BIC penalty for each hyperparameter tuned by grid search (bw, and uniform_weight if used).
#' @param min_het_frac Numeric in [0,1]. Threshold for the fraction of BAF values in \code{het_range} considered heterozygous. If the observed heterozygous fraction exceeds this value, CN=1 is excluded from \code{cn_grid}, as a meaningful proportion of heterozygous loci makes haploid (CN=1) implausible. Default 0.05.
#' @param het_range Numeric vector of length 2. BAF interval used to define heterozygous loci (default \code{c(0.2, 0.8)}). Values outside this range are treated as homozygous for the purpose of the \code{min_het_frac} filter.
#'
#' @return A list with the following elements:
#'   \item{n_obs}{Number of usable BAF observations (after filtering NA/out-of-range).}
#'   \item{best}{A list describing the best model configuration: \code{dist}, \code{bw}, \code{add_uniform}, \code{uniform_weight}, \code{best_cn}, \code{logLik}, \code{prob}, \code{BIC}, \code{n_obs}, \code{het_frac} (observed heterozygous fraction), and \code{cn_grid_used} (CN states actually evaluated after filtering).}
#'   \item{grid_results}{A data.frame summarizing all grid search results (one row per configuration).}
#'   \item{plot}{(If plot=TRUE) ggplot object for the best model.}
#'   \item{note}{(If no usable data or all models fail) diagnostic message.}
#'
#' @details
#' For each combination of distribution, bandwidth, and uniform noise settings, the function:
#'   - Builds BAF templates for each CN state using `generate_baf_template()`.
#'   - Computes the multinomial log-likelihood of the observed BAF histogram under each template using `baf_log_likelihood()`.
#'   - Selects the best CN (highest log-likelihood) for that configuration.
#'   - Computes the BIC for the configuration: \eqn{\mathrm{BIC} = -2 \log L + p \log n}, where \eqn{p} is the number of free parameters (including grid-tuned ones if `count_grid_as_params=TRUE`), and \eqn{n} is the number of usable BAF values.
#'   - Returns a summary table of all configurations and highlights the best (lowest BIC).
#'   - If `plot=TRUE`, only the best model's plot is generated and returned.
#'
#' This approach is inspired by the model selection strategy in the nQuack package (Gaynor et al. 2024), but uses fixed templates and grid search for efficiency.
#'
#' @seealso \code{compute_baf_likelihoods}, \code{generate_baf_template}, \code{baf_log_likelihood}
#'
#' @references
#' Gaynor, M., Landis, J., O'Connor, T., Laport, R., Doyle, J., Soltis, D., Ponciano, J., & Soltis, P. (2024). "nQuack: An R package for predicting ploidal level from sequence data using site-based heterozygosity." Applications in Plant Sciences, 12(4), e11606. doi:10.1002/aps3.11606
#'
#' @export
#' @author Cristiane Taniguti
select_best_baf_model <- function(
  qploidy_standardization = NULL,
  sample = NULL,
  baf_vec = NULL,
  cn_grid,
  dists = c("gaussian", "beta", "beta_binomial", "negative_binomial"),
  reflect = TRUE,
  bw_grid = c(0.02, 0.03, 0.04),
  add_uniform_grid = FALSE,
  uniform_weight_grid = c(0.01, 0.03, 0.05, 0.10, 0.15),
  M = 100,
  plot = FALSE,
  param_count = NULL,
  count_grid_as_params = TRUE,
  min_het_frac = 0.05,
  het_range = c(0.2, 0.8)
) {
  # --- Input validation ---
  if (is.null(baf_vec)) {
    # Extract baf_vec from qploidy_standardization + sample
    if (!inherits(qploidy_standardization, "qploidy_standardization"))
      stop("'qploidy_standardization' must be an object of class 'qploidy_standardization' (or provide 'baf_vec' directly).")
    if (is.null(qploidy_standardization$data) || !is.data.frame(qploidy_standardization$data))
      stop("'qploidy_standardization$data' must be a data.frame.")
    required_cols <- c("SampleName", "baf")
    missing_cols <- setdiff(required_cols, names(qploidy_standardization$data))
    if (length(missing_cols) > 0)
      stop(sprintf("'qploidy_standardization$data' is missing required columns: %s", paste(missing_cols, collapse = ", ")))

    if (!is.character(sample) || length(sample) != 1)
      stop("'sample' must be a single character string when 'qploidy_standardization' is used.")

    available_samples <- unique(qploidy_standardization$data$SampleName)
    if (!sample %in% available_samples)
      stop(sprintf(
        "Sample '%s' not found in the dataset. Available samples are:\n  %s",
        sample, paste(available_samples, collapse = ", ")
      ))

    one_sample <- qploidy_standardization$data[qploidy_standardization$data$SampleName == sample, ]
    baf_vec <- one_sample$baf
  } else {
    # baf_vec provided directly; qploidy_standardization and sample may be omitted
    if (!is.numeric(baf_vec) || length(baf_vec) == 0)
      stop("'baf_vec' must be a non-empty numeric vector.")
  }

  stopifnot(is.numeric(cn_grid) || is.integer(cn_grid))
  cn_grid <- as.integer(cn_grid)
  if (length(cn_grid) < 1 || any(cn_grid < 1L))
    stop("'cn_grid' must be a non-empty integer vector with all values >= 1.")

  stopifnot(is.numeric(het_range), length(het_range) == 2, all(is.finite(het_range)))
  if (het_range[1] >= het_range[2])
    stop("'het_range[1]' must be strictly less than 'het_range[2]'.")

  stopifnot(is.numeric(min_het_frac), length(min_het_frac) == 1, min_het_frac >= 0, min_het_frac <= 1)

  allowed_dists <- c("gaussian", "beta", "beta_binomial", "negative_binomial")
  if (!is.character(dists) || length(dists) < 1)
    stop("'dists' must be a non-empty character vector.")
  invalid_dists <- setdiff(dists, allowed_dists)
  if (length(invalid_dists) > 0)
    stop(sprintf("Invalid distribution(s) in 'dists': %s. Allowed values: %s.",
                 paste(invalid_dists, collapse = ", "), paste(allowed_dists, collapse = ", ")))

  if (!is.logical(reflect) || length(reflect) != 1)
    stop("'reflect' must be a single logical value (TRUE or FALSE).")
  if (!is.logical(plot) || length(plot) != 1)
    stop("'plot' must be a single logical value (TRUE or FALSE).")
  if (!is.logical(count_grid_as_params) || length(count_grid_as_params) != 1)
    stop("'count_grid_as_params' must be a single logical value (TRUE or FALSE).")

  if (!is.numeric(M) || length(M) != 1 || M < 2 || M != as.integer(M))
    stop("'M' must be a single integer >= 2.")
  M <- as.integer(M)

  stopifnot(is.numeric(bw_grid), all(is.finite(bw_grid)), all(bw_grid > 0))
  stopifnot(is.logical(add_uniform_grid), length(add_uniform_grid) >= 1)
  stopifnot(is.numeric(uniform_weight_grid),
            all(is.finite(uniform_weight_grid)),
            all(uniform_weight_grid >= 0),
            all(uniform_weight_grid <= 1))

  stopifnot(is.numeric(baf_vec))

  # n = number of usable observations for BIC penalty
  usable_baf <- baf_vec[!is.na(baf_vec) & baf_vec >= 0 & baf_vec <= 1]
  n_obs <- length(usable_baf)

  # Default parameter counts per distribution
  D <- length(dists)
  if (is.null(param_count)) {
    param_count <- setNames(rep.int(0L, D), dists)
  } else {
    stopifnot(is.numeric(param_count) || is.integer(param_count))
    if (is.null(names(param_count))) stop("param_count must be a *named* vector keyed by dists.")
    missing <- setdiff(dists, names(param_count))
    if (length(missing) > 0) {
      param_count <- c(param_count, setNames(rep.int(0L, length(missing)), missing))
    }
    param_count <- as.integer(param_count[dists])
  }

  # If no data, return NA-ish structure early
  if (n_obs == 0L) {
    return(list(
      n_obs = n_obs,
      best = NULL,
      grid_results = data.frame(),
      note = "No usable BAF observations (n_obs = 0)."
    ))
  }

  # Build grid of configs
  # - if add_uniform=FALSE, we only evaluate uniform_weight=0 (to avoid duplicates)
  grid <- expand.grid(
    dist = dists,
    bw = bw_grid,
    add_uniform = add_uniform_grid,
    uniform_weight = uniform_weight_grid,
    stringsAsFactors = FALSE
  )
  grid <- grid[!(grid$add_uniform == FALSE & grid$uniform_weight != uniform_weight_grid[1]), , drop = FALSE]
  # Replace remaining uniform_weight with 0 when add_uniform=FALSE
  grid$uniform_weight[grid$add_uniform == FALSE] <- 0

  # Remove exact duplicate rows (in case uniform_weight_grid[1] != 0)
  grid <- unique(grid)

  # Helpers
  bic_from_ll <- function(logLik, p, n) {
    # BIC = -2 logLik + p log(n)
    -2 * logLik + p * log(n)
  }

  # Store per-grid best CN results
  grid_results <- grid
  grid_results$best_cn <- NA_integer_
  grid_results$logLik <- NA_real_
  grid_results$prob <- NA_real_
  grid_results$BIC <- NA_real_

  # Iterate over grid
  for (g in seq_len(nrow(grid))) {
    dist_g <- grid$dist[g]
    bw_g <- grid$bw[g]
    add_u <- grid$add_uniform[g]
    uw_g <- grid$uniform_weight[g]

    # Always call with plot=FALSE for efficiency
    res <- compute_baf_likelihoods(
      baf_vec = baf_vec,
      cn_grid = cn_grid,
      M = M,
      bw = bw_g,
      plot = FALSE,
      dist = dist_g,
      reflect = reflect,
      add_uniform = add_u,
      uniform_weight = uw_g,
      het_range = het_range,
      min_het_frac = min_het_frac
    )

    ll <- res$ll_vec
    pr <- res$prob_vec

    if (all(is.na(ll))) next
    idx <- which.max(ll)

    best_cn <- cn_grid[idx]
    best_ll <- ll[idx]
    best_pr <- pr[idx]

    # Parameter penalty
    p <- param_count[dist_g]

    if (isTRUE(count_grid_as_params)) {
      # You are effectively selecting bw by search -> count +1
      p <- p + 1L
      # If uniform component is enabled, and you are also searching its weight -> +1
      if (isTRUE(add_u)) p <- p + 1L
    }

    bic <- bic_from_ll(best_ll, p = p, n = n_obs)

    grid_results$best_cn[g] <- best_cn
    grid_results$logLik[g] <- best_ll
    grid_results$prob[g] <- best_pr
    grid_results$BIC[g] <- bic
  }

  # Pick global best row (min BIC)
  if (all(is.na(grid_results$BIC))) {
    return(list(
      n_obs = n_obs,
      best = NULL,
      grid_results = grid_results,
      note = "All model evaluations returned NA likelihoods."
    ))
  }

  best_idx <- which.min(grid_results$BIC)
  best_row <- grid_results[best_idx, , drop = FALSE]

  best <- list(
    dist = best_row$dist,
    bw = best_row$bw,
    add_uniform = best_row$add_uniform,
    uniform_weight = best_row$uniform_weight,
    best_cn = best_row$best_cn,
    logLik = best_row$logLik,
    prob = best_row$prob,
    BIC = best_row$BIC,
    n_obs = n_obs,
    min_het_frac = min_het_frac,
    cn_grid_used = cn_grid
  )

  out <- list(
    n_obs = n_obs,
    best = best,
    grid_results = grid_results[order(grid_results$BIC), , drop = FALSE]
  )

  # Only generate the plot for the best model if requested
  if (plot) {
    best_plot_res <- compute_baf_likelihoods(
      baf_vec = baf_vec,
      cn_grid = cn_grid,
      M = M,
      bw = best$bw,
      plot = TRUE,
      dist = best$dist,
      reflect = reflect,
      add_uniform = best$add_uniform,
      uniform_weight = best$uniform_weight,
      min_het_frac = min_het_frac,
      het_range = het_range
    )
    out$plot <- best_plot_res$plot
  }

  if (length(unique(out$grid_results$best_cn[1:3])) > 1) {
    warn_msg <- if (!is.null(sample))
      paste0("Sample ", sample, ": Top 3 models differ on estimated CN. Inspect plot for quality control.")
    else
      "Top 3 models differ on estimated CN. Inspect plot for quality control."
    warning(warn_msg)
  }

  class(out) <- "selected_BAF_model"
  return(out)
}

#' Build a polyploid BAF comb template for a given copy number
#'
#' Constructs a discrete BAF density template on [0,1] for a given copy number
#' \code{c}, by placing kernels at genotype fractions \eqn{d/c} for
#' \eqn{d=0,...,c}. Binomial weights are used for relative peak heights.
#'
#' Optionally, the template can be extended with a uniform noise component to
#' capture "off-comb" BAF values (e.g., mapping bias, noisy loci).
#'
#' @param cn Integer. Copy number (e.g., 2 for diploid, 4 for tetraploid).
#' @param M Integer. Number of bins over [0,1]. Default is 101.
#' @param bw Numeric. Bandwidth / concentration control (interpretation depends
#'   on \code{dist}). For \code{gaussian}, it is sd on the BAF scale.
#'   For \code{beta}, it controls concentration around the mode.
#'   For \code{beta_binomial} and \code{negative_binomial}, it controls optional
#'   smoothing applied on the x-grid.
#' @param floor_eps Numeric. Small positive value added before renormalization.
#' @param dist Character. One of \code{"gaussian"}, \code{"beta"},
#'   \code{"beta_binomial"}, \code{"negative_binomial"}.
#' @param reflect Logical. If TRUE (default), apply reflection for continuous
#'   kernels to keep mass within [0,1] (useful for gaussian).
#' @param add_uniform Logical. If TRUE, add a uniform noise component to the
#'   mixture before renormalization. Default FALSE.
#' @param uniform_weight Numeric in [0,1]. Mixture weight of the uniform
#'   component when \code{add_uniform = TRUE}. Default 0.05.
#'
#' @return Numeric vector of length \code{M} summing to 1.
#'
#' @export
#' @importFrom stats dnorm dbeta dbinom dnbinom
#' 
#' @author Cristiane Taniguti
generate_baf_template <- function(cn, M = 101, bw = 0.03, floor_eps = 1e-8,
                         dist = c("gaussian","beta","beta_binomial","negative_binomial"),
                         reflect = TRUE,
                         add_uniform = FALSE,
                         uniform_weight = 0.05) {

  dist <- match.arg(dist)
  stopifnot(length(cn) == 1, is.finite(cn), cn >= 1)
  cn <- as.integer(cn)
  stopifnot(length(M) == 1, M >= 2)
  M <- as.integer(M)

  if (isTRUE(add_uniform)) {
    stopifnot(is.finite(uniform_weight), uniform_weight >= 0, uniform_weight <= 1)
  }

  x <- seq(0, 1, length.out = M)

  # genotype centers exactly at d/cn
  d <- 0:cn
  centers <- d / cn

  # relative heights (user can swap later if desired)
  wd <- dbinom(d, size = cn, prob = 0.5)
  wd <- wd / sum(wd)

  # --- helper: add reflection for continuous kernels
  add_reflection <- function(v, mu, kernel_fun) {
    if (!isTRUE(reflect)) return(kernel_fun(mu))
    kernel_fun(mu) + kernel_fun(-mu) + kernel_fun(2 - mu)
  }

  # --- continuous kernels
  if (dist == "gaussian") {
    stopifnot(is.finite(bw), bw > 0)

    kern <- function(mu) dnorm(x, mean = mu, sd = bw)
    dens <- Reduce(`+`, Map(function(mu, w) w * add_reflection(x, mu, kern), centers, wd))

  } else if (dist == "beta") {
    kappa <- max(2, round(1 / max(bw, 1e-6)))
    kern <- function(mu) {
      if (mu <= 0) {
        a <- 1; b <- 1 + kappa
      } else if (mu >= 1) {
        a <- 1 + kappa; b <- 1
      } else {
        a <- 1 + mu * kappa
        b <- 1 + (1 - mu) * kappa
      }
      dbeta(x, shape1 = a, shape2 = b)
    }
    dens <- Reduce(`+`, Map(function(mu, w) w * kern(mu), centers, wd))

  } else if (dist %in% c("beta_binomial","negative_binomial")) {

    grid_idx_for_k <- vapply(0:cn, function(k) which.min(abs(x - (k / cn))), integer(1))

    smooth_on_grid <- function(p) {
      if (!is.finite(bw) || bw <= 0) return(p)
      sd_bins <- bw / (1 / (M - 1))
      r <- max(1L, as.integer(ceiling(4 * sd_bins)))
      offs <- (-r):r
      g <- dnorm(offs, mean = 0, sd = sd_bins)
      g <- g / sum(g)
      p2 <- numeric(length(p))
      for (i in seq_along(p)) {
        j <- i + offs
        ok <- j >= 1 & j <= length(p)
        p2[i] <- sum(p[j[ok]] * g[ok])
      }
      p2
    }

    make_spike_pmf <- function(p_k) {
      p <- numeric(M)
      for (k in 0:cn) {
        p[grid_idx_for_k[k + 1]] <- p[grid_idx_for_k[k + 1]] + p_k[k + 1]
      }
      smooth_on_grid(p)
    }

    if (dist == "beta_binomial") {
      phi <- max(2, round(1 / max(bw, 1e-6)))

      one_kernel <- function(d0) {
        if (d0 == 0) {
          alpha <- 1; beta <- 1 + phi
        } else if (d0 == cn) {
          alpha <- 1 + phi; beta <- 1
        } else {
          alpha <- 1 + (d0 / cn) * phi
          beta <- 1 + (1 - d0 / cn) * phi
        }
        k <- 0:cn
        logpmf <- lchoose(cn, k) + lbeta(k + alpha, cn - k + beta) - lbeta(alpha, beta)
        p_k <- exp(logpmf - max(logpmf))
        p_k <- p_k / sum(p_k)
        make_spike_pmf(p_k)
      }

      dens <- Reduce(`+`, Map(function(d0, w) w * one_kernel(d0), d, wd))

    } else { # negative_binomial
      size <- max(1, cn)

      one_kernel <- function(d0) {
        k <- 0:cn
        if (d0 == 0) {
          p_k <- numeric(cn + 1); p_k[1] <- 1
          return(make_spike_pmf(p_k))
        }
        p <- size / (size + d0)
        p_k <- dnbinom(k, size = size, prob = p)
        p_k <- p_k / sum(p_k)
        make_spike_pmf(p_k)
      }

      dens <- Reduce(`+`, Map(function(d0, w) w * one_kernel(d0), d, wd))
    }
  }

  # --- add uniform noise component (optional)
  if (isTRUE(add_uniform) && uniform_weight > 0) {
    u <- rep(1 / M, M)
    dens <- (1 - uniform_weight) * dens + uniform_weight * u
  }

  dens <- dens + floor_eps
  dens <- dens / sum(dens)
  dens
}


#' BAF histogram log-likelihood under a template
#'
#' Computes the multinomial-style log-likelihood of observed BAF histogram counts against a template (probabilities over the same bins). Only bins with nonzero counts contribute. If all counts are zero, returns 0.
#'
#' @param counts Integer vector of length \code{M}. BAF bin counts (e.g., from \code{hist(..., plot = FALSE)$counts}).
#' @param templ Numeric vector of length \code{M}. Strictly positive probabilities summing to 1 (e.g., output of \code{generate_baf_template}).
#' @param eps Numeric. Small positive value added inside \code{log()} for extra numerical safety. Default is 1e-8.
#'
#' @return Numeric scalar. The log-likelihood \eqn{\sum_i n_i \log p_i}.
#'
#' @details
#' Only bins with \code{counts > 0} are used. If all counts are zero, returns 0. Ensure \code{templ} has no zeros (use \code{floor_eps} in \code{generate_baf_template}).
#'
#' @examples
#' templ <- generate_baf_template(4, M = 21)
#' counts <- rmultinom(1, size = 200, prob = templ)[,1]
#' baf_log_likelihood(counts, templ)
#'
#' @export
#' @author Cristiane Taniguti
baf_log_likelihood <- function(counts, templ, eps=1e-8) {
  idx <- counts > 0L
  if (!any(idx)) return(0)  # empty histogram contributes nothing
  sum(counts[idx] * log(templ[idx] + eps))
}


#' Print method for selected_BAF_model objects
#'
#' Prints a concise summary of the best BAF model selected by
#' \code{select_best_baf_model()}, including the winning configuration,
#' goodness-of-fit statistics, and a short overview of the full grid search.
#'
#' @param x An object of class \code{'selected_BAF_model'}.
#' @param ... Additional arguments (ignored).
#'
#' @method print selected_BAF_model
#'
#' @export
print.selected_BAF_model <- function(x, ...) {
  if (!inherits(x, "selected_BAF_model"))
    stop("Object is not of class 'selected_BAF_model'.")

  cat("selected_BAF_model\n")
  cat("  Usable BAF observations:", x$n_obs, "\n")

  if (!is.null(x$note)) {
    cat("  Note:", x$note, "\n")
    return(invisible(x))
  }

  b <- x$best
  cat("\nBest model\n")
  cat("  CN:               ", b$best_cn, "\n")
  cat("  Distribution:     ", b$dist, "\n")
  cat("  Bandwidth (bw):   ", b$bw, "\n")
  cat("  Uniform component:", b$add_uniform)
  if (isTRUE(b$add_uniform))
    cat(sprintf("  (weight = %.3f)", b$uniform_weight))
  cat("\n")
  cat("  log-Likelihood:   ", round(b$logLik, 3), "\n")
  cat("  Probability:      ", round(b$prob, 4), "\n")
  cat("  BIC:              ", round(b$BIC, 3), "\n")

  cat("\nGrid search summary (", nrow(x$grid_results), " configurations evaluated)\n", sep = "")
  gr <- x$grid_results[!is.na(x$grid_results$BIC), ]
  if (nrow(gr) > 0) {
    cn_tab <- sort(table(gr$best_cn), decreasing = TRUE)
    cat("  CN votes across grid:", paste(names(cn_tab), cn_tab, sep = "x", collapse = "  "), "\n")
    cat("  BIC range: [", round(min(gr$BIC), 2), ",", round(max(gr$BIC), 2), "]\n")
    top3 <- head(gr[order(gr$BIC), c("dist", "bw", "add_uniform", "best_cn", "BIC")], 3)
    cat("  Top 3 configurations:\n")
    for (i in seq_len(nrow(top3))) {
      cat(sprintf("    %d. CN=%d  dist=%-20s  bw=%.3f  uniform=%s  BIC=%.2f\n",
                  i, top3$best_cn[i], top3$dist[i], top3$bw[i],
                  top3$add_uniform[i], top3$BIC[i]))
    }
  }

  if (!is.null(x$plot)) print(x$plot)

  invisible(x)
}

