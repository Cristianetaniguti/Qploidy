#' Polyploid copy-number HMM (z + BAF) per sample and window
#'
#' Aggregates SNP-level data into windows, builds BAF histograms, and fits a
#' Hidden Markov Model that combines a Gaussian emission for depth/z and a
#' polyploid BAF “comb” emission to infer copy number (CN) per window.
#' Returns a tidy data frame with CN calls and posterior probabilities that is
#' ready for plotting (e.g., with \code{\link{plot_cn_track}}) or merging across
#' chromosomes/batches.
#'
#' @param qploidy_standarize_result An object of class \code{qploidy_standardization} as returned by \code{standardize()}. Contains standardized SNP-level data for all samples and chromosomes.
#' @param sample_id Character scalar. Sample identifier to analyze; rows are filtered as \code{SampleName == sample_id}.
#' @param chr Optional. Character or integer vector specifying chromosomes to include. If \code{NULL}, all chromosomes are used.
#' @param snps_per_window Integer. Number of SNPs per window when windows are created. Default \code{500}.
#' @param min_snps_per_window Integer. Minimum SNPs required to keep a window (windows with fewer SNPs are dropped). Default \code{100}.
#' @param cn_grid Integer vector of copy-number states to consider (e.g., \code{2:8}).
#' @param M Integer. Number of BAF histogram bins on [0,1]. Default \code{121}.
#' @param bw Numeric. Bandwidth (SD) of the Gaussian kernels used to generate the BAF comb templates. Default \code{0.03}.
#' @param max_iter Integer. Maximum EM iterations. Default \code{60}.
#' @param het_lims Numeric vector of length 2. BAF limits to consider a SNP heterozygous. Default \code{c(0,1)}.
#' @param het_weight Numeric. Quantile used to scale BAF emission weight based on heterozygote count. Default \code{0.8}.
#' @param z_range Numeric. Padding added to min/max z for initial mean estimation. Default \code{0.2}.
#' @param transition_jump Numeric. Diagonal value for transition matrix (probability to stay in same CN state). Default \code{0.995}.
#' @param exp_ploidy Integer. Expected ploidy for initialization. If \code{NULL}, median of \code{cn_grid} is used.
#' @param z_only Logical. If \code{TRUE}, fit the HMM using the z-emission only (ignores BAF). Default \code{FALSE}.
#' @param verbose Logical. If \code{TRUE}, print progress messages. Default \code{TRUE}.
#'
#' @return An object of class \code{hmm_CN}, a list with two elements:
#'   \describe{
#'     \item{result}{A \code{data.frame} with one row per window and columns:
#'       \itemize{
#'         \item \code{Sample}, \code{Chr}, \code{WindowID}, \code{Start}, \code{End}
#'         \item \code{n_snps}, \code{n_het} – counts per window
#'         \item \code{z} – mean z per window
#'         \item \code{w_baf} – weight applied to the BAF emission (0–1)
#'         \item \code{CN_call} – Viterbi copy-number call (factor/integer)
#'         \item \code{post_max} – posterior probability of the called state
#'         \item \code{post_CN<k>} – posterior probability for each state in \code{cn_grid}
#'       }
#'     }
#'     \item{params}{A list of fitted HMM parameters:
#'       \itemize{
#'         \item \code{cn_grid} – copy-number states
#'         \item \code{mu} – estimated z means per state
#'         \item \code{sigma} – estimated z standard deviation
#'         \item \code{A} – transition matrix
#'         \item \code{pi0} – initial state probabilities
#'         \item \code{bins} – number of BAF bins (M)
#'         \item \code{bw} – BAF template bandwidth
#'         \item \code{loglik} – final log-likelihood
#'       }
#'     }
#'   }
#' The \code{result} data frame is suitable for downstream plotting or merging. The \code{params} attribute contains all HMM model parameters for reproducibility and inspection.
#'
#' @details
#' \strong{Windowing.} If \code{chr} is specified, only those chromosomes are analyzed. Windows are formed within each chromosome by grouping every \code{snps_per_window} consecutive SNPs. Windows with fewer than \code{min_snps_per_window} SNPs are removed.
#'
#' \strong{Emissions.} For state \code{c} (copy number), the z-emission is \eqn{\mathcal{N}(\mu_c, \sigma^2)}. The BAF emission is a multinomial log-likelihood comparing the observed BAF histogram to a state-specific template (a “comb” with peaks at \eqn{d/c}). Windows with few heterozygotes are down-weighted via \code{w_baf} (derived from \code{n_het}).
#'
#' \strong{Numerical safety.} The implementation guards against non-finite emissions, enforces stochastic \code{A} rows, and lower-bounds \code{sigma}.
#'
#' @section Required helpers:
#' This function calls \code{\link{baf_template}}, \code{\link{baf_ll}}, \code{\link{logsumexp}}, and \code{\link{viterbi}} which must be available in your package/namespace.
#'
#' @examples
#' \dontrun{
#' # Suppose qploidy_standarize_result is from standardize(),
#' # with columns: SampleName, Chr, Position, baf, z
#' res <- hmm_estimate_CN(
#'   qploidy_standarize_result = qploidy_standarize_result,
#'   sample_id = "ASample",
#'   chr = c("Chr1", "Chr2"),
#'   snps_per_window = 500,
#'   cn_grid = 2:6
#' )
#' head(res$result)
#' res$params
#' }
#'
#' @seealso \code{\link{baf_template}}, \code{\link{baf_ll}},
#'   \code{\link{viterbi}}, \code{\link{plot_cn_track}}
#' @importFrom stats dnorm quantile sd
#' @importFrom graphics hist
#' @importFrom utils tail
#' @export
hmm_estimate_CN <- function(
    qploidy_standarize_result,
    sample_id,
    chr = NULL,
    snps_per_window = 500,
    min_snps_per_window = 20,
    cn_grid = 2:8,
    M = 100,
    bw = 0.03,
    max_iter = 60,
    het_lims = c(0,1), # baf limits to consider a SNP heterozygous
    het_weight = 0.8, # increase this value to reduce the weight of baf when few hets
    z_range = 0.2, # increase this value if you think that extreme ploidy tested are unlikely
    transition_jump = 0.995, # decrease this value if you think there changes in CN is likely
    exp_ploidy = NULL, # if NULL, median(cn_grid) is used
    z_only = FALSE,
    verbose = TRUE
) {

  # --- input checks ---
  if (!is(qploidy_standarize_result, "qploidy_standardization")) {
    stop("Input must be a qploidy_standardization object as returned by standardize().")
  }

  if (!is.character(sample_id) || length(sample_id) != 1 || nchar(sample_id) == 0) {
    stop("sample_id must be a non-empty character scalar.")
  }

  if (is.null(exp_ploidy)) exp_ploidy <- median(cn_grid)

  df <- as.data.frame(qploidy_standarize_result$data)
  if (nrow(df) == 0) stop("Input data is empty. No SNPs found for any sample.")

  # subset to one sample
  if (verbose) cat("Subsetting by sample and chromosomes.\n")
  d <- df[df[["SampleName"]] == sample_id, , drop = FALSE]
  if (nrow(d) == 0) stop(sprintf("No data found for sample '%s'.", sample_id))

  # subset to chromosomes
  if (!is.null(chr)) {
    if (is.numeric(chr)) {
      chrs <- unique(d[["Chr"]])[chr]
      d <- d[which(d[["Chr"]] %in% chrs), , drop = FALSE]
    } else {
      d <- d[which(d[["Chr"]] %in% chr), , drop = FALSE]
    }
    if (nrow(d) == 0) stop("No data found for specified chromosomes after filtering.")
  }

  # --- build windows ---
  d <- d[order(d[["Chr"]], d[["Position"]]), ]

  # equal-SNP windows within chromosome
  d$.__w__ <- with(d, ave(seq_len(nrow(d)), d[["Chr"]],
                          FUN = function(k) ceiling(seq_along(k) / snps_per_window)))
  win_col <- ".__w__"

  # summarize per window
  if (verbose) cat("Summarizing by window.\n")
  agg <- within(d, {
    is_het <- baf > het_lims[1] & baf < het_lims[2]
  })

  win_df <- do.call(rbind, by(agg, list(agg[["Chr"]], agg[[win_col]]), function(x) {
    if (is.null(x) || nrow(x) == 0) return(NULL)
    data.frame(
      Chr        = x[1, "Chr"],
      WindowID   = x[1, win_col],
      Start      = min(x[["Position"]]),
      End        = max(x[["Position"]]),
      n_snps     = nrow(x),
      n_het      = sum(x$is_het, na.rm = TRUE),
      z_mean     = mean(x[["z"]], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(win_df) <- NULL
  if (is.null(win_df) || nrow(win_df) == 0) stop("No windows could be formed. Check input data and window parameters.")

  # drop very small windows - mostly for the last window
  keep <- win_df$n_snps >= min_snps_per_window
  if (!any(keep)) stop("All windows dropped by min_snps_per_window filter. Try lowering min_snps_per_window or check input data.")
  if (!all(keep)) win_df <- win_df[keep, ]
  if (nrow(win_df) == 0) stop("All windows dropped by min_snps_per_window filter. Try lowering min_snps_per_window or check input data.")

  # order windows first by chr then start
  o <- order(win_df$Chr, win_df$Start)
  win_df <- win_df[o, ]
  W <- nrow(win_df) # W = number of windows

  # extract BAF vectors per window (list) and z vector
  baf_list <- vector("list", W)
  for (i in seq_len(W)) {
    wrows <- d[["Chr"]] == win_df$Chr[i] &
      d[[win_col]]  == win_df$WindowID[i] &
      d[["Position"]] >= win_df$Start[i] &
      d[["Position"]] <= win_df$End[i]
    baf_list[[i]] <- d[wrows, "baf"]
    if (length(baf_list[[i]]) == 0 || all(is.na(baf_list[[i]]))) {
      warning(sprintf("Window %d (Chr %s, WindowID %s) has no valid BAF values.",
                   i, win_df$Chr[i], win_df$WindowID[i]))
    }
  }
  z <- win_df$z_mean
  if (length(z) == 0 || all(is.na(z))) stop("No valid z values found in any window.")

  # BAF histograms
  breaks <- seq(0, 1, length.out = M + 1)
  hist_counts <- matrix(0L, nrow = W, ncol = M)
  for (i in seq_len(W)) {
    if (length(baf_list[[i]]) == 0 || all(is.na(baf_list[[i]]))) {
      warning(sprintf("Window %d (Chr %s, WindowID %s) has no valid BAF values for histogram.",
                   i, win_df$Chr[i], win_df$WindowID[i]))
    }
    hist_counts[i, ] <- hist(baf_list[[i]], breaks = breaks, plot = FALSE)$counts
  }

  if(verbose) cat("Building BAF distributions templates.\n")
  # Templates for BAF emissions
  # This creates a list of BAF templates, one per CN state in cn_grid according
  # to a multi normal distribution with means on each ploidy peak and sd = bw
  # The values are spread over M bins
  templates <- lapply(cn_grid, function(c) baf_template(c, M=M, bw=bw))
  names(templates) <- as.character(cn_grid)
  # normalize and floor templates to avoid log(0)
  templates <- lapply(templates, function(t) { t <- pmax(t, 1e-8); t / sum(t) })

  if(verbose) cat("Counting heterozygous.\n")
  # BAF weights from heterozygote counts
  # Lower the number of heterozygotes in the window, lower the weight of the BAF emission
  # If there are no heterozygous, only z-score will be considered
  HET_N <- win_df$n_het
  if (any(HET_N > 0)) {
    ref_q <- suppressWarnings(quantile(HET_N[HET_N>0], het_weight, na.rm=TRUE))
    if (!is.finite(ref_q) || ref_q <= 0) ref_q <- max(HET_N, na.rm=TRUE)
    if (!is.finite(ref_q) || ref_q <= 0) ref_q <- 1
    w_baf <- sqrt(pmin(1, HET_N / ref_q))  # sqrt: temper influence
  } else {
    w_baf <- rep(0, W)
  }
  w_baf[!is.finite(w_baf)] <- 0

  # --- HMM setup ---
  K <- length(cn_grid) # The number of HMM states are the number of ploidies tested
  state_ids <- as.character(cn_grid)

  if(verbose) cat("Setting transition matrices.\n")
  # transitions
  # A is K×K. Before normalizing, every off-diagonal entry is 0.001 (rare jumps), and the diagonal is 0.995 (strong tendency to stay in the same CN state).
  A <- matrix(1e-3, K, K); diag(A) <- transition_jump
  # A <- A / rowSums(A) makes each row sum to 1 (a valid Markov matrix). After this, values are almost unchanged (just scaled so each row sums exactly to 1).
  A <- A / rowSums(A)
  # To consider: If is possible to change this matrix to favor specific small jumps, like if changing from 2->4 is more likely than 2->6
  # pi0 is the initial state distribution at the first window; here it’s uniform (no prior preference for any CN at the start).
  pi0 <- rep(1/K, K)

  if(verbose) cat("Defining z-score distribution templates.\n")

  # z mean init (monotone ramp)
  # z where is the mean z per window
  # this subtracts 0.2 and add 0.2 to the min and max of z, and splits it in K values
  # The small padding ±0.2 keeps edge states from starting exactly at the extremes, which helps EM avoid collapsing to boundary values.
  # Because CN is ordered (2 < 3 < 4 …) and z increases with CN, this guarantees an ordered starting guess for means.
  # the maximum value + 0.2 will be referring to the highest CN state - what I am not sure if it is a correct assumption
  # the mean z value will reffer to the expected ploidy provided
  z_mean  <- mean(z, na.rm = TRUE)
  z_lo    <- min(z, na.rm = TRUE) - z_range
  z_hi    <- max(z, na.rm = TRUE) + z_range

  cmin <- min(cn_grid)
  cmax <- max(cn_grid)

  # Sanity: clamp exp_ploidy to grid if needed
  if (!exp_ploidy %in% cn_grid) {
    warning("exp_ploidy not in cn_grid; using nearest grid value.")
    exp_ploidy <- cn_grid[which.min(abs(cn_grid - exp_ploidy))]
  }

  # Choose a single linear step so extremes fit within [z_lo, z_hi]
  # We need:
  #   z_mean + step*(cmin - exp_ploidy) <= z_lo
  #   z_mean + step*(cmax - exp_ploidy) >= z_hi
  # Solve for step and take the max magnitude to satisfy both.
  step_lo <- if (exp_ploidy > cmin) (z_mean - z_lo) / (exp_ploidy - cmin) else 0
  step_hi <- if (exp_ploidy < cmax) (z_hi  - z_mean) / (cmax - exp_ploidy) else 0
  step    <- max(step_lo, step_hi, 1e-6)   # avoid zero step

  # Monotone, baseline-centered initialization:
  mu_vec <- z_mean + step * (as.numeric(cn_grid) - exp_ploidy)
  mu     <- setNames(mu_vec, as.character(cn_grid))  # names must match state_ids

  # If state_ids are strings of cn_grid, this aligns. If not, reorder:
  mu <- mu[state_ids]

  # sig is the (shared) standard deviation of the z emission across states.
  #It starts at the sample SD of z, with a safety floor of 0.1 to avoid zero/near-zero variance that would blow up log-likelihoods.
  sig <- sd(z, na.rm = TRUE); if (!is.finite(sig) || sig <= 1e-6) sig <- 0.1

  # --- EM loop ---
  if(verbose) cat("Starting EM.\n")

  ll_hist <- numeric(max_iter)
  for (iter in 1:max_iter) {
    # emissions
    ll_em <- matrix(NA_real_, nrow=W, ncol=K, dimnames=list(NULL, state_ids))
    for (k in seq_len(K)) {
      c <- cn_grid[k]
      llz <- dnorm(z, mean=mu[as.character(c)], sd=sig, log=TRUE) # Values from a normal distribution for z considering ploidy/state mean
      if (z_only) {
        ll_em[,k] <- llz
      } else {
        templ <- templates[[as.character(c)]]
        llb <- apply(hist_counts, 1, baf_ll, templ=templ)
        ll_em[,k] <- llz + w_baf * llb # here is how the emission likelihood is being calculated
      }
    }
    if (!all(is.finite(ll_em))) {
      bad_w <- which(!is.finite(rowSums(ll_em)))[1]
      bad_k <- which(!is.finite(ll_em[bad_w, ]))
      stop(sprintf("Non-finite emission at window %d, states: %s.",
                   bad_w, paste(colnames(ll_em)[bad_k], collapse=", ")))
    }

    logA <- log(A); logpi0 <- log(pi0)
    log_alpha <- matrix(-Inf, W, K); log_beta <- matrix(0, W, K)

    # forward
    log_alpha[1, ] <- logpi0 + ll_em[1, ]
    for (i in 2:W) {
      for (k in 1:K) {
        log_alpha[i,k] <- ll_em[i,k] + logsumexp(log_alpha[i-1, ] + logA[,k])
      }
    }
    # backward
    for (i in (W-1):1) {
      for (k in 1:K) {
        log_beta[i,k] <- logsumexp(logA[k, ] + ll_em[i+1, ] + log_beta[i+1, ])
      }
    }
    loglik <- logsumexp(log_alpha[W, ])
    ll_hist[iter] <- loglik

    # posteriors
    log_gamma <- log_alpha + log_beta
    log_gamma <- sweep(log_gamma, 1, apply(log_gamma, 1, logsumexp), "-")
    gamma <- exp(log_gamma)

    # pairwise
    xi_sum <- matrix(0, K, K)
    for (i in 1:(W-1)) {
      M_ij <- outer(log_alpha[i, ], log_beta[i+1, ], "+") +
        logA + matrix(ll_em[i+1, ], K, K, byrow=TRUE)
      M_ij <- M_ij - logsumexp(as.vector(M_ij))
      xi_sum <- xi_sum + exp(M_ij)
    }

    # M-step
    pi0 <- gamma[1, ] / sum(gamma[1, ])
    A <- xi_sum / pmax(rowSums(xi_sum), 1e-12)
    A[!is.finite(A)] <- 0
    A <- sweep(A, 1, pmax(rowSums(A), 1e-12), "/")
    A <- pmax(A, 1e-12); A <- sweep(A, 1, rowSums(A), "/")

    # update z means
    for (k in 1:K) {
      w <- gamma[,k]
      mu[k] <- sum(w * z) / pmax(sum(w), 1e-12)
    }
    # pooled sigma
    sig <- sqrt(sum(gamma * (matrix(z, W, K) - rep(mu, each=W))^2) /
                  pmax(sum(gamma), 1e-12))
    sig <- max(sig, 1e-3)

    if (iter > 4 && is.finite(ll_hist[iter]) && is.finite(ll_hist[iter-1]) &&
        abs(ll_hist[iter] - ll_hist[iter-1]) < 1e-4) break
  }

  if(verbose) cat("Finished EM.\n")

  if(verbose) cat("Getting path with Viterbi.\n")
  vit_path <- viterbi(ll_em, log(A), log(pi0))
  cn_call <- cn_grid[vit_path]

  # max posterior per window
  post_max <- apply(gamma, 1, max)

  # --- tidy output for plots/tables ---
  # posterior columns as wide format: post_CN2, post_CN3, ...
  post_df <- as.data.frame(gamma)
  names(post_df) <- paste0("post_CN", cn_grid)

  if(verbose) cat("Just organizing the results.\n")

  result <- cbind(
    data.frame(
      Sample    = sample_id,
      Chr       = win_df$Chr,
      WindowID  = win_df$WindowID,
      Start     = win_df$Start,
      End       = win_df$End,
      n_snps    = win_df$n_snps,
      n_het     = win_df$n_het,
      z         = z,
      w_baf     = if(z_only) 0 else w_baf,
      CN_call   = cn_call,
      post_max  = post_max,
      stringsAsFactors = FALSE
    ),
    post_df
  )

  params <- list(
    cn_grid = cn_grid,
    mu = mu,
    sigma = sig,
    A = A,
    pi0 = pi0,
    bins = M,
    bw = bw,
    loglik = tail(ll_hist[is.finite(ll_hist)], 1)
  )

  if(verbose) cat("Done!\n")

  return(structure(list(result =result, params = params), class = "hmm_CN"))
}

#' Run hmm_estimate_CN in parallel for multiple samples (using parLapply)
#'
#' Runs copy-number HMM for each sample in a user-defined vector (or all samples if "all" is specified).
#' Returns a combined data.frame with results for all samples.
#'
#' @param qploidy_standarize_result An object of class qploidy_standardization as returned by standardize().
#' @param sample_ids Character vector of sample IDs to analyze, or "all" for all samples in the data.
#' @param n_cores Number of cores to use (default: 2).
#' @param ... Additional arguments passed to hmm_estimate_CN (e.g., chr, snps_per_window, etc).
#' @return A data.frame with results for all samples, as returned by hmm_estimate_CN$result, combined.
#'
#' @importFrom parallel makeCluster parLapply stopCluster clusterExport
#'
#' @export
hmm_estimate_CN_multi <- function(qploidy_standarize_result,
                                  sample_ids = "all",
                                  n_cores = 2,
                                  ...) {
  # sanity check
  if (!inherits(qploidy_standarize_result, "qploidy_standardization")) {
    stop("Input must be a qploidy_standardization object as returned by standardize().")
  }

  df <- as.data.frame(qploidy_standarize_result$data)
  all_samples <- unique(df$SampleName)

  if (identical(sample_ids, "all")) {
    sample_ids <- all_samples
  } else {
    sample_ids <- intersect(sample_ids, all_samples)
    if (length(sample_ids) == 0) stop("No valid sample IDs found in input data.")
  }

  # capture dots ONCE
  dots <- list(...)

  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)

  clusterEvalQ(cl, {
    ## make errors surface promptly
    options(warn = 1)
    ## packages that define hmm_estimate_CN, classes, and helpers
    library(methods)    # important for S4 on PSOCK clusters
    library(Qploidy)
    # library(dplyr)
    # library(tidyr)
    NULL
  })

  # make sure workers can see hmm_estimate_CN (and any helpers it needs)
  clusterExport(cl,
                          varlist = c("worker", "hmm_estimate_CN", "baf_template", "baf_ll", "logsumexp", "viterbi"),
                          envir = environment()
  )

  results_list <- parLapply(cl, sample_ids, worker,
                                      qploidy_standarize_result, dots)

  parameters <- lapply(results_list, function(x) x$params)
  results_list <- lapply(results_list, function(x) x$result)
  results_list <- Filter(Negate(is.null), results_list)
  if (length(results_list) == 0) stop("No results returned for any sample.")

  combined <- do.call(rbind, results_list)
  rownames(combined) <- NULL
  return(structure(list(result =combined, params = parameters[[1]]), class = "hmm_CN"))
}


