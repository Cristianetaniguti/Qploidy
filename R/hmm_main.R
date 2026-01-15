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
#' @param segment_zscore Logical. If \code{TRUE}, segment z-scores using changepoint detection to define windows instead of fixed SNP counts. Default \code{TRUE}.
#' @param snps_per_window Integer. Number of SNPs per window when windows are created. Default \code{500}.
#' @param reflect Logical. If TRUE (default), apply reflection for continuous
#'   kernels to keep mass within [0,1] (useful for gaussian).
#' @param add_uniform Logical. If TRUE, add a uniform noise component to the
#'   mixture before renormalization. Default FALSE.
#' @param uniform_weight Numeric in [0,1]. Mixture weight of the uniform
#'   component when \code{add_uniform = TRUE}. Default 0.05.
#' @param min_snps_per_window Integer. Minimum SNPs required to keep a window (windows with fewer SNPs are dropped). Default \code{100}.
#' @param cn_grid Integer vector of copy-number states to consider (e.g., \code{2:8}).
#' @param M Integer. Number of BAF histogram bins on [0,1]. Default \code{121}.
#' @param max_iter Integer. Maximum EM iterations. Default \code{60}.
#' @param het_lims Numeric vector of length 2. BAF limits to consider a SNP heterozygous. Default \code{c(0,1)}.
#' @param het_quantile Numeric. Quantile used to scale BAF emission weight based on heterozygote count. Default \code{0.8}.
#' @param baf_weight Numeric. Overall weight applied to BAF emission (0–1). Default \code{1}.
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
#' This function calls \code{\link{generate_baf_template}}, \code{\link{baf_log_likelihood}}, \code{\link{logsumexp}}, and \code{\link{viterbi}} which must be available in your package/namespace.
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
#' @seealso \code{\link{generate_baf_template}}, \code{\link{baf_log_likelihood}},
#'   \code{\link{viterbi}}, \code{\link{plot_cn_track}}
#' @importFrom stats dnorm quantile sd
#' @importFrom graphics hist
#' @importFrom utils tail
#' @export
hmm_estimate_CN <- function(
    qploidy_standarize_result,
    sample_id,
    chr = NULL,
    reflect = TRUE,
    add_uniform = FALSE,
    uniform_weight = 0.05,
    segment_zscore = TRUE,
    snps_per_window = 500,
    min_snps_per_window = 20,
    cn_grid = 2:8,
    M = 100,
    max_iter = 60,
    het_lims = c(0,1), # baf limits to consider a SNP heterozygous
    het_quantile = 0.8, # increase this value to reduce the weight of baf when few hets
    baf_weight = 1,
    z_range = NULL, 
    transition_jump = 0.995, # decrease this value if you think there changes in CN is likely
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

  df <- as.data.frame(qploidy_standarize_result$data)
  if (nrow(df) == 0) stop("Input data is empty. No SNPs found for any sample.")

  # subset to one sample
  if (verbose) cat("\nSubsetting by sample and chromosomes...\n")
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

  # Remove markers with missing baf and z
  rm <- which(is.na(d[["baf"]]) & is.na(d[["z"]]))
  if(length(rm) > 0) d <- d[-rm, , drop = FALSE]

  # If z_range is not provided, estimate from data
  if (is.null(z_range) || (length(z_range) == 1 && is.na(z_range))) {
    z_range <- (1/length(cn_grid)) * (max(d$z, na.rm = TRUE) - min(d$z, na.rm = TRUE))
    if (verbose) cat(sprintf("    Estimated z_range from data: %f\n", z_range))
  }

  # --- set expected ploidy ---
  # Calculate expected ploidy using sample-level BAF distribution
  if(verbose) cat("Selecting best BAF model to set expected ploidy based on sample-level BAF distribution...\n")
  selected_model <- select_best_baf_model(d$baf,
                                          cn_grid= cn_grid,
                                          M = M,
                                          reflect = reflect)

  if(verbose) {
    cat("  Best BAF model summary:\n")
    cat(sprintf(
      "    CN: %s\n    Distribution: %s\n    Bandwidth: %.4f\n    Uniform noise: %s (weight: %.2f)\n    Log-likelihood: %.2f\n    BIC: %.2f\n    Probability: %.3f\n    Observations: %d\n",
      selected_model$best$best_cn,
      selected_model$best$dist,
      selected_model$best$bw,
      ifelse(selected_model$best$add_uniform, "yes", "no"),
      selected_model$best$uniform_weight,
      selected_model$best$logLik,
      selected_model$best$BIC,
      selected_model$best$prob,
      selected_model$best$n_obs
    ))
  }

  # --- build windows ---
  if (verbose) cat("Building windows...\n")
  d <- d[order(d[["Chr"]], d[["Position"]]), ]

  # Segmented z-score
  if (segment_zscore){
    if (verbose) cat("  Segmenting z-scores to define windows.\n")
    d <- add_changepoint_windows(d, minseglen = min_snps_per_window)

  } else {
    # simple fixed-size windows
    if (verbose) cat("  Building fixed-size windows.\n")
    # equal-SNP windows within chromosome
    d$.__w__ <- with(
      d,
      ave(
        seq_len(nrow(d)),
        d[["Chr"]],
        FUN = function(k) {
          n <- length(k)

          # initial windowing (same idea as before)
          w <- ceiling(seq_along(k) / snps_per_window)

          # if there is a remainder AND we have more than one full window,
          # merge the last (small) window into the previous one
          if (n > snps_per_window && n %% snps_per_window != 0L) {
            last <- max(w)
            w[w == last] <- last - 1L
          }
          w
        }
      )
    )
  }
  win_col <- ".__w__"

  # summarize per window
  if (verbose) cat("Summarizing data by window...\n")
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

  if(any(is.nan(win_df$z_mean))) {
    warning("Some windows have no z-score values.")
  }

  if (is.null(win_df) || nrow(win_df) == 0) stop("No windows could be formed. Check input data and window parameters.")

  # drop very small windows - mostly for the last window
  keep <- win_df$n_snps >= min_snps_per_window
  if (!any(keep)) stop("All windows dropped by min_snps_per_window filter. Try lowering min_snps_per_window or check input data.")

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

  # --- BAF histograms and likelihoods ---
  breaks <- seq(0, 1, length.out = M + 1)
  hist_counts <- matrix(0L, nrow = W, ncol = M)
  for (i in seq_len(W)) {
    if (length(baf_list[[i]]) == 0 || all(is.na(baf_list[[i]]))) {
      warning(sprintf("Window %d (Chr %s, WindowID %s) has no valid BAF values for histogram.",
                      i, win_df$Chr[i], win_df$WindowID[i]))
    }
    hist_counts[i, ] <- hist(baf_list[[i]], breaks = breaks, plot = FALSE)$counts
  }

  # Generate BAF likelihoods and probabilities per window
  # Uses parameters from selected_model
  if(!z_only){
    if(verbose) cat("Generating BAF likelihoods per window...\n")
    baf_results <- lapply(baf_list, function(baf_vec) compute_baf_likelihoods(baf_vec,
                                                                              cn_grid,
                                                                              M = M,
                                                                              bw = selected_model$best$bw,
                                                                              plot = FALSE,
                                                                              dist = selected_model$best$dist,
                                                                              reflect = reflect,
                                                                              add_uniform = selected_model$best$add_uniform,
                                                                              uniform_weight = selected_model$best$uniform_weight))

    ll_baf_matrix <- do.call(rbind, lapply(baf_results, function(res) res$ll_vec))
    prob_baf_matrix <- do.call(rbind, lapply(baf_results, function(res) res$prob_vec))

    if(verbose) cat("Counting heterozygous to determine BAF weights...\n")
    # BAF weights from heterozygote counts
    # Lower the number of heterozygotes in the window, lower the weight of the BAF emission
    # If there are no heterozygous, only z-score will be considered
    HET_N <- win_df$n_het
    if (any(HET_N > 0)) {
      ref_q <- suppressWarnings(quantile(HET_N[HET_N>0], het_quantile, na.rm=TRUE))
      if (!is.finite(ref_q) || ref_q <= 0) ref_q <- max(HET_N, na.rm=TRUE)
      if (!is.finite(ref_q) || ref_q <= 0) ref_q <- 1
      w_baf <- sqrt(pmin(1, HET_N / ref_q))  # sqrt: temper influence
    } else {
      w_baf <- rep(0, W)
    }
    w_baf[!is.finite(w_baf)] <- 0
    w_baf <- w_baf * baf_weight
  } else {
    ll_baf_matrix <- matrix(0, nrow = W, ncol = length(cn_grid))
    prob_baf_matrix <- matrix(0, nrow = W, ncol = length(cn_grid))
    w_baf <- rep(0, W)
  }

  # --- HMM by chromosome ---
  results_list <- list()
  params_list <- list()
  for (chr_name in unique(win_df$Chr)) {
    if (verbose) cat(sprintf("\nFitting HMM for chromosome %s\n", chr_name))

    # Split by chromosome
    chr_idx <- which(win_df$Chr == chr_name)
    W_chr <- length(chr_idx)
    win_df_chr <- win_df[chr_idx, ]
    z_chr <- z[chr_idx]
    baf_list_chr <- baf_list[chr_idx]
    ll_baf_matrix_chr <- ll_baf_matrix[chr_idx, , drop = FALSE]
    prob_baf_matrix_chr <- prob_baf_matrix[chr_idx, , drop = FALSE]
    w_baf_chr <- w_baf[chr_idx]

    # HMM setup for this chromosome
    K_chr <- length(cn_grid) # The number of HMM states are the number of ploidies tested
    state_ids_chr <- as.character(cn_grid)

    # transitions
    # A is K×K. Before normalizing, every off-diagonal entry is 0.001 (rare jumps), and the diagonal is 0.995 (strong tendency to stay in the same CN state).
    if(verbose) cat("  Setting transition matrices.\n")
    A_chr <- matrix(1e-3, K_chr, K_chr); diag(A_chr) <- transition_jump
    # A <- A / rowSums(A) makes each row sum to 1 (a valid Markov matrix). After this, values are almost unchanged (just scaled so each row sums exactly to 1).
    A_chr <- A_chr / rowSums(A_chr)
    # To consider: If is possible to change this matrix to favor specific small jumps, like if changing from 2->4 is more likely than 2->6
    # pi0 is the initial state distribution at the first window; here it’s uniform (no prior preference for any CN at the start).
    pi0_chr <- rep(1/K_chr, K_chr)

    if(verbose) cat("  Defining z-score distribution templates.\n")
    # z mean init (monotone ramp)
    # z where is the mean z per window
    # this subtracts 0.2 and add 0.2 to the min and max of z, and splits it in K values
    # The small padding ±0.2 keeps edge states from starting exactly at the extremes, which helps EM avoid collapsing to boundary values.
    # Because CN is ordered (2 < 3 < 4 …) and z increases with CN, this guarantees an ordered starting guess for means.
    # the maximum value + 0.2 will be referring to the highest CN state - what I am not sure if it is a correct assumption
    # the mean z value will reffer to the expected ploidy provided
    z_mean_chr  <- mean(z_chr, na.rm = TRUE)
    z_lo_chr    <- min(z_chr, na.rm = TRUE) - z_range
    z_hi_chr    <- max(z_chr, na.rm = TRUE) + z_range
    cmin_chr <- min(cn_grid)
    cmax_chr <- max(cn_grid)

    # Recompute chromosome-level expected ploidy
    chromosome_level_ploidy <- select_best_baf_model(unlist(baf_list_chr)[-which(is.na(unlist(baf_list_chr)))],
                                                     cn_grid,
                                                     M = M,
                                                     bw = selected_model$best$bw,
                                                     plot = FALSE,
                                                     dist = selected_model$best$dist,
                                                     reflect = reflect,
                                                     add_uniform = selected_model$best$add_uniform,
                                                     uniform_weight = selected_model$best$uniform_weight)

    exp_ploidy_chr <- chromosome_level_ploidy$best$best_cn
    # If for some reason exp_ploidy_chr is NULL or non-finite, fall back to sample-level exp_ploidy
    if (is.null(exp_ploidy_chr) || !is.finite(exp_ploidy_chr)) {
      exp_ploidy_chr <- exp_ploidy
    }

    # Choose a single linear step so extremes fit within [z_lo, z_hi]
    # We need:
    #   z_mean + step*(cmin - exp_ploidy) <= z_lo
    #   z_mean + step*(cmax - exp_ploidy) >= z_hi
    # Solve for step and take the max magnitude to satisfy both.
    step_lo_chr <- if (exp_ploidy_chr > cmin_chr) (z_mean_chr - z_lo_chr) / (exp_ploidy_chr - cmin_chr) else 0
    step_hi_chr <- if (exp_ploidy_chr < cmax_chr) (z_hi_chr  - z_mean_chr) / (cmax_chr - exp_ploidy_chr) else 0
    step_chr    <- max(step_lo_chr, step_hi_chr, 1e-6)

    # Monotone, baseline-centered initialization:
    mu_vec_chr <- z_mean_chr + step_chr * (as.numeric(cn_grid) - exp_ploidy_chr)
    mu_chr     <- setNames(mu_vec_chr, as.character(cn_grid))
    # If state_ids are strings of cn_grid, this aligns. If not, reorder:
    mu_chr <- mu_chr[state_ids_chr]

    # sig is the (shared) standard deviation of the z emission across states.
    # It starts at the sample SD of z, with a safety floor of 0.1 to avoid zero/near-zero variance that would blow up log-likelihoods.
    sig_chr <- sd(z_chr, na.rm = TRUE); if (!is.finite(sig_chr) || sig_chr <= 1e-6) sig_chr <- 0.1

    # EM loop for this chromosome
    # --- EM loop ---
    ll_hist_chr <- numeric(max_iter)
    # If there is only one window, skip HMM and assign CN by BAF likelihood only
    if (W_chr == 1) {
      if(verbose) cat("  Only one window found for this chromosome. Skipping EM and assigning CN by BAF likelihood only.\n")
      cn_call_chr <- cn_grid[apply(ll_baf_matrix_chr, 1, which.max)]
      post_max_chr <- apply(prob_baf_matrix_chr, 1, max)
      post_df_chr <- as.data.frame(prob_baf_matrix_chr)
      names(post_df_chr) <- paste0("post_CN", cn_grid)
      result_chr <- cbind(
        data.frame(
          Sample    = sample_id,
          Chr       = win_df_chr$Chr,
          WindowID  = win_df_chr$WindowID,
          Start     = win_df_chr$Start,
          End       = win_df_chr$End,
          n_snps    = win_df_chr$n_snps,
          n_het     = win_df_chr$n_het,
          z         = z_chr,
          w_baf     = if(z_only) 0 else 1,  # when only one window, baf weight is always 1
          CN_call   = cn_call_chr,
          post_max  = post_max_chr,
          stringsAsFactors = FALSE
        ),
        post_df_chr
      )
      params_chr <- list(
        cn_grid = cn_grid,
        mu = NA,
        sigma = NA,
        A = NA,
        pi0 = NA,
        bins = M,
        bw = selected_model$best$bw,
        loglik = NA
      )
      results_list[[chr_name]] <- result_chr
      params_list[[chr_name]] <- params_chr
      next
    }
    # Otherwise, run EM/HMM as usual
    if(verbose) cat("  Starting EM.\n")
    for (iter in 1:max_iter) {
      # Emissions
      ll_em_chr <- matrix(NA_real_, nrow=W_chr, ncol=K_chr, dimnames=list(NULL, state_ids_chr))
      for (k in seq_len(K_chr)) {
        c <- cn_grid[k]
        llz <- dnorm(z_chr, mean=mu_chr[as.character(c)], sd=sig_chr, log=TRUE)
        if(any(is.nan(llz))) llz[which(is.nan(llz))] <- 0
        if (z_only) {
          ll_em_chr[,k] <- llz
        } else {
          llb <- ll_baf_matrix_chr[,k]
          ll_em_chr[, k] <- llz + w_baf_chr * llb
        }
      }

      if (!all(is.finite(ll_em_chr))) {
        bad_w <- which(!is.finite(rowSums(ll_em_chr)))[1]
        bad_k <- which(!is.finite(ll_em_chr[bad_w, ]))
        stop(sprintf("Non-finite emission at window %d, states: %s.",
                     bad_w, paste(colnames(ll_em_chr)[bad_k], collapse=", ")))
      }

      logA <- log(A_chr); logpi0 <- log(pi0_chr)
      log_alpha <- matrix(-Inf, W_chr, K_chr); log_beta <- matrix(0, W_chr, K_chr)

      # forward
      log_alpha[1, ] <- logpi0 + ll_em_chr[1, ]
      for (i in 2:W_chr) {
        for (k in 1:K_chr) {
          log_alpha[i,k] <- ll_em_chr[i,k] + logsumexp(log_alpha[i-1, ] + logA[,k])
        }
      }
      # backward
      for (i in (W_chr-1):1) {
        for (k in 1:K_chr) {
          log_beta[i,k] <- logsumexp(logA[k, ] + ll_em_chr[i+1, ] + log_beta[i+1, ])
        }
      }
      loglik <- logsumexp(log_alpha[W_chr, ])
      ll_hist_chr[iter] <- loglik

      # E-step: posteriors
      log_gamma <- log_alpha + log_beta
      log_gamma <- sweep(log_gamma, 1, apply(log_gamma, 1, logsumexp), "-")
      gamma <- exp(log_gamma)

      # pairwise
      xi_sum <- matrix(0, K_chr, K_chr)
      for (i in 1:(W_chr-1)) {
        M_ij <- outer(log_alpha[i, ], log_beta[i+1, ], "+") +
          logA + matrix(ll_em_chr[i+1, ], K_chr, K_chr, byrow=TRUE)
        M_ij <- M_ij - logsumexp(as.vector(M_ij))
        xi_sum <- xi_sum + exp(M_ij)
      }

      # M-step: update parameters
      pi0_chr <- gamma[1, ] / sum(gamma[1, ])
      A_chr <- xi_sum / pmax(rowSums(xi_sum), 1e-12)
      A_chr[!is.finite(A_chr)] <- 0
      A_chr <- sweep(A_chr, 1, pmax(rowSums(A_chr), 1e-12), "/")
      A_chr <- pmax(A_chr, 1e-12); A_chr <- sweep(A_chr, 1, rowSums(A_chr), "/")

      # update mu and sigma
      mu_chr <- numeric(K_chr)
      for (k in 1:K_chr) {
        w <- gamma[,k]
        mu_chr[k] <- sum(w * z_chr) / pmax(sum(w), 1e-12)
      }
      mu_chr <- setNames(mu_chr, as.character(cn_grid))

      # update shared sigma
      sig_chr <- sqrt(sum(gamma * (matrix(z_chr, W_chr, K_chr) - rep(mu_chr, each=W_chr))^2) /
                        pmax(sum(gamma), 1e-12))
      sig_chr <- max(sig_chr, 1e-3)

      # Convergence check
      if (iter > 4 && is.finite(ll_hist_chr[iter]) && is.finite(ll_hist_chr[iter-1]) &&
          abs(ll_hist_chr[iter] - ll_hist_chr[iter-1]) < 1e-4) break
    }

    if(verbose) cat(sprintf("  EM converged in %d iterations. Final log-likelihood: %.2f\n", iter, ll_hist_chr[iter]))

    # Decode Viterbi path
    vit_path_chr <- viterbi(ll_em_chr, log(A_chr), log(pi0_chr))
    cn_call_chr <- cn_grid[vit_path_chr]

    # max posterior per window
    post_max_chr <- apply(gamma, 1, max)

    # Prepare output
    post_df_chr <- as.data.frame(gamma)
    names(post_df_chr) <- paste0("post_CN", cn_grid)
    result_chr <- cbind(
      data.frame(
        Sample    = sample_id,
        Chr       = win_df_chr$Chr,
        WindowID  = win_df_chr$WindowID,
        Start     = win_df_chr$Start,
        End       = win_df_chr$End,
        n_snps    = win_df_chr$n_snps,
        n_het     = win_df_chr$n_het,
        z         = z_chr,
        w_baf     = if(z_only) 0 else w_baf_chr,
        CN_call   = cn_call_chr,
        post_max  = post_max_chr,
        stringsAsFactors = FALSE
      ),
      post_df_chr
    )
    params_chr <- list(
      cn_grid = cn_grid,
      mu = mu_chr,
      sigma = sig_chr,
      A = A_chr,
      pi0 = pi0_chr,
      bins = M,
      bw = selected_model$best$bw,
      loglik = tail(ll_hist_chr[is.finite(ll_hist_chr)], 1)
    )
    results_list[[chr_name]] <- result_chr
    params_list[[chr_name]] <- params_chr
  }
  # Combine results
  combined_result <- do.call(rbind, results_list)
  rownames(combined_result) <- NULL
  combined_params <- params_list
  if(verbose) cat("\nDone!\n")
  return(structure(list(result = combined_result, params = combined_params), class = "hmm_CN"))
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
                varlist = c("worker", "hmm_estimate_CN", "generate_baf_template", "baf_log_likelihood", "logsumexp", "viterbi"),
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


