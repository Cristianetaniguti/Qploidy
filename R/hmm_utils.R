##' @keywords internal
# Suppress global variable warnings for non-standard evaluation in dplyr/ggplot2
globalVariables(c(
  "Sample", "Chr", "Start", "End", "CN_call", "prob_call", "w_baf", "Mid",
  "n_vlines", "n_regions", "region_id"
))

#' Numerically stable log-sum-exp
#'
#' Computes the logarithm of the sum of exponentials, \eqn{\log\left(\sum_i e^{x_i}\right)}, in a numerically stable way by subtracting the maximum element before exponentiating. This avoids overflow/underflow for large or small values in \code{x}.
#'
#' @param x Numeric vector. Values to exponentiate and sum.
#'
#' @return Numeric scalar. The stable computation of \code{log(sum(exp(x)))}.
#'
#' @details
#' This function uses the identity \eqn{\log\sum_i e^{x_i} = m + \log\sum_i e^{x_i - m}}, where \eqn{m=\max_i x_i}, to ensure numerical stability.
#'
#'
#' @keywords internal
logsumexp <- function(x) {
  m <- max(x); m + log(sum(exp(x - m)))
}


#' Viterbi decoding for a first-order HMM in log-space
#'
#' Computes the most likely state sequence (Viterbi path) for a Hidden Markov Model given per-position emission log-likelihoods and log transition/prior probabilities. Operates in log-space for numerical stability.
#'
#' @param ll_em Numeric matrix (W x K). Emission log-likelihoods, where \code{W} is the number of positions/windows and \code{K} is the number of states. Column order must match \code{logA}.
#' @param logA Numeric matrix (K x K). Log transition probabilities, where \code{logA[i, j]} is \eqn{\log p(s_t=j | s_{t-1}=i)}. Each row should represent a valid log-probability distribution.
#' @param logpi0 Numeric vector of length \code{K}. Initial state log-probabilities at position 1.
#'
#' @return Integer vector of length \code{W}. 1-based indices of the most likely states at each position (Viterbi path).
#'
#' @details
#' Stores backpointers to reconstruct the optimal path after dynamic programming. All calculations are performed in log-space for stability.
#'
#' @examples
#' set.seed(1)
#' W <- 5; K <- 3
#' ll_em <- matrix(rnorm(W*K), W, K)
#' A <- matrix(1/K, K, K); diag(A) <- 0.8; A <- A / rowSums(A)
#' logA <- log(A)
#' logpi0 <- log(rep(1/K, K))
#' vpath <- viterbi(ll_em, logA, logpi0)
#' vpath
#'
#' @export
viterbi <- function(ll_em, logA, logpi0) {
  W <- nrow(ll_em); K <- ncol(ll_em)
  delta <- matrix(-Inf, W, K)
  psi <- matrix(NA_integer_, W, K)
  delta[1, ] <- logpi0 + ll_em[1, ]
  for (i in 2:W) {
    for (k in 1:K) {
      tmp <- delta[i-1, ] + logA[,k]
      psi[i, k] <- which.max(tmp)
      delta[i, k] <- ll_em[i, k] + max(tmp)
    }
  }
  path <- integer(W)
  path[W] <- which.max(delta[W, ])
  for (i in (W-1):1) path[i] <- psi[i+1, path[i+1]]
  path
}

#' Internal worker for parallel HMM CN estimation
#'
#' Runs hmm_estimate_CN for a single sample within a parallel loop, forwarding arguments.
#' Returns the result or NULL on error (with a warning).
#'
#' @param sid Character. Sample identifier.
#' @param obj Standardized input object for copy-number estimation (usually of class 'qploidy_standardization').
#' @param dots Named list of additional arguments to pass to hmm_estimate_CN.
#'
#' @return List as returned by hmm_estimate_CN, or NULL if an error occurs.
#'
#' @keywords internal
#' @noRd
worker <- function(sid, obj, dots, data = NULL, geno.pos = NULL, use_values = c("BAF", "zscore")) {
  collected_warnings <- character(0)
  result <- withCallingHandlers(
    tryCatch({
      if (!is.null(obj)) {
        do.call(hmm_estimate_CN,
                c(list(qploidy_standarize_result = obj, sample_id = sid, use_values = use_values), dots))
      } else {
        do.call(hmm_estimate_CN,
                c(list(qploidy_standarize_result = NULL, sample_id = sid, data = data, geno.pos = geno.pos, use_values = use_values), dots))
      }
    }, error = function(e) {
      collected_warnings <<- c(
        collected_warnings,
        sprintf("Sample '%s' failed: %s", sid, conditionMessage(e))
      )
      NULL
    }),
    warning = function(w) {
      collected_warnings <<- c(
        collected_warnings,
        sprintf("[Sample '%s'] %s", sid, conditionMessage(w))
      )
      invokeRestart("muffleWarning")
    }
  )
  list(result = result, warnings = collected_warnings)
}

#' Summarize copy number mode and posterior probability
#'
#' Summarizes the mode of copy number calls and the mean maximum posterior probability per sample,
#' chromosome, or chromosome arm. Handles input as a data.frame or hmm_CN object.
#'
#' @param df Data.frame or hmm_CN object containing windowed CN calls and posterior probabilities.
#' @param level Character. Summarization level: 'sample', 'chromosome', or 'chromosome-arm'.
#' @param centromeres Optional. Named numeric vector or data.frame with centromere positions (required for chromosome-arm).
#' @param cn_col Character. Column name for CN calls (default: 'CN_call').
#' @param post_col Character. Column name for maximum posterior probability (default: 'post_max').
#'
#' @return Data.frame with columns for sample, chromosome (and arm if requested), CN mode, mean max posterior, and window count.
#'
#'
#' @export
summarize_cn_mode <- function(df,
                              level = c("sample", "chromosome", "chromosome-arm"),
                              centromeres = NULL,
                              cn_col = "CN_call",
                              post_col = "post_max") {
  level <- match.arg(level)

  # unwrap if hmm_CN-like object with $by_window
  dat <- if (is.data.frame(df)) df else if (!is.null(df$by_window) && is.data.frame(df$by_window)) df$by_window else
    stop("Input must be a data.frame or an hmm_CN-like object with a data.frame at `$by_window`.")

  # required columns
  req_cols <- c("Sample", "Chr", "Start", "End", cn_col, post_col)
  miss <- setdiff(req_cols, names(dat))
  if (length(miss)) stop("Missing required columns in data: ", paste(miss, collapse = ", "))

  if (!is.numeric(dat[[cn_col]]) && !is.integer(dat[[cn_col]])) {
    stop(sprintf("Column '%s' must be numeric/integer.", cn_col))
  }
  if (!is.numeric(dat[[post_col]])) {
    stop(sprintf("Column '%s' must be numeric.", post_col))
  }

  dat <- as.data.frame(dat)

  # chromosome-arm handling → create chrID.1 (Chr) and chrID.2 (1=p, 2=q)
  if (level == "chromosome-arm") {
    if (is.null(centromeres))
      stop("For level='chromosome-arm', provide `centromeres` (named numeric vector or data.frame with Chr, Centromere).")

    if (is.data.frame(centromeres)) {
      if (!all(c("Chr", "Centromere") %in% names(centromeres)))
        stop("centromeres data.frame must have columns: Chr, Centromere")
      cm <- setNames(centromeres$Centromere, centromeres$Chr)
    } else if (is.numeric(centromeres) && !is.null(names(centromeres))) {
      cm <- centromeres
    } else {
      stop("centromeres must be a named numeric vector or a data.frame with Chr and Centromere.")
    }

    mid <- (dat$Start + dat$End) / 2
    has_cm <- dat$Chr %in% names(cm)

    # Assign arm index: 1 = p (left), 2 = q (right)
    chrID.2 <- rep(NA_integer_, nrow(dat))
    chrID.2[has_cm] <- ifelse(mid[has_cm] < cm[dat$Chr[has_cm]], 1L, 2L)

    if (any(!has_cm)) {
      missing_chr <- unique(dat$Chr[!has_cm])
      warning("No centromere provided for: ", paste(missing_chr, collapse = ", "),
              ". Rows for these chromosomes will be dropped.")
    }

    keep <- has_cm & !is.na(chrID.2)
    dat <- dat[keep, , drop = FALSE]
    chrID.2 <- chrID.2[keep]

    dat$chrID.1 <- dat$Chr
    dat$chrID.2 <- chrID.2
  }

  # grouping variables
  group_vars <- c("Sample")
  if (level %in% c("chromosome", "chromosome-arm")) group_vars <- c(group_vars, "Chr")
  if (level == "chromosome-arm") group_vars <- c(group_vars, "chrID.1", "chrID.2")

  # summarize (base R)
  split_idx <- interaction(dat[group_vars], drop = TRUE, lex.order = TRUE)
  grouped <- split(seq_len(nrow(dat)), split_idx)

  out <- lapply(grouped, function(idx) {
    g <- dat[idx, , drop = FALSE]
    row <- list(
      Sample        = g$Sample[1],
      CN_mode       = mode(g[[cn_col]]),
      mean_max_prob = mean(g[[post_col]], na.rm = TRUE),
      n_windows     = nrow(g)
    )
    if ("Chr" %in% group_vars)     row$Chr      <- g$Chr[1]
    if ("chrID.1" %in% group_vars) row$chrID.1  <- g$chrID.1[1]
    if ("chrID.2" %in% group_vars) row$chrID.2  <- g$chrID.2[1]
    as.data.frame(row, stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, out)

  # nice column ordering
  want <- c("Sample",
            if (level %in% c("chromosome", "chromosome-arm")) "Chr",
            if (level == "chromosome-arm") c("chrID.1", "chrID.2"),
            "CN_mode", "mean_max_prob", "n_windows")
  res <- res[want]

  rownames(res) <- NULL
  res
}

#' Merge CN summary with area-based ploidy estimates
#'
#' Merges summarized HMM CN results with area-based ploidy estimates for each sample and chromosome (or arm).
#' Useful for combining HMM and area-based results for reporting or downstream analysis.
#'
#' @param hmm_summarized Data.frame from summarize_cn_mode.
#' @param qploidy_area_ploidy_estimation Object containing area-based ploidy matrices.
#' @param level Character. Merge level: 'chromosome', 'chromosome-arm', or 'sample'.
#'
#' @return Data.frame with merged columns: sample, chromosome, HMM CN mode, area-based CN, confidence metrics, and window count.
#'
#' @importFrom dplyr rename
#' @importFrom tidyr pivot_longer as_tibble
#'
#' @export
merge_cn_summary_with_estimates <- function(hmm_summarized,
                                            qploidy_area_ploidy_estimation,
                                            level = c("chromosome", "chromosome-arm", "sample")) {
  level <- match.arg(level)

  if (!inherits(qploidy_area_ploidy_estimation, "qploidy_area_ploidy_estimation"))
    stop("`qploidy_area_ploidy_estimation` must be a qploidy_area_ploidy_estimation object.")

  ploidy_mat <- qploidy_area_ploidy_estimation$ploidy
  diff_mat   <- qploidy_area_ploidy_estimation$diff_first_second

  if (is.null(rownames(ploidy_mat)) || is.null(rownames(diff_mat)))
    stop("Matrices in `qploidy_area_ploidy_estimation` must have rownames = sample IDs.")

  mat_long <- function(M, value_name) {
    as_tibble(M, rownames = "Sample") |>
      pivot_longer(cols = -Sample, names_to = "Chr", values_to = value_name)
  }

  if (level == "chromosome") {
    est_ploidy_long <- mat_long(ploidy_mat, "CN_area")
    est_diff_long   <- mat_long(diff_mat,   "area_diff_prob")
    est_long <- merge(est_ploidy_long, est_diff_long, by = c("Sample", "Chr"), all = TRUE, sort = FALSE)

    # hmm_summarized is expected to have CN_mode, mean_max_prob, n_windows from summarize_cn_mode()
    out <- merge(hmm_summarized, est_long, by = c("Sample", "Chr"), all.x = TRUE, sort = FALSE)

  } else if (level == "chromosome-arm") {
    if (!all(c("Sample", "chrID.1", "chrID.2") %in% names(hmm_summarized)))
      stop("For level='chromosome-arm', `hmm_summarized` must include: Sample, chrID.1, chrID.2.")
    hmm_summarized$ChrArm <- paste(hmm_summarized$chrID.1, hmm_summarized$chrID.2, sep = ".")

    est_ploidy_long <- mat_long(ploidy_mat, "CN_area") |>
      rename(ChrArm = Chr)
    est_diff_long   <- mat_long(diff_mat,   "area_diff_prob") |>
      rename(ChrArm = Chr)

    est_long <- merge(est_ploidy_long, est_diff_long, by = c("Sample", "ChrArm"), all = TRUE, sort = FALSE)
    out <- merge(hmm_summarized, est_long, by.x = c("Sample", "ChrArm"), by.y = c("Sample", "ChrArm"),
                 all.x = TRUE, sort = FALSE)

    out$chrID.1 <- NULL
    out$chrID.2 <- NULL
    out$ChrArm  <- NULL

  } else { # level == "sample"
    # Aggregate matrices per sample
    CN_area         <- apply(ploidy_mat, 1, function(x) mode(as.numeric(x)))
    area_diff_prob  <- rowMeans(diff_mat, na.rm = TRUE)

    est_sample <- tibble(
      Sample         = names(CN_area),
      CN_area        = as.numeric(CN_area),
      area_diff_prob = as.numeric(area_diff_prob)
    )

    out <- merge(hmm_summarized, est_sample, by = "Sample", all.x = TRUE, sort = FALSE)
    if (!("Chr" %in% names(out))) out$Chr <- NA_character_
    out <- out[c("Sample", "Chr", setdiff(names(out), c("Sample", "Chr")))]
  }

  # ---- Rename to requested schema ----
  # From summarize_cn_mode(): CN_mode -> CN_HMM; mean_max_prob -> HMM_max_prob; n_windows -> HMM_n_windows
  rename_map <- c(CN_mode = "CN_HMM",
                  mean_max_prob = "HMM_max_prob",
                  n_windows = "HMM_n_windows")
  for (old in names(rename_map)) {
    if (old %in% names(out)) names(out)[names(out) == old] <- rename_map[[old]]
  }

  # Ensure all columns exist and ordered as requested
  want <- c("Sample", "Chr", "CN_HMM", "CN_area", "HMM_max_prob", "area_diff_prob", "HMM_n_windows")
  missing <- setdiff(want, names(out))
  for (m in missing) out[[m]] <- NA
  out <- out[want]

  rownames(out) <- NULL
  out
}

##' Initialize monotonic z-score means for HMM ploidy states
##'
##' Computes a monotonic ramp of z-score means (mu) for each ploidy state, centered on the expected ploidy, with padding to avoid boundary collapse.
##' Used for initializing HMM emission parameters in ploidy estimation.
##'
##' @param z Numeric vector. Window-level z-scores.
##' @param cn_grid Integer vector. Copy-number states to consider.
##' @param exp_ploidy Numeric. Expected ploidy value (center of ramp).
##' @param z_range Numeric. Padding added to min/max z for ramp initialization. If NULL, estimated from z interquartile range.
##' @param verbose Logical. If TRUE, prints estimated z_range. Default FALSE.
##' @param z_range_out Logical. If TRUE, expand range (min - z_range, max + z_range). If FALSE, reduce range (min + z_range, max - z_range).
##'
##' @return Named numeric vector of z-score means (mu) for each ploidy state in cn_grid.
##'
##' @details
##' The ramp is constructed as: mu_c = z_mean + step * (c - exp_ploidy), where step is chosen so the lowest and highest ploidy states fit within the padded z range.
##' Ensures monotonic initialization for EM fitting.
##' @keywords internal
##' @export
define_z_limits <- function(z, z_window, cn_grid, exp_ploidy, z_range = NULL, verbose = FALSE, z_range_out = TRUE) {
  state_ids <- as.character(cn_grid)

  if (is.null(z_range) || (length(z_range) == 1 && is.na(z_range))) {
    z_range <- (1/length(cn_grid)) * (as.numeric(quantile(z, probs = 0.75)) - as.numeric(quantile(z, probs = 0.25)))
    vmsg("Estimated z_range from data: %f", verbose = verbose, level = 2, type = ">>", z_range)
  }
  z_mean <- mean(z_window, na.rm = TRUE)
  if (z_range_out) {
    z_lo <- min(z_window, na.rm = TRUE) - z_range
    z_hi <- max(z_window, na.rm = TRUE) + z_range
  } else {
    z_lo <- min(z_window, na.rm = TRUE) + z_range
    z_hi <- max(z_window, na.rm = TRUE) - z_range
  }
  cmin <- min(cn_grid)
  cmax <- max(cn_grid)

  step_lo <- if (exp_ploidy > cmin) (z_mean - z_lo) / (exp_ploidy - cmin) else 0
  step_hi <- if (exp_ploidy < cmax) (z_hi  - z_mean) / (cmax - exp_ploidy) else 0
  step    <- max(step_lo, step_hi, 1e-6)

  mu_vec <- z_mean + step * (as.numeric(cn_grid) - exp_ploidy)
  mu     <- setNames(mu_vec, as.character(cn_grid))
  mu <- mu[state_ids]
  return(mu)
}


#' Update a multi-sample hmm_CN_multi object with a new single-sample hmm_CN result
#'
#' Replaces or adds the results for a given sample in a multi-sample HMM CN object.
#' If the sample already exists, its data and parameters are replaced; otherwise, the new sample is appended.
#'
#' @param hmm_CN_multi An object of class 'hmm_CN_multi' (list with by_window, by_marker, params_samples).
#' @param hmm_CN An object of class 'hmm_CN' (single-sample result with by_window, by_marker, params).
#' @param rm_sample Logical. If TRUE, removes the sample from hmm_CN_multi without adding the new hmm_CN results (useful for cleanup). Default FALSE.
#'
#' @return An updated hmm_CN_multi object with the new or replaced sample's results.
#' @details
#' This function is useful for incrementally building or updating a multi-sample HMM CN result object
#' as new samples are processed. It ensures that only one entry per sample is present in each component.
#'
#' @export
update_hmm_CN_multi <- function(hmm_CN_multi, hmm_CN, rm_sample = FALSE){

  if(!inherits(hmm_CN_multi, "hmm_CN") || !all(c("by_window", "by_marker", "params_samples") %in% names(hmm_CN_multi)))
    stop("hmm_CN multi must be of class 'hmm_CN' with components: by_window, by_marker, params_samples")

  if(!inherits(hmm_CN, "hmm_CN") || !all(c("by_window", "by_marker", "params") %in% names(hmm_CN)))
    stop("hmm_CN must be of class 'hmm_CN' with components: by_window, by_marker, params")

  hmm_CN_multi_new <- hmm_CN_multi

  sample <- unique(hmm_CN$by_window$Sample)
  sample1 <- unique(hmm_CN$by_marker$SampleName)
  if(any(sample != sample1)) stop("Sample in hmm by window and by marker differ")

  if(any(hmm_CN_multi$by_window$Sample %in% sample)){
    idx <- which(hmm_CN_multi$by_window$Sample %in% sample)
    hmm_CN_multi_new$by_window <- hmm_CN_multi$by_window[-idx,]
  }
  if(any(hmm_CN_multi$by_marker$SampleName %in% sample)){
    idx <- which(hmm_CN_multi$by_marker$SampleName %in% sample)
    hmm_CN_multi_new$by_marker <- hmm_CN_multi$by_marker[-idx,]
  }
  if(any(names(hmm_CN_multi_new$params_samples) %in% sample)){
    idx <- which(names(hmm_CN_multi_new$params_samples) %in% sample)
    hmm_CN_multi_new$params_samples <- hmm_CN_multi_new$params_samples[-idx]
  }

  if(rm_sample) return(hmm_CN_multi_new)

  if(length(sample) == 1) {
    params <- list(hmm_CN$params)
    names(params) <- sample
  } else params <- hmm_CN$params

  by_window <- bind_rows(list(hmm_CN_multi_new$by_window, hmm_CN$by_window))
  by_marker <- bind_rows(list(hmm_CN_multi_new$by_marker, hmm_CN$by_marker))

  idx <- which(colnames(by_window) == "post_max")
  idx1 <- order(colnames(by_window)[(idx + 1):ncol(by_window)])
  by_window <- by_window[,c(1:idx, idx+idx1)]

  hmm_CN_multi_new$params_samples <- c(hmm_CN_multi_new$params_samples, params)
  hmm_CN_multi_new$by_window <- by_window
  hmm_CN_multi_new$by_marker <- by_marker

  return(hmm_CN_multi_new)
}
