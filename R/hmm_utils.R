# Suppress global variable warnings for non-standard evaluation in dplyr/ggplot2
globalVariables(c(
  "Sample", "Chr", "Start", "End", "CN_call", "prob_call", "w_baf", "Mid"
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

#' Build a polyploid BAF comb template for a given copy number
#'
#' Constructs a discrete BAF density template on [0,1] for a given copy number \code{c}, by placing Gaussian kernels at genotype fractions \eqn{d/c} for \eqn{d=0,...,c}. Reflection kernels keep mass within [0,1]. Binomial weights are used for relative peak heights.
#'
#' @param c Integer. Copy number (e.g., 2 for diploid, 4 for tetraploid).
#' @param M Integer. Number of histogram bins over [0,1]. Default is 101.
#' @param bw Numeric. Kernel bandwidth (standard deviation) for the Gaussians. Default is 0.03.
#' @param floor_eps Numeric. Small positive value added to the template before renormalization to ensure strictly positive probabilities. Default is 1e-8.
#'
#' @return Numeric vector of length \code{M} summing to 1. The BAF template for copy number \code{c}.
#'
#' @details
#' Peaks are centered at \eqn{d/c}. Binomial weights \code{dbinom(0:c, c, 0.5)} provide relative mass across clusters. Boundary reflections (\eqn{\mu,-\mu,2-\mu}) prevent leakage outside [0,1].
#'
#' @examples
#' t4 <- baf_template(4, M = 121, bw = 0.03)
#' sum(t4)            # 1
#' which.max(t4)      # near BAF = 0.5 for c = 4
#'
#' @export
baf_template <- function(c, M=101, bw=0.03, floor_eps=1e-8) {
  x <- seq(0, 1, length.out=M)
  centers <- 0:c / c
  wd <- dbinom(0:c, size=c, prob=0.5)
  # Gaussian kernels with reflection to keep mass in [0,1]
  kfun <- function(mu) dnorm(x, mu, bw) + dnorm(x, -mu, bw) + dnorm(x, 2 - mu, bw)
  dens <- Reduce(`+`, Map(function(mu,w) w * kfun(mu), centers, wd))
  dens <- dens + floor_eps              # strictly positive
  dens <- dens / sum(dens)              # renormalize
  dens
}

#' BAF histogram log-likelihood under a template
#'
#' Computes the multinomial-style log-likelihood of observed BAF histogram counts against a template (probabilities over the same bins). Only bins with nonzero counts contribute. If all counts are zero, returns 0.
#'
#' @param counts Integer vector of length \code{M}. BAF bin counts (e.g., from \code{hist(..., plot = FALSE)$counts}).
#' @param templ Numeric vector of length \code{M}. Strictly positive probabilities summing to 1 (e.g., output of \code{baf_template}).
#' @param eps Numeric. Small positive value added inside \code{log()} for extra numerical safety. Default is 1e-8.
#'
#' @return Numeric scalar. The log-likelihood \eqn{\sum_i n_i \log p_i}.
#'
#' @details
#' Only bins with \code{counts > 0} are used. If all counts are zero, returns 0. Ensure \code{templ} has no zeros (use \code{floor_eps} in \code{baf_template}).
#'
#' @examples
#' templ <- baf_template(4, M = 21)
#' counts <- rmultinom(1, size = 200, prob = templ)[,1]
#' baf_ll(counts, templ)
#'
#' @export
baf_ll <- function(counts, templ, eps=1e-8) {
  idx <- counts > 0L
  if (!any(idx)) return(0)  # empty histogram contributes nothing
  sum(counts[idx] * log(templ[idx] + eps))
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

#' Plot copy-number segments per window with posterior shading and BAF
#'
#' Plots a genome track for each window, showing the called copy number (CN) as a horizontal segment and coloring by the posterior probability of that call. Facets are split by chromosome. The top panel shows z-scores per window, colored by BAF weight; the bottom panel shows CN segments colored by posterior probability. The new top panel shows BAF values per SNP, colored by BAF weight, for the selected sample and chromosomes.
#'
#' @param hmm_CN An object of class \code{hmm_CN} (output from \code{hmm_estimate_CN}), containing a \code{result} data frame with one row per window and columns: \code{Sample}, \code{Chr}, \code{Start}, \code{End}, \code{CN_call}, \code{post_CN*}, etc.
#' @param qploidy_standarize_result An object of class \code{qploidy_standardization} (output from \code{standardize}), used to extract BAF values for the BAF panel.
#' @param sample_id Character scalar. Which sample from \code{hmm_CN$Sample} to display. Defaults to the first unique value in \code{hmm_CN$Sample}.
#' @param cn_min,cn_max Numeric scalars. Y-axis limits for CN. Defaults span the min/max of \code{CN_call}.
#' @param show_window_lines Logical. If TRUE, show dashed vertical lines at window boundaries.
#' @param include_first_in_chr Logical. If TRUE, include the first window line in each chromosome.
#' @param line_color,line_alpha,line_width,line_linetype Appearance settings for window boundary lines.
#' @param heights Numeric vector of length 2 or 3. Relative heights of the BAF, z, and CN panels.
#'
#' @return A \code{ggplot} object (from ggpubr::ggarrange). Print to render, or add layers/scales as needed.
#'
#' @details
#' Posterior columns are detected by the prefix "post_CN" and matched to \code{CN_call} values, so the function is agnostic to the specific CN grid. The BAF panel shows per-SNP BAF values for the selected sample and chromosomes, colored by the BAF weight for the corresponding region. The z panel shows window z-scores, colored by BAF weight. The CN panel shows copy number segments colored by posterior probability.
#'
#' @section Expected columns:
#' The function assumes posterior columns named exactly \code{post_CN<k>} for each copy-number state \code{k}. If your column naming differs, rename them before calling this function.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(ggplot2)
#' toy <- data.frame(
#'   Sample   = "S1",
#'   Chr      = c("chr1","chr1","chr1"),
#'   Start    = c(1, 1e6, 2e6),
#'   End      = c(1e6-1, 2e6-1, 3e6-1),
#'   CN_call  = c(2,3,2),
#'   post_CN2 = c(0.95, 0.05, 0.9),
#'   post_CN3 = c(0.05, 0.94, 0.1)
#' )
#' # qploidy_standarize_result should be a qploidy_standardization object with $data containing BAF values
#' plot_cn_track(toy, qploidy_standarize_result, sample_id = "S1")
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @importFrom magrittr %>%
#' @importFrom ggpubr ggarrange
#'
#' @export
plot_cn_track <- function(hmm_CN,
                          qploidy_standarize_result,  # for BAF plot
                          sample_id = NULL,
                          cn_min = NULL,
                          cn_max = NULL,
                          show_window_lines = FALSE,
                          include_first_in_chr = FALSE,
                          line_color = "grey60",
                          line_alpha = 0.6,
                          line_width = 0.3,
                          line_linetype = "dashed",
                          heights = c(2, 2.5)) {

  stopifnot(inherits(hmm_CN, "hmm_CN"))
  stopifnot(inherits(qploidy_standarize_result, "qploidy_standardization"))

  df <- hmm_CN$result

  # defaults BEFORE filtering
  if (is.null(sample_id)) sample_id <- unique(df$Sample)[1]
  x <- filter(df, Sample == sample_id)
  if (nrow(x) == 0) stop("No rows for sample_id = ", sample_id)

  if (is.null(cn_min)) cn_min <- min(x$CN_call, na.rm = TRUE)
  if (is.null(cn_max)) cn_max <- max(x$CN_call, na.rm = TRUE)

  # compute P(CN_call) per window
  cn_states <- as.integer(gsub("^post_CN", "", grep("^post_CN", names(x), value = TRUE)))
  post_mat  <- as.matrix(x[, paste0("post_CN", cn_states), drop = FALSE])
  idx       <- match(x$CN_call, cn_states)
  x$prob_call <- post_mat[cbind(seq_len(nrow(x)), idx)]

  # coordinates
  x$Start <- as.numeric(x$Start)
  x$End   <- as.numeric(x$End)
  x$Mid   <- (x$Start + x$End)/2

  # dashed verticals (one per unique Start, per Chr)
  vlines <- NULL
  if (show_window_lines) {
    vlines <- x |>
      distinct(Chr, Start) |>
      arrange(Chr, Start)
    if (!include_first_in_chr) {
      vlines <- vlines |>
        group_by(Chr) |>
        filter(Start != min(Start)) |>
        ungroup()
    }
  }

  data_sample2 <- qploidy_standarize_result$data %>% filter(SampleName == sample_id & Chr %in% unique(x$Chr))

  # Per-chromosome genomic range (based on BAF; you can also combine with Start/End if you prefer)
  chrom_ghost <- data_sample2 %>%
    group_by(Chr) %>%
    summarise(Position = range(Position, na.rm = TRUE), .groups = "drop")
  # This gives 2 rows per Chr: min and max Position

  # -------- top panel: z dots (fill = w_baf) --------
  p_z <- ggplot(x, aes(Mid, z)) +
    # force x scale to span full chromosome range
    geom_blank(
      data = chrom_ghost,
      inherit.aes = FALSE,
      aes(x = Position, y = 0)
    ) +
    # connect points within each chromosome
    geom_line(aes(group = Chr, alpha = w_baf), color = "grey50", linewidth = 0.6) +
    scale_alpha(range = c(0.2, 0.9), guide = "none") +
    # points colored by w_baf
    geom_point(aes(color = w_baf), size = 1.8, alpha = 0.9) +
    { if (show_window_lines)
      geom_vline(data = vlines, aes(xintercept = Start),
                 color = line_color, alpha = line_alpha,
                 linewidth = line_width, linetype = line_linetype)
      else NULL } +
    facet_wrap(~ Chr, scales = "free_x", nrow = 1) +
    scale_color_distiller(palette = "RdBu", direction = -1, limits = c(0, 1), name = "BAF weight")+
    labs(x = NULL, y = "z", title = sample_id) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      strip.background   = element_rect(fill = "grey95"),
      legend.position    = "none",
      axis.title.x       = element_blank(),
      axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)
    )

  # -------- bottom panel: CN segments (color = P(CN call)) --------
  p_cn <- ggplot(x) +
    geom_blank(
      data = chrom_ghost,
      inherit.aes = FALSE,
      aes(x = Position, y = 0)
    ) +
    geom_segment(aes(x = Start, xend = End, y = CN_call, yend = CN_call, color = prob_call),
                 linewidth = 2, lineend = "butt") +
    { if (show_window_lines)
      geom_vline(data = vlines, aes(xintercept = Start),
                 color = line_color, alpha = line_alpha,
                 linewidth = line_width, linetype = line_linetype) else NULL } +
    facet_wrap(~ Chr, scales = "free_x", nrow = 1) +
    scale_y_continuous(breaks = seq(cn_min, cn_max, by = 1), minor_breaks = NULL) +
    scale_color_viridis_c(name = "P(CN call)", limits = c(0, 1)) +
    labs(x = "Genomic position (bp)", y = "Copy number") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      strip.background   = element_rect(fill = "grey95"),
      legend.position    = "bottom",
      axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)
    )

  # -------- BAF -----

  data_tagged <- data_sample2 %>%
    group_by(Chr) %>%
    group_modify(function(df, key) {
      chr_starts <- vlines %>%
        filter(Chr == key$Chr[1]) %>%
        arrange(Start) %>%
        pull(Start)
      # regions: [-Inf, s1), [s1, s2), ..., [sN, +Inf)
      brks <- c(-Inf, chr_starts, Inf)
      df$region_id <- cut(df$Position, breaks = brks, labels = FALSE, right = FALSE)
      df
    }) %>%
    ungroup()

  ## 2) Build a region-to-w_baf table
  regions_tbl <- vlines %>%
    arrange(Chr, Start) %>%
    count(Chr, name = "n_vlines") %>%
    mutate(n_regions = n_vlines + 1L) %>%
    select(Chr, n_regions) %>%
    rowwise() %>%
    mutate(region_id = list(seq_len(n_regions))) %>%
    unnest(region_id) %>%
    ungroup()

  # If your x$w_baf is for a single chromosome (like chr1) and matches the
  # region order, attach it like this:
  # (If you have per-Chr w_baf values already in a data.frame, join that instead.)
  w_baf_tbl <- regions_tbl %>%
    group_by(Chr) %>%
    mutate(w_baf = x$w_baf[region_id]) %>%  # ensure the vector aligns: length == n_regions
    ungroup()

  ## 3) Join w_baf onto points and plot
  plot_df <- data_tagged %>%
    left_join(w_baf_tbl, by = c("Chr", "region_id"))

  p_baf <- ggplot(plot_df, aes(x = Position, y = baf, color = w_baf)) +
    geom_blank(
      data = chrom_ghost,
      inherit.aes = FALSE,
      aes(x = Position, y = 0)
    ) +
    geom_point(alpha = 0.7, size = 1) +
    { if (show_window_lines)
      geom_vline(data = vlines, aes(xintercept = Start),
                 color = line_color, alpha = line_alpha,
                 linewidth = line_width, linetype = line_linetype) else NULL } +
    facet_wrap(~Chr, scales = "free_x", nrow=1) +
    theme_bw() +
    ylab("BAF") +
    scale_color_distiller(palette = "RdBu", direction = -1, limits = c(0, 1), name = "BAF weight")+
    theme(
      axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
      legend.position = "top",
      text = element_text(size = 12),
      axis.title.x       = element_blank()
    )

  # -------- stack with ggpubr (no patchwork needed) --------
  # align = "v" keeps x-axes aligned; widths are matched automatically
  ggarrange(
    p_baf, p_z, p_cn,
    ncol = 1, nrow = 3,
    heights = heights,
    align = "v"
  )
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
worker <- function(sid, obj, dots) {
  tryCatch({
    res <- do.call(hmm_estimate_CN,
                   c(list(obj, sample_id = sid), dots))
    res
  }, error = function(e) {
    warning(sprintf("Sample '%s' failed: %s", sid, conditionMessage(e)))
    NULL
  })
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

  # unwrap if hmm_CN-like object with $result
  dat <- if (is.data.frame(df)) df else if (!is.null(df$result) && is.data.frame(df$result)) df$result else
    stop("Input must be a data.frame or an hmm_CN-like object with a data.frame at `$result`.")

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
