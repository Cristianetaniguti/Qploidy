#' Numerically stable log-sum-exp
#'
#' Computes \eqn{\log\left(\sum_i e^{x_i}\right)} in a numerically stable way
#' by subtracting the maximum element before exponentiating.
#'
#' @param x A numeric vector.
#'
#' @return A single numeric value: \code{log(sum(exp(x)))} computed stably.
#'
#' @details
#' Directly computing \code{log(sum(exp(x)))} can overflow when \code{x}
#' contains large positive values. This function uses the identity
#' \deqn{\log\sum_i e^{x_i} = m + \log\sum_i e^{x_i - m}, \quad m=\max_i x_i}
#' to avoid overflow/underflow.
#'
#' @examples
#' logsumexp(c(-1000, 0, 1))  # finite result
#'
#' @keywords internal
logsumexp <- function(x) {
  m <- max(x); m + log(sum(exp(x - m)))
}

#' Build a polyploid BAF comb template for a given copy number
#'
#' Constructs a discrete BAF density template on \eqn{[0,1]} for copy number
#' \code{c}, by placing Gaussian kernels at genotype fractions \eqn{d/c},
#' \eqn{d=0,\dots,c}. Reflection kernels keep mass within \eqn{[0,1]}.
#'
#' @param c Integer copy number (e.g., 2 for diploid, 4 for tetraploid).
#' @param M Integer, number of histogram bins over \eqn{[0,1]}.
#' @param bw Numeric, kernel bandwidth (standard deviation) of the Gaussians.
#' @param floor_eps Small positive numeric added to the template before
#'   renormalization to ensure strictly positive probabilities.
#'
#' @return A numeric vector of length \code{M} summing to 1: the BAF template
#'   for copy number \code{c}.
#'
#' @details
#' Peaks are centered at \eqn{d/c}. Binomial weights \code{dbinom(0:c, c, 0.5)}
#' provide a reasonable default relative mass across clusters. Boundary
#' reflections (\eqn{\mu,-\mu,2-\mu}) prevent leakage outside \eqn{[0,1]}.
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
#' Computes a multinomial-style log-likelihood of observed BAF histogram counts
#' against a template (probabilities over the same bins). Empty histograms
#' contribute zero to the likelihood.
#'
#' @param counts Integer vector of length \eqn{M} with BAF bin counts
#'   (e.g., from \code{hist(..., plot = FALSE)$counts}).
#' @param templ Numeric vector of length \eqn{M} with strictly positive
#'   probabilities that sum to 1 (e.g., output of \code{\link{baf_template}}).
#' @param eps Small positive numeric added inside \code{log()} for extra
#'   numerical safety.
#'
#' @return A single numeric value: the log-likelihood \eqn{\sum_i n_i \log p_i}.
#'
#' @details
#' Only bins with \code{counts > 0} are used. If all counts are zero, returns 0.
#' Ensure \code{templ} has no zeros (use \code{floor_eps} in \code{baf_template}).
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
#' Computes the most likely state sequence (Viterbi path) for a Hidden Markov
#' Model given per-position emission log-likelihoods and log transition/prior
#' probabilities.
#'
#' @param ll_em A \eqn{W \times K} numeric matrix of emission log-likelihoods,
#'   where \code{W} is the number of positions/windows and \code{K} is the
#'   number of states. Column order must match \code{logA}.
#' @param logA A \eqn{K \times K} numeric matrix of log transition probabilities
#'   where \code{logA[i, j]} is \eqn{\log p(s_t=j \mid s_{t-1}=i)}. Each row
#'   should represent a valid \emph{log}-probability distribution.
#' @param logpi0 Numeric length-\code{K} vector of initial state log
#'   probabilities at position 1.
#'
#' @return An integer vector of length \code{W} with the 1-based indices of the
#'   most likely states at each position (Viterbi path).
#'
#' @details
#' This implementation operates in log-space for numerical stability. It stores
#' backpointers to reconstruct the optimal path after dynamic programming.
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

#' Plot copy-number segments per window with posterior shading
#'
#' Draws a genome track where each window is plotted as a horizontal segment at
#' its called copy number (CN) and colored by the posterior probability of that
#' call. Panels (facets) are split by chromosome.
#'
#' @param hmm_CN A data frame with one row per window containing at least the
#'   columns:
#'   \itemize{
#'     \item \code{Sample}: sample identifier (character/factor)
#'     \item \code{Chr}: chromosome/contig name
#'     \item \code{Start}, \code{End}: genomic window bounds (numeric or integer)
#'     \item \code{CN_call}: integer copy-number call for the window
#'     \item \code{post_CN*}: posterior probabilities for each CN state, with
#'       column names of the form \code{post_CN2}, \code{post_CN3}, … matching
#'       the CN grid used during HMM inference
#'   }
#' @param sample_id Character scalar. Which sample from \code{hmm_CN$Sample} to
#'   display. Defaults to the first unique value in \code{hmm_CN$Sample}.
#' @param cn_min,cn_max Numeric scalars giving the y-axis limits for CN. By
#'   default they span the minimum/maximum of \code{hmm_CN$CN_call}.
#'
#' @return A \code{ggplot} object. Print it to render, or add layers/scales as
#'   needed.
#'
#' @details
#' The color mapped to each segment is the posterior probability corresponding to
#' the called state in that window. Posterior columns are detected by the prefix
#' \code{"post_CN"} and matched to \code{CN_call} values, so the function is
#' agnostic to the specific CN grid (e.g., 2:6, 2:8).
#'
#' @section Expected columns:
#' The function assumes posterior columns named exactly \code{post_CN<k>} for
#' each copy-number state \code{k}. If your column naming differs, rename them
#' before calling this function.
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
#' plot_cn_track(toy, sample_id = "S1")
#' }
#'
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @export
plot_cn_track <- function(hmm_CN,
                          sample_id = NULL,
                          cn_min = NULL,
                          cn_max = NULL,
                          show_window_lines = FALSE,
                          include_first_in_chr = FALSE,
                          line_color = "grey60",
                          line_alpha = 0.6,
                          line_width = 0.3,
                          line_linetype = "dashed",
                          heights = c(1, 1.2)) {

  stopifnot(inherits(hmm_CN, "hmm_CN"))
  df <- hmm_CN$result

  # defaults BEFORE filtering
  if (is.null(sample_id)) sample_id <- unique(df$Sample)[1]
  x <- dplyr::filter(df, Sample == sample_id)
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
      dplyr::distinct(Chr, Start) |>
      dplyr::arrange(Chr, Start)
    if (!include_first_in_chr) {
      vlines <- vlines |>
        dplyr::group_by(Chr) |>
        dplyr::filter(Start != min(Start)) |>
        dplyr::ungroup()
    }
  }

  # -------- top panel: z dots (fill = w_baf) --------
  p_z <- ggplot(x, aes(Mid, z)) +
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
      legend.position    = "right",
      axis.title.x       = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)
    )

  # -------- bottom panel: CN segments (color = P(CN call)) --------
  p_cn <- ggplot(x) +
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
      legend.position    = "right",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)
    )

  # -------- stack with ggpubr (no patchwork needed) --------
  # align = "v" keeps x-axes aligned; widths are matched automatically
  ggpubr::ggarrange(
    p_z, p_cn,
    ncol = 1, nrow = 2,
    heights = heights,
    align = "v"
  )
}
