if(getRversion() >= "2.15.1") utils::globalVariables(c(
  "BAF", "Density", "Type", "label", "CN"
))

#' Compute BAF Template Likelihoods for Multiple CN States
#'
#' Given a list of BAF vectors (one per window) and a vector of tested copy numbers (cn_grid),
#' this function creates BAF templates for each CN, computes the likelihood of each window's BAF
#' histogram under each template, and returns a matrix of likelihoods (windows × CN states).
#'
#' @param baf_list List of numeric vectors, each containing BAF values for a window.
#' @param cn_grid Integer vector of copy-number states to test (e.g., 2:6).
#' @param M Integer. Number of BAF histogram bins on [0,1]. Default 100.
#' @param bw Numeric. Bandwidth (SD) of the Gaussian kernels used to generate the BAF comb templates. Default 0.03.
#' @param dist Character. Distribution type for template peaks: "gaussian" (default), "beta", "beta-binomial", or "negative binomial".
#' @param plot Logical. If TRUE, return a third item: a list of ggplot objects, one per window, showing observed and template distributions with likelihood and probability text for each CN. Default FALSE.
#'
#' @return A list containing:
#'   \item{ll_matrix}{Matrix of log-likelihoods (n_windows × n_cn), where each entry is the log-likelihood of the window's BAF histogram under the CN template.}
#'   \item{prob_matrix}{Matrix of probabilities (n_windows × n_cn), where each entry is the softmax-transformed probability of the window's BAF histogram under the CN template.}
#'   \item{plots}{List of ggplot objects, one per window, showing observed and template distributions with likelihood and probability text for each CN (if plot = TRUE).}
#' @details Uses the baf_template and baf_ll helpers from Qploidy.
#' @importFrom graphics hist
#' @export
multi_distribution_BAF <- function(baf_list, cn_grid, M = 100, bw = 0.03, plot = FALSE, dist="gaussian") {
  W <- length(baf_list)
  K <- length(cn_grid)
  # Create histogram counts for each window
  breaks <- seq(0, 1, length.out = M + 1)
  hist_counts <- matrix(0L, nrow = W, ncol = M)
  for (i in seq_len(W)) {
    if (length(baf_list[[i]]) == 0 || all(is.na(baf_list[[i]]))) {
      hist_counts[i, ] <- 0
    } else {
      hist_counts[i, ] <- hist(baf_list[[i]], breaks = breaks, plot = FALSE)$counts
    }
  }
  # Create BAF templates for each CN
  templates <- lapply(cn_grid, function(c) baf_template(c, M = M, bw = bw, dist = dist))
  names(templates) <- as.character(cn_grid)
  templates <- lapply(templates, function(t) { t <- pmax(t, 1e-8); t / sum(t) })
  # Compute likelihoods: windows × CN
  ll_matrix <- matrix(NA_real_, nrow = W, ncol = K, dimnames = list(NULL, as.character(cn_grid)))
  for (k in seq_len(K)) {
    templ <- templates[[as.character(cn_grid[k])]]
    ll_matrix[, k] <- apply(hist_counts, 1, baf_ll, templ = templ)
  }
  # Transform log-likelihoods to probabilities using softmax
  prob_matrix <- t(apply(ll_matrix, 1, function(x) {
    exp_x <- exp(x - max(x, na.rm = TRUE)) # numerical stability
    exp_x / sum(exp_x, na.rm = TRUE)
  }))

  plot_list <- NULL
  if (plot) {
    plot_list <- vector("list", W)
    for (i in seq_len(W)) {
      obs_hist <- hist_counts[i, ]
      obs_density <- obs_hist / sum(obs_hist)
      df <- data.frame(
        BAF = seq(0, 1, length.out = M),
        Observed = obs_density
      )
      # Add template densities for each CN
      for (k in seq_len(K)) {
        df[[paste0("CN", cn_grid[k])]] <- templates[[as.character(cn_grid[k])]]
      }
      # Melt for ggplot
      df_long <- tidyr::pivot_longer(df, cols = -BAF, names_to = "Type", values_to = "Density")
      # Add likelihood and probability text for each CN to legend
      legend_labels <- c("Observed", paste0("CN", cn_grid, " (logLik=", sprintf("%.2f", ll_matrix[i, ]), ", P=", sprintf("%.2f", prob_matrix[i, ]), ")"))
      names(legend_labels) <- c("Observed", paste0("CN", cn_grid))
      p <- ggplot(df_long, aes(x = BAF, y = Density, fill = Type)) +
        geom_area(data = subset(df_long, Type != "Observed"), position = "identity", alpha = 0.4) +
        geom_area(data = subset(df_long, Type == "Observed"), position = "identity", alpha = 0.4, fill = "black") +
        geom_line(data = subset(df_long, Type == "Observed"), aes(x = BAF, y = Density), color = "black", size = 1) +
        scale_fill_manual(values = c("Observed" = "black", setNames(RColorBrewer::brewer.pal(min(8, K), "Set1"), paste0("CN", cn_grid))), labels = legend_labels) +
        theme_bw() +
        labs(y = "Density") +
        theme(legend.position = "top")
      plot_list[[i]] <- p
    }
  }
  # Return both matrices in a list, and plot_list if requested
  out <- list(
    ll_matrix = ll_matrix,
    prob_matrix = prob_matrix
  )
  if (!is.null(plot_list)) out$plots <- plot_list
  return(out)
}
