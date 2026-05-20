#' Polyploid copy-number HMM (z + BAF) per sample and window
#'
#' Aggregates SNP-level data into windows, builds BAF histograms, and fits a
#' Hidden Markov Model that combines a Gaussian emission for depth/z and a
#' polyploid BAF “comb” emission to infer copy number (CN) per window.
#' Returns a tidy data frame with CN calls and posterior probabilities that is
#' ready for plotting (e.g., with \code{\link{plot_cn_track}}) or merging across
#' chromosomes/batches.
#'
#' @param qploidy_standarize_result An object of class \code{qploidy_standardization} as returned by \code{standardize()}, or \code{NULL}. When \code{NULL}, raw data must be supplied via \code{data} and \code{geno.pos}.
#' @param data Optional. A \code{data.frame} with columns \code{MarkerName}, \code{SampleName}, \code{X}, \code{Y}, \code{R}, and \code{ratio}. Required when \code{qploidy_standarize_result} is \code{NULL}.
#' @param geno.pos Optional. A \code{data.frame} with columns \code{MarkerName}, \code{Chromosome}, and \code{Position}. Required when \code{qploidy_standarize_result} is \code{NULL}.
#' @param use_values Character vector of length 2 specifying which columns to use as the BAF-like and depth-like signals. When \code{qploidy_standarize_result} is provided, the first element can be \code{"BAF"} (uses column \code{baf}) or \code{"ratio"} (uses column \code{ratio}), and the second element can be \code{"zscore"} (uses column \code{z}) or \code{"R"} (uses column \code{R}). All four combinations are accepted: \code{c("BAF", "zscore")} (default), \code{c("BAF", "R")}, \code{c("ratio", "zscore")}, \code{c("ratio", "R")}. When \code{qploidy_standarize_result} is \code{NULL}, the only accepted value is \code{c("ratio", "R")}.
#' @param sample_id Character scalar. Sample identifier to analyze; rows are filtered as \code{SampleName == sample_id}.
#' @param chr Optional. Character or integer vector specifying chromosomes to include. If \code{NULL}, all chromosomes are used.
#' @param segment_zscore Logical. If \code{TRUE}, segment z-scores using changepoint detection to define windows instead of fixed SNP counts. Default \code{TRUE}.
#' @param snps_per_window Integer. Number of SNPs per window when windows are created. Default \code{500}.
#' @param reflect Logical. If TRUE (default), apply reflection for continuous
#'   kernels to keep mass within [0,1] (useful for gaussian).
#' @param add_uniform Logical. If TRUE, add a uniform noise component to the
#'   mixture before renormalization. Default FALSE.
#' @param min_snps_per_window Integer. Minimum SNPs required to keep a window. If NULL, a dynamic value is chosen based on chromosome size (see code for details). Windows with fewer SNPs are dropped. Default: NULL (dynamic), or user-specified value.
#' @param cn_grid Integer vector of copy-number states to consider (e.g., \code{1:4}). Unlikely values will be discarded during estimation if their z means are not monotonic with ploidy.
#' @param M Integer. Number of BAF histogram bins on [0,1]. Default \code{121}.
#' @param max_iter Integer. Maximum EM iterations. Default \code{60}.
#' @param het_quantile Numeric. Quantile used to scale BAF emission weight based on heterozygote count. Default \code{0.8}.
#' @param baf_weight Numeric. Overall weight applied to BAF emission (0–1). Default \code{1}.
#' @param z_range Numeric. Padding added to min/max z for initial mean estimation. Default \code{0.2}.
#' @param transition_jump Numeric. Diagonal value for transition matrix (probability to stay in same CN state). Default \code{0.995}.
#' @param initial_prob Numeric. Initial probability for the best CN state in the initial state distribution (pi0). Default \code{0.15}. Sets the prior probability for the expected ploidy (or best CN from BAF model) at the first window; remaining probability is distributed uniformly across other states. If the best CN is not found, pi0 is uniform across all states.
#' @param z_only Logical. If \code{TRUE}, fit the HMM using the z-emission only (ignores BAF). Default \code{FALSE}.
#' @param verbose Logical. If \code{TRUE}, print progress messages. Default \code{TRUE}.
#' @param exp_ploidy Numeric. Expected ploidy value. If \code{NA} or \code{NULL}, it is set to the best CN from the BAF model. Default \code{NA}.
#' @param rm_outliers Logical. If \code{TRUE}, remove outliers from z-scores before HMM fitting. Default \code{TRUE}.
#' @param outlier_alpha Numeric. Alpha threshold for outlier removal in z-scores. Default 0.05.
#' @param selected_model Optional. An object of class \code{selected_BAF_model}. If provided, it is used instead of running \code{select_best_baf_model}.
#' @param dists Character vector. Distribution families to test in BAF model selection (e.g., c("gaussian", "beta", ...)).
#' @param bw_grid Numeric vector. Bandwidth (SD or concentration) values to try for kernel smoothing in BAF model selection. Default: c(0.02, 0.03, 0.04).
#' @param add_uniform_grid Logical vector. Whether to include a uniform noise component in the BAF template. Default: c(FALSE, TRUE).
#' @param uniform_weight_grid Numeric vector. Mixture weights for the uniform component (when enabled). Default: c(0.01, 0.03, 0.05, 0.10, 0.15).
#' @param param_count Optional named integer vector. Number of free parameters per distribution (for BIC penalty in BAF model selection). If NULL, defaults to 0 for all.
#' @param count_grid_as_params Logical. If TRUE (default), adds +1 to BIC penalty for each hyperparameter tuned by grid search (bw, and uniform_weight if used).
#' @param correct_scale Logical. If TRUE (default), the BAF log-likelihood is corrected by the number of markers with valid BAF values in each window, so that windows with different numbers of markers contribute equally to the combined emission. Prevents windows with many markers from dominating the HMM via the BAF term alone.
#' @param min_het_frac Numeric in [0,1]. Threshold for the fraction of BAF values in \code{het_range} considered heterozygous. If the observed heterozygous fraction exceeds this value, CN=1 is excluded from \code{cn_grid} during BAF model selection and per-window likelihood computation, as a meaningful proportion of heterozygous loci makes haploid (CN=1) implausible. Default \code{0.05}.
#' @param het_range Numeric vector of length 2. BAF interval used to define heterozygous loci (default \code{c(0.2, 0.8)}). Used by the \code{min_het_frac} filter and by \code{plot_heterozygosity}.
#' @param dosage_threshold Numeric in [0,1]. Minimum posterior probability required for a dosage call to be counted as a heterozygote when computing the BAF emission weight (\code{w_baf}) per window. Markers whose maximum dosage posterior falls below this threshold are excluded from the heterozygote count. Default \code{0.6}.
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
#' head(res$by_window)
#' res$params
#' }
#'
#' @seealso \code{\link{generate_baf_template}}, \code{\link{baf_log_likelihood}},
#'   \code{\link{viterbi}}, \code{\link{plot_cn_track}}
#' @importFrom stats dnorm quantile sd
#' @importFrom graphics hist
#' @importFrom utils tail
#' @importFrom dplyr filter
#'
#' @export
hmm_estimate_CN <- function(
    qploidy_standarize_result = NULL,
    sample_id,
    data = NULL,
    geno.pos = NULL,
    use_values = c("BAF", "zscore"),
    chr = NULL,
    reflect = TRUE,
    add_uniform = FALSE,
    segment_zscore = TRUE,
    snps_per_window = 50,
    min_snps_per_window = NULL,
    cn_grid = 1:4,
    M = 100,
    max_iter = 60,
    het_quantile = 0.8, # increase this value to reduce the weight of baf when few hets
    baf_weight = 0.5,
    z_range = NULL,
    transition_jump = 0.995, # decrease this value if you think there changes in CN is likely
    initial_prob = 0.95, # Initial probability for the best CN state in the initial state distribution (pi0). Default 0.15. Sets the prior probability for the expected ploidy (or best CN from BAF model) at the first window; remaining probability is distributed uniformly across other states. If the best CN is not found, pi0 is uniform across all states.
    z_only = FALSE,
    verbose = TRUE,
    exp_ploidy = NA,
    rm_outliers = TRUE,
    outlier_alpha = 0.05,
    selected_model = NULL,
    dists = c("gaussian", "beta", "beta_binomial", "negative_binomial"),
    bw_grid = c(0.02, 0.03, 0.04),
    add_uniform_grid = c(FALSE, TRUE),
    uniform_weight_grid = c(0.01, 0.03, 0.05, 0.10, 0.15),
    param_count = NULL,
    count_grid_as_params = TRUE,
    correct_scale = TRUE,
    min_het_frac = 0.05,
    het_range = c(0.2,0.8),
    dosage_threshold = 0.6
) {

  # --- input checks ---
  vmsg("Preparing inputs and applying initial filters", verbose = verbose, level = 0, type = ">>")

  if (!is.character(sample_id) || length(sample_id) != 1 || nchar(sample_id) == 0) {
    stop("sample_id must be a non-empty character scalar.")
  }

  # Validate and resolve use_values
  valid_uv_std <- list(c("BAF", "zscore"), c("BAF", "R"), c("ratio", "zscore"), c("ratio", "R"))
  valid_uv_raw <- c("ratio", "R")
  if (!is.null(qploidy_standarize_result)) {
    if (!is(qploidy_standarize_result, "qploidy_standardization")) {
      stop("Input must be a qploidy_standardization object as returned by standardize().")
    }
    if (!any(sapply(valid_uv_std, identical, use_values))) {
      stop("When qploidy_standarize_result is provided, use_values must be one of: c(\"BAF\", \"zscore\"), c(\"BAF\", \"R\"), c(\"ratio\", \"zscore\"), c(\"ratio\", \"R\").")
    }
    baf_col <- if (use_values[1] == "BAF") "baf" else "ratio"
    z_col   <- if (use_values[2] == "zscore") "z" else "R"
    df <- as.data.frame(qploidy_standarize_result$data)
    if (nrow(df) == 0) stop("Input data is empty. No SNPs found for any sample.")
    d <- filter(df, SampleName == sample_id)
    if (nrow(d) == 0) stop(sprintf("No data found for sample '%s'.", sample_id))
  } else {
    if (!identical(use_values, valid_uv_raw)) {
      stop("When qploidy_standarize_result is NULL, use_values must be c(\"ratio\", \"R\").")
    }
    if (is.null(data) || is.null(geno.pos)) {
      stop("When qploidy_standarize_result is NULL, both 'data' and 'geno.pos' must be provided.")
    }
    req_data <- c("MarkerName", "SampleName", "ratio", "R")
    req_gp   <- c("MarkerName", "Chromosome", "Position")
    miss_d  <- setdiff(req_data, names(data))
    miss_gp <- setdiff(req_gp,   names(geno.pos))
    if (length(miss_d))  stop(paste("'data' is missing columns:",   paste(miss_d,  collapse = ", ")))
    if (length(miss_gp)) stop(paste("'geno.pos' is missing columns:", paste(miss_gp, collapse = ", ")))
    baf_col <- "ratio"
    z_col   <- "R"
    d <- data[data$SampleName == sample_id, , drop = FALSE]
    if (nrow(d) == 0) stop(sprintf("No data found for sample '%s'.", sample_id))
    gp <- geno.pos[, req_gp]
    colnames(gp)[colnames(gp) == "Chromosome"] <- "Chr"
    d <- merge(d, gp, by = "MarkerName", all.x = FALSE)
    if (nrow(d) == 0) stop("No markers remain after merging data with geno.pos.")
  }

  # subset to chromosomes
  if (!is.null(chr)) {
    chrs <- if (is.numeric(chr)) unique(d$Chr)[chr] else chr
    d <- filter(d, Chr %in% chrs)
    if (nrow(d) == 0) stop("No data found for specified chromosomes after filtering.")
  }

  # Replace R == 0 with NA (zero depth is uninformative and distorts z/R distributions)
  if ("R" %in% names(d)) d$R[d$R == 0] <- NA

  # Remove outliers from depth/z scores
  if (rm_outliers) {
    n_na <- sum(is.na(d[[z_col]]))
    if (z_col != "z") names(d)[names(d) == z_col] <- "z"
    d <- rm_outlier(d, z = TRUE, alpha = outlier_alpha)
    if (z_col != "z") names(d)[names(d) == "z"] <- z_col
    vmsg("Outliers removed from %s scores: %s", verbose = verbose, level = 1, type = ">>", z_col, sum(is.na(d[[z_col]])) - n_na)
  }

  # Remove markers with missing BAF/ratio and z/R
  rm <- which(is.na(d[[z_col]]))
  if (length(rm) > 0) d <- d[-rm, , drop = FALSE]

  vmsg("Inputs good to go!", verbose = verbose, level = 1, type = ">>")

  # --- set expected ploidy ---
  # Calculate expected ploidy using sample-level BAF distribution
  vmsg("Setting internal parameters", verbose = verbose, level = 0, type = ">>")

  if (!is.null(selected_model)) {
    if (!inherits(selected_model, "selected_BAF_model")) {
      stop("selected_model argument must be of class 'selected_BAF_model'.")
    }
    # Use user-provided selected_model
    vmsg("Using user-provided selected_model object", verbose = verbose, level = 1, type = ">>")
  } else {
    selected_model <- select_best_baf_model(baf_vec = d[[baf_col]],
                                            sample = sample_id,
                                            cn_grid= cn_grid,
                                            dists = dists,
                                            M = M,
                                            reflect = reflect,
                                            bw_grid = bw_grid,
                                            add_uniform_grid = add_uniform_grid,
                                            uniform_weight_grid = uniform_weight_grid,
                                            param_count = param_count,
                                            count_grid_as_params = count_grid_as_params,
                                            min_het_frac = min_het_frac,
                                            het_range = het_range)
  }

  # Use exp_ploidy argument if provided, otherwise use selected_model$best$best_cn
  if (is.null(exp_ploidy) || (length(exp_ploidy) == 1 && is.na(exp_ploidy))) {
    exp_ploidy <- selected_model$best$best_cn
    vmsg("exp_ploidy not provided, using best CN from BAF model: %s", verbose = verbose, level = 1, type = ">>", exp_ploidy)
  } else {
    vmsg("exp_ploidy provided by user: %s", verbose = verbose, level = 1, type = ">>", exp_ploidy)
  }

  vmsg("Parameters ready", verbose = verbose, level = 1, type = ">>")
  # --- build windows ---
  vmsg("Defining windows", verbose = verbose, level = 0, type = ">>")
  d[["Position"]] <- as.numeric(d[["Position"]])
  d <- d[order(d[["Chr"]], d[["Position"]]), ]

  if(is.null(min_snps_per_window)){
    # If min_snp_per_window is not defined, the minimum will be set based on the smaller chromosome number
    floor <- 5 # hard default
    frac <- 0.15 # hard default
    mk_by_chrom <- d %>% group_by(Chr) %>% summarize(n = n())
    x <- min(mk_by_chrom$n)
    m <- max(floor, floor(x * frac))
    min_snps_per_window <- min(m, floor(x / 2))
  }

  # Segmented z-score
  if (segment_zscore){
    vmsg("Using z-scores changepoint detection to define windows", verbose = verbose, level = 1, type = ">>")
    if (z_col != "z") names(d)[names(d) == z_col] <- "z"
    d <- add_changepoint_windows(dat = d, minseglen = min_snps_per_window)
    if (z_col != "z") names(d)[names(d) == "z"] <- z_col

  } else {
    # simple fixed-size windows
    vmsg("Using user-defined fixed-size intervals to define windows", verbose = verbose, level = 1, type = ">>")
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
  vmsg("Summarizing data by window", verbose = verbose, level = 1, type = ">>")

  win_df <- do.call(rbind, by(d, list(d[["Chr"]], d[[win_col]]), function(x) {
    if (is.null(x) || nrow(x) == 0) return(NULL)
    data.frame(
      Chr        = x[1, "Chr"],
      WindowID   = x[1, win_col],
      Start      = min(x[["Position"]]),
      End        = max(x[["Position"]]),
      n_snps     = nrow(x),
      z_mean     = mean(x[[z_col]], na.rm = TRUE)
    )
  }))
  rownames(win_df) <- NULL
  if(any(is.nan(win_df$z_mean))) {
    warning("Some windows have no z-score values.")
  }
  if (is.null(win_df) || nrow(win_df) == 0) stop("No windows could be formed. Check input data and window parameters.")

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
    baf_list[[i]] <- d[[baf_col]][wrows]
    if (length(baf_list[[i]]) == 0 || all(is.na(baf_list[[i]]))) {
      warning(sprintf("Window %d (Chr %s, WindowID %s) has no valid BAF values.",
                      i, win_df$Chr[i], win_df$WindowID[i]))
    }
  }
  z <- win_df$z_mean
  if (length(z) == 0 || all(is.na(z))) stop("No valid z values found in any window.")

  vmsg("Windows ready", verbose = verbose, level = 1, type = ">>")

  keep <- win_df$n_snps >= min_snps_per_window
  if (!any(keep)) stop("All windows dropped by min_snps_per_window filter. Try lowering min_snps_per_window or check input data.")
  # n_het is now filled later using vectorized dosages
  # Handle single-window case: assign CN by BAF likelihood only, skip HMM/EM
  if (sum(keep) == 1) {
    vmsg("Only one window remains after filtering. Assigning CN by BAF likelihood only", verbose = verbose, level = 1, type = ">>")
    # Use BAF likelihoods to assign CN
    ll_baf_matrix <- compute_baf_likelihoods(baf_list[[1]],
                                             cn_grid,
                                             M = M,
                                             bw = selected_model$best$bw,
                                             plot = FALSE,
                                             dist = selected_model$best$dist,
                                             reflect = reflect,
                                             add_uniform = selected_model$best$add_uniform,
                                             uniform_weight = selected_model$best$uniform_weight,
                                             min_het_frac = min_het_frac,
                                             het_range = het_range)

    cn_call <- cn_grid[which.max(ll_baf_matrix$ll_vec)]
    post_max <- rep(1, 1)
    post_df <- as.data.frame(matrix(0, nrow=1, ncol=length(cn_grid)))
    names(post_df) <- paste0("post_CN", cn_grid)
    post_df[1, which.max(ll_baf_matrix$prob_vec)] <- 1
    # n_het is now calculated using dosages
    dosages <- mapply(function(x, y) call_BAF_dosages(x,
                                                      cn = y,
                                                      bw = selected_model$best$bw,
                                                      plot = FALSE,
                                                      dist = selected_model$best$dist,
                                                      add_uniform = selected_model$best$add_uniform,
                                                      uniform_weight = selected_model$best$uniform_weight), baf_list[keep], cn_call, SIMPLIFY = FALSE)
    n_het <- sum(dosages[[1]]$dosage != 0 & dosages[[1]]$dosage != cn_call[1], na.rm = TRUE)
    result <- cbind(
      data.frame(
        Sample    = sample_id,
        Chr       = win_df$Chr[keep],
        WindowID  = win_df$WindowID[keep],
        Start     = win_df$Start[keep],
        End       = win_df$End[keep],
        n_snps    = win_df$n_snps[keep],
        n_het     = n_het,
        z         = win_df$z_mean[keep],
        w_baf     = 1,
        CN_call   = cn_call,
        post_max  = post_max,
        stringsAsFactors = FALSE
      ),
      post_df
    )
    params <- list(
      cn_grid = cn_grid,
      distribution = selected_model$best$dist,
      mu = NA,
      sigma = NA,
      A = NA,
      pi0 = NA,
      bins = M,
      bw = selected_model$best$bw,
      loglik = NA,
      z_range = z_range,
      het_quantile = het_quantile,
      baf_weight = baf_weight,
      transition_jump = transition_jump,
      z_only = z_only,
      exp_ploidy = exp_ploidy,
      rm_outliers = rm_outliers,
      outlier_alpha = outlier_alpha,
      segment_zscore = segment_zscore,
      snps_per_window = snps_per_window,
      min_snps_per_window = min_snps_per_window,
      add_uniform = add_uniform,
      uniform_weight = selected_model$best$uniform_weight
    )
    window_map <- match(d$.__w__, result$WindowID)
    d$w_baf    <- result$w_baf[window_map]
    d$CN_call  <- result$CN_call[window_map]
    d$post_max <- result$post_max[window_map]
    vmsg("Done!", verbose = verbose, level = 1, type = ">>")
    return(structure(list(by_window = result, by_marker = d, params = params), class = "hmm_CN"))
  }

  # Generate BAF likelihoods and probabilities per window
  # Uses parameters from selected_model
  if(!z_only){
    vmsg("Generating BAF likelihoods by window", verbose = verbose, level = 0, type = ">>")
    grid1 <- unique(c(1,cn_grid)) # always test 1 for LOH loci

    baf_results <- lapply(baf_list, function(baf_vec) compute_baf_likelihoods(baf_vec,
                                                                              grid1,
                                                                              M = M,
                                                                              bw = selected_model$best$bw,
                                                                              plot = FALSE,
                                                                              dist = selected_model$best$dist,
                                                                              reflect = reflect,
                                                                              add_uniform = selected_model$best$add_uniform,
                                                                              uniform_weight = selected_model$best$uniform_weight,
                                                                              het_range = het_range,
                                                                              min_het_frac = min_het_frac))

    ll_baf_matrix <- do.call(rbind, lapply(baf_results, function(res) res$ll_vec))
    ploidies_temp <- apply(ll_baf_matrix, 1, which.max)
    ploidies_temp <- grid1[ploidies_temp]

    colnames(ll_baf_matrix) <- paste0("CN", grid1)
    if(!any(cn_grid == 1)) {
      ll_baf_matrix <- ll_baf_matrix[, -which(colnames(ll_baf_matrix) == "CN1"), drop=FALSE]
    }
    vmsg("BAF likelihoods generated", verbose = verbose, level = 1, type = ">>")

    vmsg("Generating BAF weights by window", verbose = verbose, level = 0, type = ">>")

    vmsg("Calling dosages", verbose = verbose, level = 1, type = ">>")
    dosages <- mapply(function(x, y) call_BAF_dosages(x,
                                                      cn = y,
                                                      bw = selected_model$best$bw,
                                                      plot = FALSE,
                                                      dist = selected_model$best$dist,
                                                      add_uniform = selected_model$best$add_uniform,
                                                      uniform_weight = selected_model$best$uniform_weight),
                      baf_list, ploidies_temp, SIMPLIFY = FALSE)

    # For each window, count heterozygotes: dosage != 0 & dosage != ploidies_temp & dosage_prob > dosage_threshold
    vmsg("Counting heterozygotes per window", verbose = verbose, level = 1, type = ">>")

    n_het_before_thresh <- vapply(seq_along(dosages), function(i) {
      sum(dosages[[i]]$data$dosage != 0 &
            dosages[[i]]$data$dosage != ploidies_temp[i], na.rm = TRUE)
    }, integer(1))

    n_het_window <- vapply(seq_along(dosages), function(i) {
      sum(dosages[[i]]$data$dosage != 0 &
            dosages[[i]]$data$dosage != ploidies_temp[i] &
            dosages[[i]]$data$max_prob > dosage_threshold, na.rm = TRUE)
    }, integer(1))

    # Warn about windows where all heterozygotes were discarded by the dosage_threshold filter
    prob_filtered_idx <- which(n_het_before_thresh > 0 & n_het_window == 0)
    if (length(prob_filtered_idx) > 0) {
      prob_filtered_labels <- paste0(
        win_df$Chr[prob_filtered_idx], ":", win_df$WindowID[prob_filtered_idx]
      )
      vmsg(
        paste0(
          "Warning: %d window(s) had all heterozygous markers discarded due to low dosage ",
          "probability (< dosage_threshold = %.2f), resulting in BAF weight = 0 for those windows. ",
          "Only z-score will contribute to CN calling there.\n  Affected windows (Chr:WindowID): %s"
        ),
        verbose = verbose, level = 2, type = ">>",
        length(prob_filtered_idx),
        dosage_threshold,
        paste(prob_filtered_labels, collapse = ", ")
      )
    }

    # Count number of markers with BAF values in each window
    n_baf <- sapply(baf_list, function(x) if(any(is.na(x))) length(x[-which(is.na(x))]) else length(x))

    # BAF weights from heterozygote counts
    # Lower the number of heterozygotes in the window, lower the weight of the BAF emission
    # If there are no heterozygous, only z-score will be considered
    HET_N <- n_het_window
    if (any(HET_N > 0)) {
      ref_q <- suppressWarnings(quantile(HET_N[HET_N>0], het_quantile, na.rm=TRUE))
      if (!is.finite(ref_q) || ref_q <= 0) ref_q <- max(HET_N, na.rm=TRUE)
      if (!is.finite(ref_q) || ref_q <= 0) ref_q <- 1
      w_baf <- sqrt(pmin(1, HET_N / ref_q))
    } else {
      w_baf <- rep(0, W)
    }
    w_baf[!is.finite(w_baf)] <- 0
    if(round(exp_ploidy, 0) == 1) w_baf <- rep(1, W) # if expected ploidy is 1, we don't expect heterozygotes, so BAF likelihood should be taken fully into account regardless of heterozygote count
    w_baf <- w_baf * baf_weight
  } else {
    n_het_window <- NA
    ll_baf_matrix <- matrix(0, nrow = W, ncol = length(cn_grid))
    prob_baf_matrix <- matrix(0, nrow = W, ncol = length(cn_grid))
    w_baf <- rep(0, W)
    if(round(exp_ploidy, 0) == 1) {
      w_baf <- rep(1, W) # if expected ploidy is 1, we don't expect heterozygotes, so BAF likelihood should be taken fully into account regardless of heterozygote count
      w_baf <- w_baf * baf_weight
    }
  }
  vmsg("BAF weights defined", verbose = verbose, level = 1, type = ">>")

  # --- HMM for all chromosomes together ---
  # fit a single HMM across all windows (all chromosomes)
  vmsg("Setting initial and transitions HMM matrices", verbose = verbose, level = 0, type = ">>")

  # Setup for all windows
  K <- length(cn_grid)
  state_ids <- as.character(cn_grid)

  # Transition matrix and initial state distribution
  # A is K×K. Before normalizing, every off-diagonal entry is 0.001 (rare jumps), and the diagonal is 0.995 (strong tendency to stay in the same CN state).
  A <- matrix(1e-3, K, K); diag(A) <- transition_jump
  # A <- A / rowSums(A) makes each row sum to 1 (a valid Markov matrix). After this, values are almost unchanged (just scaled so each row sums exactly to 1).
  A <- A / rowSums(A)
  # To consider: If is possible to change this matrix to favor specific small jumps, like if changing from 2->4 is more likely than 2->6
  # pi0 is the initial state distribution at the first window; here it’s uniform (no prior preference for any CN at the start).
  # Set pi0 to favor exp_ploidy (if provided) or best CN from BAF model
  pi0 <- rep((1 - initial_prob) / (K - 1), K)
  if (!is.null(exp_ploidy) && !is.na(exp_ploidy)) {
    best_cn <- as.numeric(exp_ploidy)
  } else {
    best_cn <- as.numeric(selected_model$best$best_cn)
  }
  best_idx <- which(as.numeric(cn_grid) == best_cn)
  if (length(best_idx) == 1) {
    pi0[best_idx] <- 0.85
  } else {
    pi0 <- rep(1/K, K) # fallback to uniform if best CN not found
  }

  vmsg("Initial and transition HMM matrices ready", verbose = verbose, level = 1, type = ">>")


  vmsg("Defining z-score distribution templates", verbose = verbose, level = 0, type = ">>")
  mu <- define_z_limits(z = d[[z_col]], z_window = z, cn_grid = cn_grid,
                        exp_ploidy = exp_ploidy, z_range = z_range, verbose = verbose)

  # sig is the (shared) standard deviation of the z emission across states.
  # It starts at the sample SD of z, with a safety floor of 0.1 to avoid zero/near-zero variance that would blow up log-likelihoods.
  sig <- sd(z, na.rm = TRUE); if (!is.finite(sig) || sig <= 1e-6) sig <- 0.1
  W <- length(z)
  vmsg("Initial z-score mean by state: %s", verbose = verbose, level = 2, type = ">>", paste0(round(mu,3), collapse = ", "))
  vmsg("Initial z-score SD: %s", verbose = verbose, level = 2, type = ">>", sig)

  # --- EM loop ---
  vmsg("Starting EM loop", verbose = verbose, level = 0, type = ">>")

  rm_res <- em_hmm_cn(cn_grid, mu, K, state_ids, sig, z,
                      z_only, ll_baf_matrix, n_baf, w_baf,
                      correct_scale, A, pi0, W, max_iter, verbose)

  list2env(rm_res, envir = environment())
  vmsg("Updated z-score mean by state: %s", verbose = verbose, level = 2, type = ">>", paste(sprintf("CN%d: %.3f", cn_grid, mu), collapse = "; "))
  vmsg("Updated z-score SD: %s", verbose = verbose, level = 2, type = ">>", round(sig,3))

  # If z mean is not from the lowest to the highest follow lower ploidy to higher ploidy
  # It means that user tested unlikely ploidies, in this case, modify cn_grid and run again
  idx <- 0
  while(any(mu != sort(mu)) & idx < 10) {
    # Avoid infinite loop
    idx <- idx + 1
    # Identify valid ploidies - for lower ploidies than the expected should have lower mu, higher ploidies than expect should have higher mu
    # If not, the ploidy is not valid and should be removed from the grid, and the estimation should be rerun
    exp_idx <- which(as.numeric(names(mu)) == as.numeric(exp_ploidy))
    mu_exp <- mu[exp_idx]
    # For lower ploidies, keep only those with strictly decreasing mu
    lower_ploidies <- which(as.numeric(names(mu)) < as.numeric(exp_ploidy))
    keep_lower <- lower_ploidies[order(-as.numeric(names(mu)[lower_ploidies]))] # descending order
    if (length(keep_lower) > 0) {
      last_mu <- mu_exp
      valid_lower <- c()
      for (idx in keep_lower) {
        if (mu[idx] < last_mu) {
          valid_lower <- c(valid_lower, idx)
          last_mu <- mu[idx]
        }
      }
      lower_idx <- valid_lower
    } else {
      lower_idx <- integer(0)
    }
    # For higher ploidies, keep only those with strictly increasing mu
    higher_ploidies <- which(as.numeric(names(mu)) > as.numeric(exp_ploidy))
    keep_higher <- higher_ploidies[order(as.numeric(names(mu)[higher_ploidies]))]
    if (length(keep_higher) > 0) {
      last_mu <- mu_exp
      valid_higher <- c()
      for (idx in keep_higher) {
        if (mu[idx] > last_mu) {
          valid_higher <- c(valid_higher, idx)
          last_mu <- mu[idx]
        }
      }
      higher_idx <- valid_higher
    } else {
      higher_idx <- integer(0)
    }
    valid_idx <- c(lower_idx, exp_idx, higher_idx)
    valid_idx <- sort(valid_idx)
    cn_grid <- cn_grid[valid_idx]
    mu <- define_z_limits(z = d[[z_col]], z_window = z, cn_grid = cn_grid,
                          exp_ploidy = exp_ploidy, z_range_out = FALSE, verbose = verbose) # redefine mu with the new cn_grid, but without z_range to avoid changing the limits too much and keep the same order of the ploidies, which is already checked in the previous steps
    K <- length(cn_grid)
    ll_baf_matrix <- ll_baf_matrix[,valid_idx]
    pi0 <- pi0[valid_idx]
    state_ids <- state_ids[valid_idx]
    A <- A[valid_idx, valid_idx]
    vmsg("Some ploidies were removed due to non-monotonic z means", verbose = verbose, level = 2, type = ">>")
    vmsg("Rerunning EM with updated cn_grid", verbose = verbose, level = 1, type = ">>")
    rm_res <- em_hmm_cn(cn_grid, mu, K, state_ids, sig, z,
                        z_only, as.matrix(ll_baf_matrix), n_baf, w_baf,
                        correct_scale, as.matrix(A), pi0, W, max_iter, verbose)
    list2env(rm_res, envir = environment())
    vmsg("Updated z-score mean by state: %s", verbose = verbose, level = 2, type = ">>", paste(sprintf("CN%d: %.3f", cn_grid, mu), collapse = "; "))
    vmsg("Updated z-score SD: %s", verbose = verbose, level = 2, type = ">>", sprintf("%.3f", sig))
  }

  vmsg("EM loop complete", verbose = verbose, level = 1, type = ">>")

  vmsg("Decoding Viterbi path", verbose = verbose, level = 0, type = ">>")
  # Decode Viterbi path
  vit_path <- viterbi(ll_em, log(A), log(pi0))
  cn_call <- cn_grid[vit_path]

  # max posterior per window
  post_max <- apply(gamma, 1, max)

  vmsg("Viterbi path decoded", verbose = verbose, level = 1, type = ">>")

  # Prepare output
  vmsg("Preparing output", verbose = verbose, level = 0, type = ">>")
  post_df <- as.data.frame(gamma)
  names(post_df) <- paste0("post_CN", cn_grid)
  result <- cbind(
    data.frame(
      Sample    = sample_id,
      Chr       = win_df$Chr,
      WindowID  = win_df$WindowID,
      Start     = win_df$Start,
      End       = win_df$End,
      n_snps    = win_df$n_snps,
      n_het     = n_het_window,
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
    distribution = selected_model$best$dist,
    mu = mu,
    sigma = sig,
    A = A,
    pi0 = pi0,
    bins = M,
    bw = selected_model$best$bw,
    loglik = ll_hist[length(ll_hist)],
    z_range = z_range,
    het_quantile = het_quantile,
    baf_weight = baf_weight,
    transition_jump = transition_jump,
    z_only = z_only,
    exp_ploidy = exp_ploidy,
    rm_outliers = rm_outliers,
    outlier_alpha = outlier_alpha,
    segment_zscore = segment_zscore,
    snps_per_window = snps_per_window,
    min_snps_per_window = min_snps_per_window,
    add_uniform = add_uniform,
    uniform_weight = selected_model$best$uniform_weight
  )

  if (!is.null(result) && !is.null(d)) {
    map_df <- result[, c("Chr", "WindowID", "w_baf", "CN_call", "post_max")]
    names(map_df)[names(map_df) == "WindowID"] <- ".__w__"
    d <- merge(d, map_df, by = c("Chr", ".__w__"), all.x = TRUE, sort = FALSE)
    d <- d[order(match(seq_len(nrow(d)), as.integer(rownames(d)))), ]
    rownames(d) <- NULL
  }

  vmsg("Done!", verbose = verbose, level = 1, type = ">>")

  return(structure(list(by_window = result, by_marker = d, params = params), class = "hmm_CN"))
}

#' Run hmm_estimate_CN in parallel for multiple samples (using parLapply)
#'
#' Runs copy-number HMM for each sample in a user-defined vector (or all samples if "all" is specified).
#' Returns a combined data.frame with results for all samples.
#'
#' @param qploidy_standarize_result An object of class qploidy_standardization as returned by standardize().
#' @param sample_ids Character vector of sample IDs to analyze, or "all" for all samples in the data.
#' @param n_cores Number of cores to use (default: 2).
#' @param parallel_type Character. Parallel backend to use: \code{"FORK"}, \code{"PSOCK"}, or \code{"auto"} (default). \code{"auto"} selects \code{"FORK"} on Unix/macOS (faster; workers inherit the parent environment automatically) and \code{"PSOCK"} on Windows (the only option available there). Use \code{"PSOCK"} explicitly if you need cross-platform reproducibility or are debugging worker crashes.
#' @param ... Additional arguments passed to hmm_estimate_CN (e.g., chr, snps_per_window, etc).
#' @return A data.frame with results for all samples, as returned by hmm_estimate_CN$by_window, combined.
#'
#' @importFrom parallel makeCluster parLapply stopCluster clusterExport
#' @importFrom dplyr bind_rows
#'
#' @export
hmm_estimate_CN_multi <- function(qploidy_standarize_result = NULL,
                                  data = NULL,
                                  geno.pos = NULL,
                                  use_values = c("BAF", "zscore"),
                                  sample_ids = "all",
                                  n_cores = 2,
                                  parallel_type = "auto",
                                  ...) {
  # sanity check: require either a standardization object or raw data + geno.pos
  if (!is.null(qploidy_standarize_result)) {
    if (!inherits(qploidy_standarize_result, "qploidy_standardization")) {
      stop("Input must be a qploidy_standardization object as returned by standardize().")
    }
    all_samples <- unique(as.data.frame(qploidy_standarize_result$data)$SampleName)
  } else {
    if (is.null(data) || is.null(geno.pos)) {
      stop("When qploidy_standarize_result is NULL, both 'data' and 'geno.pos' must be provided.")
    }
    req_data <- c("MarkerName", "SampleName", "ratio", "R")
    req_gp   <- c("MarkerName", "Chromosome", "Position")
    miss_d  <- setdiff(req_data, names(data))
    miss_gp <- setdiff(req_gp,   names(geno.pos))
    if (length(miss_d))  stop(paste("'data' is missing columns:",   paste(miss_d,  collapse = ", ")))
    if (length(miss_gp)) stop(paste("'geno.pos' is missing columns:", paste(miss_gp, collapse = ", ")))
    if (!identical(use_values, c("ratio", "R"))) {
      stop("When qploidy_standarize_result is NULL, use_values must be c(\"ratio\", \"R\").")
    }
    all_samples <- unique(data$SampleName)
  }

  # resolve parallel backend
  allowed_types <- c("auto", "FORK", "PSOCK")
  if (!parallel_type %in% allowed_types)
    stop(sprintf("'parallel_type' must be one of: %s.", paste(allowed_types, collapse = ", ")))
  if (parallel_type == "auto")
    parallel_type <- if (.Platform$OS.type == "windows") "PSOCK" else "FORK"

  if (identical(sample_ids, "all")) {
    sample_ids <- all_samples
  } else {
    sample_ids <- intersect(sample_ids, all_samples)
    if (length(sample_ids) == 0) stop("No valid sample IDs found in input data.")
  }

  # capture dots ONCE
  dots <- list(...)

  cl <- makeCluster(n_cores, type = parallel_type)
  on.exit(stopCluster(cl), add = TRUE)

  if (parallel_type == "PSOCK") {
    # PSOCK workers are blank R processes: must reload packages and export symbols
    clusterEvalQ(cl, {
      options(warn = 1)
      library(methods)
      library(Qploidy)
      NULL
    })
    clusterExport(cl,
                  varlist = c("worker", "hmm_estimate_CN", "generate_baf_template",
                              "baf_log_likelihood", "logsumexp", "viterbi"),
                  envir = environment()
    )
  }
  # FORK workers inherit everything from the parent process — no export needed

  results_list <- parLapply(cl, sample_ids, worker,
                            qploidy_standarize_result, dots, data, geno.pos, use_values)

  # Re-issue warnings collected inside workers on the main session
  for (worker_res in results_list) {
    for (w in worker_res$warnings) warning(w, call. = FALSE)
  }
  results_list <- lapply(results_list, `[[`, "result")

  parameters <- lapply(results_list, function(x) x$params)
  names(parameters) <- sample_ids
  by_window <- lapply(results_list, function(x) x$by_window)
  by_marker <- lapply(results_list, function(x) x$by_marker)

  by_window <- Filter(Negate(is.null), by_window)
  if (length(by_window) == 0) stop("No results returned for any sample.")

  by_window <- bind_rows(by_window)
  by_marker <- bind_rows(by_marker)

  idx <- which(colnames(by_window) == "post_max")
  idx1 <- order(colnames(by_window)[(idx + 1):ncol(by_window)])
  by_window <- by_window[,c(1:idx, idx+idx1)]

  rownames(by_window) <- NULL
  return(structure(list(by_window =by_window, by_marker = by_marker, params_samples = parameters), class = "hmm_CN"))
}


#' Print method for hmm_CN objects
#'
#' Prints a summary of the HMM copy-number estimation results, including key parameters.
#' If the object contains a params_samples list (multi-sample), prints the number of samples.
#' Otherwise, prints details from the params list.
#'
#' @param x An object of class 'hmm_CN'.
#' @param ... Additional arguments (ignored).
#'
#' @method print hmm_CN
#'
#' @export
print.hmm_CN <- function(x, ...) {
  if (!inherits(x, "hmm_CN")) {
    stop("Object is not of class 'hmm_CN'.")
  }
  if (!is.null(x$params_samples)) {
    cat("hmm_CN multi-sample result\n")
    cat("  Number of samples:", length(x$params_samples), "\n")
    invisible(x)
    return()
  }
  params <- x$params
  cat("hmm_CN result\n")
  cat("  Copy-number grid:", paste(params$cn_grid, collapse=", "), "\n")
  cat("  Expected ploidy:", params$exp_ploidy, "\n")
  cat("  Minimum SNPs per window:", params$min_snps_per_window, "\n")
  cat("  Initial state probabilities (pi0):", paste(round(params$pi0, 3), collapse=", "), "\n")
  cat("  Estimated z means per CN:", paste(round(params$mu, 3), collapse=", "), "\n")
  cat("  Estimated z mean:", mean(x$by_window$z , na.rm=TRUE), "\n")
  cat("  Estimated z sigma:", round(params$sigma, 3), "\n")
  cat("  BAF Emission distribution:", if(!is.null(params$distribution)) params$distribution else "(not specified)", "\n")
  cat("  Final log-likelihood:", params$loglik, "\n")
  # Print range of CN_call values
  if (!is.null(x$by_window) && !is.null(x$by_window$CN_call)) {
    cn_vals <- x$by_window$CN_call
    cat("  Range of CN_call (window-level):", paste(range(as.numeric(cn_vals), na.rm=TRUE), collapse=" - "), "\n")
  }
  invisible(x)
}

#' Write hmm_CN object to three CSV files
#'
#' Writes the three main components of an hmm_CN object (by_window, by_marker, params) to CSV files.
#' The user must provide a file prefix; files will be named <prefix>_by_window.csv, <prefix>_by_marker.csv, <prefix>_params.csv.
#' @param hmm_CN An object of class 'hmm_CN'.
#' @param prefix File prefix for output files (character scalar).
#'
#' @importFrom data.table fwrite
#' @export
write_hmm_CN <- function(hmm_CN, prefix) {
  stopifnot(inherits(hmm_CN, "hmm_CN"))
  stopifnot(is.character(prefix) && length(prefix) == 1)
  fwrite(hmm_CN$by_window, paste0(prefix, "_by_window.csv.gz"), row.names = FALSE, compress = "gzip")
  fwrite(hmm_CN$by_marker, paste0(prefix, "_by_marker.csv.gz"), row.names = FALSE, compress = "gzip")
  # Save params or params_samples as RDS, depending on which is present
  if (!is.null(hmm_CN$params_samples)) {
    saveRDS(list(params_samples = hmm_CN$params_samples), file = paste0(prefix, "_params.rds"))
  } else if (!is.null(hmm_CN$params)) {
    saveRDS(hmm_CN$params, file = paste0(prefix, "_params.rds"))
  } else {
    stop("hmm_CN object must have either 'params' or 'params_samples'.")
  }
}

#' Read hmm_CN object from three files
#'
#' Reads the three main components of an hmm_CN object (by_window, by_marker, params or params_samples) from files.
#' The user must provide the three file paths explicitly.
#' @param by_window_file Path to the by_window CSV file.
#' @param by_marker_file Path to the by_marker CSV file.
#' @param params_file Path to the params RDS file.
#' @importFrom data.table fread
#'
#' @return An object of class 'hmm_CN'.
#' @export
read_hmm_CN <- function(by_window_file, by_marker_file, params_file) {
  by_window <- fread(by_window_file, data.table = FALSE)
  by_marker <- fread(by_marker_file, data.table = FALSE)

  # Input checks for by_window columns
  required_cols_window <- c(
    "Sample",      # Identifier for the sample
    "Chr",         # Chromosome name
    "WindowID",    # Unique identifier for the genomic window
    "Start",       # Start position of the window in base pairs
    "End",         # End position of the window in base pairs
    "n_snps",      # Number of SNPs within the window
    "n_het",       # Number of heterozygous SNPs within the window
    "z",           # Z-score for the window
    "w_baf",       # B-allele frequency HMM weight for the window
    "CN_call",     # Copy number call for the window
    "post_max"     # Maximum posterior probability for the copy number call
  )
  missing_cols_window <- setdiff(required_cols_window, colnames(by_window))
  if (length(missing_cols_window) > 0) {
    stop(sprintf(
      "by_window file is missing required column(s): %s",
      paste(missing_cols_window, collapse = ", ")
    ))
  }
  # Check for at least one post_CN column
  post_cn_cols <- grep("^post_CN[0-9]+$", colnames(by_window), value = TRUE)
  if (length(post_cn_cols) == 0) {
    stop("by_window file must contain at least one 'post_CN[number]' column (e.g., post_CN1, post_CN2, ...)")
  }

  # Input checks for by_marker columns
  required_cols_marker <- c(
    "Chr",        # Chromosome name
    ".__w__",          # Window
    "MarkerName", # Name of the marker
    "SampleName", # Identifier for the sample
    "X",          # X-coordinate for the marker
    "Y",          # Y-coordinate for the marker
    "R",          # Intensity value for the marker
    "ratio",      # Ratio value for the marker
    "geno",       # Genotype call for the marker
    "baf",        # B-allele frequency for the marker
    "z",          # Z-score for the marker
    "Position",   # Genomic position of the marker
    "outlier",    # Outlier status for the marker
    "w_baf",      # Weighted B-allele frequency for the marker
    "CN_call",    # Copy number call for the marker
    "post_max"    # Maximum posterior probability for the copy number call
  )
  missing_cols_marker <- setdiff(required_cols_marker, colnames(by_marker))
  if (length(missing_cols_marker) > 0) {
    stop(sprintf(
      "by_marker file is missing required column(s): %s",
      paste(missing_cols_marker, collapse = ", ")
    ))
  }

  # Or make it flexible to handle both local files and URLs
  if (grepl("^http", params_file)) {
    params_obj <- readRDS(gzcon(url(params_file)))
  } else {
    params_obj <- readRDS(params_file)
  }
  # Multi-sample: params_obj is a list with params_samples, or just the params_samples list itself
  required_params <- c(
    "cn_grid", "distribution", "mu", "sigma", "A", "pi0", "bins", "bw", "loglik", "z_range",
    "het_quantile", "baf_weight", "transition_jump", "z_only", "exp_ploidy", "rm_outliers",
    "outlier_alpha", "segment_zscore", "snps_per_window", "min_snps_per_window", "add_uniform", "uniform_weight"
  )

  check_params_list <- function(param_list) {
    missing <- setdiff(required_params, names(param_list))
    if (length(missing) > 0) {
      stop(sprintf(
        "params_file is missing required item(s) in a parameter list: %s",
        paste(missing, collapse = ", ")
      ))
    }
    TRUE
  }

  if (is.list(params_obj) && (!is.null(params_obj$params_samples) || (is.null(names(params_obj)) && all(sapply(params_obj, function(x) is.list(x) && all(required_params %in% names(x))))))) {
    params_samples <- if (!is.null(params_obj$params_samples)) params_obj$params_samples else params_obj
    # Check all internal lists
    lapply(params_samples, check_params_list)
    structure(list(by_window = by_window, by_marker = by_marker, params_samples = params_samples), class = "hmm_CN")
  } else if (is.list(params_obj) && !is.null(names(params_obj)) && all(required_params %in% names(params_obj))) {
    check_params_list(params_obj)
    structure(list(by_window = by_window, by_marker = by_marker, params = params_obj), class = "hmm_CN")
  } else {
    stop("params_file does not contain a recognizable params or params_samples object with all required items.")
  }
}
