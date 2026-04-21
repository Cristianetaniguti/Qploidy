#' Export hmm_dosage_calls to VCF
#'
#' Converts a hmm_dosage_calls object (data.frame) to a VCF file with GT, CN, AD, BAF, Z, post_max_CN, and post_max_dosage in FORMAT.
#'
#' @param hmm_dosage_calls A data.frame with columns: MarkerName, SampleName, Chr, Position, X, Y, CN_call, post_max_CN, dosage, post_max_dosage, BAF, z
#' @param file Path to output VCF file
#'
#' @importFrom data.table data.table setDT dcast fwrite setcolorder
#'
#' @details The VCF INFO field will contain the mode of CN_call values for each marker. The FORMAT field will include GT (genotype), CN (copy number call), AD (allelic depths), BAF, Z, PMC (posterior max CN), and PMD (posterior max dosage). Genotypes are assigned based on the dosage and CN_call, with 1 representing the alternate allele and 0 representing the reference allele.
#'
#' @return The path to the written VCF file (invisible)
#'
#' @export
export_VCF <- function(hmm_dosage_calls, file) {
  stopifnot(is.data.frame(hmm_dosage_calls))
  required_cols <- c("MarkerName", "SampleName", "Chr", "Position", "X", "Y", "CN_call", "post_max_CN", "dosage", "post_max_dosage", "baf", "z")
  missing_cols <- setdiff(required_cols, names(hmm_dosage_calls))
  if (length(missing_cols) > 0) stop(sprintf("Missing columns: %s", paste(missing_cols, collapse=", ")))

  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install the 'data.table' package for fast VCF export.")
  setDT(hmm_dosage_calls)
  samples <- unique(hmm_dosage_calls$SampleName)

  # Precompute FORMAT fields for all rows (vectorized, no vapply)
  hmm_dosage_calls[, GT := ifelse(!is.na(CN_call) & !is.na(dosage),
    sapply(seq_len(.N), function(i) {
      cn <- as.integer(CN_call[i]); dosage <- as.integer(dosage[i])
      if (is.na(cn) || is.na(dosage)) return("./.")
      gt_vec <- sort(c(rep(1, dosage), rep(0, cn - dosage)))
      paste(gt_vec, collapse="/")
    }), "./.")]
  hmm_dosage_calls[, AD := ifelse(is.na(X) | is.na(Y), ".,.", paste0(round(as.numeric(X),2), ",", round(as.numeric(Y),2)))]
  hmm_dosage_calls[, BAF := ifelse(is.na(baf), ".", format(round(as.numeric(baf),2), nsmall=2))]
  hmm_dosage_calls[, Z := ifelse(is.na(z), ".", format(round(as.numeric(z),2), nsmall=2))]
  hmm_dosage_calls[, PMC := ifelse(is.na(post_max_CN), ".", format(round(as.numeric(post_max_CN),3), nsmall=3))]
  hmm_dosage_calls[, PMD := ifelse(is.na(post_max_dosage), ".", format(round(as.numeric(post_max_dosage),3), nsmall=3))]
  hmm_dosage_calls[, CN_call := ifelse(is.na(CN_call), ".", as.character(CN_call))]
  hmm_dosage_calls[, FORMAT := paste(GT, CN_call, AD, BAF, Z, PMC, PMD, sep=":")]

  # Use data.table for wide format, but process in chunks if needed
  vcf_dt <- dcast(
    hmm_dosage_calls,
    Chr + Position + MarkerName ~ SampleName,
    value.var = "FORMAT",
    fill = "./:.:.,.:.:.:.:."
  )

  # INFO: CN is the mode of all CN_call values for this marker
  vcf_dt[, INFO := {
    cn_vals <- hmm_dosage_calls[.BY, on=.(Chr, Position, MarkerName)]$CN_call
    mode_cn <- as.integer(names(sort(table(cn_vals), decreasing=TRUE))[1])
    paste0("CN=", mode_cn)
  }, by=.(Chr, Position, MarkerName)]
  vcf_dt[, `:=`(REF = ".", ALT = ".", QUAL = ".", FILTER = ".", FORMAT = "GT:CN:AD:BAF:Z:PMC:PMD")]

  # Set column order for VCF
  setcolorder(vcf_dt, c("Chr", "Position", "MarkerName", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples))

  # Prepare VCF header
  vcf_header <- c(
    "##fileformat=VCFv4.2",
    "##INFO=<ID=CN,Number=1,Type=Integer,Description=Copy number call>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype based on dosage and CN_call>",
    "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=Copy number call>",
    "##FORMAT=<ID=AD,Number=2,Type=Integer,Description=Allelic depths: X,Y>",
    "##FORMAT=<ID=BAF,Number=1,Type=Float,Description=BAF value>",
    "##FORMAT=<ID=Z,Number=1,Type=Float,Description=Z value>",
    "##FORMAT=<ID=PMC,Number=1,Type=Float,Description=Posterior probability of CN_call>",
    "##FORMAT=<ID=PMD,Number=1,Type=Float,Description=Posterior probability of dosage call>"
  )
  vcf_col_header <- paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", paste(samples, collapse="\t"))
  vcf_header <- c(vcf_header, vcf_col_header)

  # Write header and body efficiently
  writeLines(vcf_header, file)
  fwrite(vcf_dt, file = file, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE, showProgress = FALSE)
  invisible(file)
}
#' Assign BAF Dosages Using Selected Model
#'
#' Given a vector of BAF values and a selected_BAF_model object (from select_best_baf_model),
#' assigns the most likely dosage (0:CN) to each BAF value, and returns the dosage probabilities.
#'
#' @param baf_vec Numeric vector of BAF values (in [0,1]).
#' @param selected_model Object of class 'selected_BAF_model' (output of select_best_baf_model). If NULL, user must provide bw, dist, add_uniform, uniform_weight, and cn.
#' @param bw Bandwidth (SD or concentration) for kernel smoothing (required if selected_model is NULL).
#' @param dist Distribution type for template peaks: "gaussian", "beta", "beta_binomial", or "negative_binomial" (required if selected_model is NULL).
#' @param add_uniform Logical. Whether to add a uniform noise component (required if selected_model is NULL).
#' @param uniform_weight Numeric in [0,1]. Mixture weight for uniform component (required if selected_model is NULL and add_uniform is TRUE).
#' @param cn Integer. Copy number (required if selected_model is NULL).
#' @param plot Logical. If TRUE, generates a plot of BAF values colored by maximum probability.
#'
#' @importFrom ggplot2 ggplot aes geom_density geom_point geom_vline scale_color_gradient labs theme_bw
#'
#'
#' @return A data.frame with columns: BAF, dosage (integer), max_prob (probability), and a matrix of probabilities for each possible dosage. If plot=TRUE, returns a list with data and plot.
#' @export

call_BAF_dosages <- function(baf_vec, selected_model = NULL, bw = NULL, dist = NULL,
                             add_uniform = NULL, uniform_weight = NULL, cn = NULL, plot = FALSE) {
  stopifnot(is.numeric(baf_vec))
  if (!is.null(selected_model)) {
    stopifnot(inherits(selected_model, "selected_BAF_model"))
    best <- selected_model$best
    if (is.null(best)) stop("selected_model$best is NULL; no model selected.")
    if(is.null(cn)) cn <- as.integer(best$best_cn)
    if(length(cn) > 1) stop("Provide a single CN value")
    bw <- best$bw
    dist <- best$dist
    add_uniform <- best$add_uniform
    uniform_weight <- best$uniform_weight
  } else {
    # Require all parameters if selected_model is NULL
    if (is.null(bw) || is.null(dist) || is.null(add_uniform) || is.null(cn)) {
      stop("If selected_model is NULL, you must provide bw, dist, add_uniform, and cn.")
    }
    if (isTRUE(add_uniform) && is.null(uniform_weight)) {
      stop("If add_uniform is TRUE, you must provide uniform_weight.")
    }
    if (!isTRUE(add_uniform)) uniform_weight <- 0
  }
  M <- 100
  dens <- generate_baf_template(
    cn = cn,
    M = M,
    bw = bw,
    dist = dist,
    reflect = TRUE,
    add_uniform = add_uniform,
    uniform_weight = uniform_weight
  )
  peak_pos <- seq(0, 1, length.out = cn + 1)
  prob_mat <- sapply(baf_vec, function(baf) {
    probs <- numeric(length = cn + 1)
    for (d in 0:cn) {
      if (dist == "gaussian") {
        probs[d + 1] <- dnorm(baf, mean = d / cn, sd = bw)
      } else if (dist == "beta") {
        kappa <- max(2, round(1 / max(bw, 1e-6)))
        mu <- d / cn
        if (mu <= 0) {
          a <- 1; b <- 1 + kappa
        } else if (mu >= 1) {
          a <- 1 + kappa; b <- 1
        } else {
          a <- 1 + mu * kappa
          b <- 1 + (1 - mu) * kappa
        }
        probs[d + 1] <- dbeta(baf, shape1 = a, shape2 = b)
      } else {
        probs[d + 1] <- exp(-((baf - d / cn)^2) / (2 * bw^2))
      }
    }
    if (isTRUE(add_uniform) && uniform_weight > 0) {
      probs <- (1 - uniform_weight) * probs + uniform_weight * (1 / (cn + 1))
    }
    probs <- probs / sum(probs)
    probs
  })
  prob_mat <- t(prob_mat)
  colnames(prob_mat) <- paste0("dosage_", 0:cn)
  max_idx <- max.col(prob_mat, ties.method = "first") - 1
  max_prob <- apply(prob_mat, 1, max)
  out <- data.frame(
    BAF = baf_vec,
    dosage = max_idx,
    max_prob = max_prob
  )
  out <- cbind(out, prob_mat)
  if (plot) {
    df_plot <- out
    # Create rainbow palette for dosages
    n_dosages <- length(unique(df_plot$dosage))
    rainbow_cols <- rainbow(n_dosages)
    names(rainbow_cols) <- sort(unique(df_plot$dosage))
    # Function to blend color with gray by max_prob using base R only
    blend_color <- function(dosage, max_prob) {
      base_col <- rainbow_cols[as.character(dosage)]
      gray_col <- "#000000" # darker gray for higher contrast
      # Convert hex to RGB
      hex_to_rgb <- function(hex) {
        hex <- gsub("#", "", hex)
        r <- strtoi(substr(hex, 1, 2), 16L)
        g <- strtoi(substr(hex, 3, 4), 16L)
        b <- strtoi(substr(hex, 5, 6), 16L)
        c(r, g, b)
      }
      rgb_base <- hex_to_rgb(base_col)
      rgb_gray <- hex_to_rgb(gray_col)
      rgb_blend <- round(max_prob * rgb_base + (1 - max_prob) * rgb_gray)
      sprintf("#%02X%02X%02X", rgb_blend[1], rgb_blend[2], rgb_blend[3])
    }
    df_plot$color_blend <- mapply(blend_color, df_plot$dosage, df_plot$max_prob)
    plot_obj <- ggplot(df_plot, aes(x = BAF, y = after_stat(density))) +
      geom_density(data = df_plot, aes(x = BAF), color = NA, fill = "#848d98", alpha = 0.3) +
      geom_point(aes(y = 0, color = color_blend), size = 6) +
      geom_vline(xintercept = peak_pos, color = "black", linetype = "dashed", size = 0.7) +
      scale_color_identity() +
      labs(title = paste0("BAF Dosage Assignment (CN=", cn, ", dist=", dist, ")"),
           x = "BAF",
           y = "Density",
           color = "Dosage (rainbow, gray=low prob)") +
      theme_bw()
    out <- list(data = out, plot = plot_obj)
  } else {
    out <- list(data = out)
  }
  class(out) <- c("BAF_dosage_calls", class(out))
  return(out)
}


#' Assign BAF dosages to markers in hmm_CN object
#'
#' For each marker in hmm_CN$by_marker, assigns dosage and dosage probability using call_BAF_dosages
#' with the CN_call for that marker and either the selected_model or custom parameters.
#'
#' @param hmm_CN An object of class 'hmm_CN' (output of hmm_estimate_CN)
#' @param selected_model Optional. An object of class 'selected_BAF_model'. If not provided, custom arguments must be supplied.
#' @param bw Bandwidth for kernel smoothing (required if selected_model is not provided)
#' @param dist Distribution family (required if selected_model is not provided)
#' @param add_uniform Logical, add uniform noise (required if selected_model is not provided)
#' @param uniform_weight Numeric, weight for uniform component (required if selected_model is not provided)
#' @param ... Additional arguments passed to call_BAF_dosages
#' @return A data.frame with columns: MarkerName, SampleName, Chr, Position, X, Y, CN_call, post_max_CN, dosage, post_max_dosage
#' @export
call_hmm_dosages <- function(hmm_CN, selected_model = NULL, bw = NULL, dist = NULL, add_uniform = NULL, uniform_weight = NULL, ...) {
  if (!inherits(hmm_CN, "hmm_CN")) stop("Input must be an object of class 'hmm_CN'.")
  d <- hmm_CN$by_marker
  if (is.null(d)) stop("hmm_CN object must have a by_marker data.frame.")
  required_cols <- c("MarkerName", "SampleName", "Chr", "Position", "X", "Y", "baf", "CN_call", "post_max")
  missing_cols <- setdiff(required_cols, names(d))
  if (length(missing_cols) > 0) stop(sprintf("Missing columns in by_marker: %s", paste(missing_cols, collapse=", ")))

  # Prepare output columns
  out <- d[, c("MarkerName", "SampleName", "Chr", "Position", "X", "Y", "CN_call", "post_max")]
  names(out)[names(out) == "post_max"] <- "post_max_CN"

  # For each marker, assign dosage and dosage probability
  baf_vec <- d$baf
  cn_vec <- as.numeric(as.character(d$CN_call))
  # Use selected_model if provided, else require all custom args
  if (!is.null(selected_model)) {
    if (!inherits(selected_model, "selected_BAF_model")) stop("selected_model must be of class 'selected_BAF_model'.")
    bw_val <- selected_model$best$bw
    dist_val <- selected_model$best$dist
    add_uniform_val <- selected_model$best$add_uniform
    uniform_weight_val <- selected_model$best$uniform_weight
  } else {
    if (is.null(bw) || is.null(dist) || is.null(add_uniform) || is.null(uniform_weight)) {
      stop("If selected_model is not provided, must supply bw, dist, add_uniform, and uniform_weight.")
    }
    bw_val <- bw
    dist_val <- dist
    add_uniform_val <- add_uniform
    uniform_weight_val <- uniform_weight
  }

  d_call <- split(d, d$CN_call)

  res <- mapply(function(baf, cn) {
      r <- call_BAF_dosages(baf$baf, cn=cn, bw=bw_val, plot=FALSE, dist=dist_val,
                       add_uniform=add_uniform_val, uniform_weight=uniform_weight_val)

    cbind(baf, dosage = r$data$dosage, post_max_dosage = r$data$max_prob)
  }, d_call, as.numeric(names(d_call)), SIMPLIFY = FALSE)

  d_dosages <- do.call(rbind, res)

  out <- d_dosages[,c("MarkerName", "SampleName", "Chr", "Position","X", "Y", "baf","z", "CN_call", "post_max", "dosage", "post_max_dosage")]

  colnames(out)[colnames(out) == "post_max"] <- "post_max_CN"
  out <- out[order(out$SampleName, out$Chr, out$Position),]
  class(out) <- c("hmm_dosage_calls", class(out))

  out
}
