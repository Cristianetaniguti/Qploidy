#' Segment ordered genomic z-scores into changepoint-based windows
#'
#' This function applies 1D changepoint detection on an ordered vector of
#' z-scores (typically one sample) to identify genomic segments with distinct
#' mean values. Segmentation is performed independently per chromosome and
#' markers are assigned to integer-valued windows representing changepoint
#' segments.
#'
#' The input data are first ordered by chromosome and genomic position for
#' segmentation, but the original row order is preserved in the output.
#'
#' If no changepoints are detected for a chromosome, all markers on that
#' chromosome are assigned to a single window.
#'
#' @param dat A `data.frame` containing marker-level data ordered or unordered.
#'   Must include z-scores, chromosome labels, and genomic positions.
#'
#' @param minseglen Integer giving the minimum number of markers allowed per
#'   segment. Larger values produce fewer, broader segments.
#'
#' @param verbose Logical indicating whether to print progress messages.
#'
#' @param plot Logical. If TRUE, returns a list with the segmented data and a ggplot object visualizing z-scores, window means, and changepoint boundaries per chromosome. If FALSE (default), returns only the segmented data.
#'
#' @details
#' Segmentation is performed independently for each chromosome using
#' mean-based changepoint detection via the \pkg{changepoint} package.
#'
#' Window indices restart at 1 within each chromosome. Markers assigned the
#' same window index belong to the same changepoint-defined genomic segment.
#'
#' This function is commonly used for detecting regional shifts in signal
#' intensity such as copy number variation proxies, allele balance deviations,
#' or other genome-wide z-score profiles.
#'
#' @return
#' If `plot = FALSE`, a `data.frame` identical to `dat` but with an additional
#' integer column specified by `window_col` indicating changepoint window
#' membership for each marker. If `plot = TRUE`, a list with two elements:
#'   \describe{
#'     \item{windows}{The segmented data.frame as above.}
#'     \item{plot}{A ggplot object visualizing z-scores, window means, and changepoint boundaries per chromosome.}
#'   }
#'
#' @details
#' If `plot = TRUE`, the plot shows:
#'   - z-score points for each marker (per chromosome)
#'   - horizontal lines for the mean z-score in each window
#'   - dashed vertical lines at changepoint boundaries
#'   - one facet per chromosome
#' The plot is suitable for visual inspection of segmentation quality and changepoint locations.
#'
#' @seealso
#' \code{\link[changepoint]{cpt.mean}},
#' \code{\link[DNAcopy]{segment}}
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   Chr = rep("chr1", 10),
#'   Position = seq(1, 1e6, length.out = 10),
#'   z = c(rep(0, 4), rep(1.2, 3), rep(0, 3))
#' )
#'
#' out <- add_changepoint_windows(df, minseglen = 2)
#' table(out$.__w__)
#' }
#'
#' @importFrom changepoint cpt.meanvar cpts
#' @importFrom dplyr group_by summarise arrange slice ungroup
#' @importFrom ggplot2 ggplot aes geom_point geom_segment geom_vline facet_wrap labs theme_bw theme element_blank element_rect element_text
#'
#' @export
add_changepoint_windows <- function(dat,
                                    minseglen = NULL,
                                    verbose = TRUE,
                                    plot = FALSE) {
  # Required columns
  required_cols <- c("z", "Chr", "Position")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns in input data: %s", paste(missing_cols, collapse=", ")))
  }
  z_col <- "z"
  chr_col <- "Chr"
  pos_col <- "Position"
  window_col <- ".__w__"

  # keep original row order, but segment on sorted genomic order within chr
  dat$.orig_row__ <- seq_len(nrow(dat))
  dat <- dat[order(dat[[chr_col]], dat[[pos_col]]), ]

  dat[[window_col]] <- NA_integer_

  # segment per chromosome
  for (chr in unique(dat[[chr_col]])) {
    idx <- which(dat[[chr_col]] == chr)
    z <- dat[[z_col]][idx]

    # handle tiny chromosomes / all-NA robustly
    if (length(idx) < 2 || all(is.na(z))) {
      dat[[window_col]][idx] <- 1L
      if (isTRUE(verbose)) {
        message(sprintf("[changepoint] %s: 1 window (insufficient data)", chr))
      }
      next
    }

    if(is.null(minseglen)){
      floor <- 5 # hard default
      frac <- 0.1 # hard default
      m <- max(floor, floor(length(z) * frac))
      minseglen <- min(m, floor(length(z) / 2))
    }

    fit <- cpt.meanvar(
      z,
      method = "PELT",
      penalty = "Manual",
      pen.value = log(length(z)),
      minseglen = minseglen,
      class = TRUE
    )

    ends <- cpts(fit)
    if (length(ends) == 0) {
      ends <- length(z)
    }
    starts <- c(1, head(ends, -1) + 1)

    w <- integer(length(z))
    for (k in seq_along(starts)) {
      w[starts[k]:ends[k]] <- k
    }
    # Ensure window indices are consecutive and start at 1 (no 0s)
    if (min(w) == 0) {
      w <- w + 1L
    }
    # Remap window indices so that the first marker is always window 1, and indices increase in genomic order
    unique_windows <- unique(w)
    new_indices <- match(w, unique_windows)
    w <- new_indices
    dat[[window_col]][idx] <- w

    if (isTRUE(verbose)) {
      message(sprintf(
        "[changepoint] %s: %d window%s detected",
        chr, length(starts), ifelse(length(starts) == 1, "", "s")
      ))
    }
  }

  # restore original row order, drop helper column
  dat <- dat[order(dat$.orig_row__), ]
  dat$.orig_row__ <- NULL

  plot_obj <- NULL
  if (plot) {
    df <- dat
    win_summary <- group_by(df, .data[[chr_col]], .data[[window_col]]) %>%
      summarise(
        start = min(.data[[pos_col]]),
        end   = max(.data[[pos_col]]),
        mean_z = mean(.data[[z_col]], na.rm = TRUE),
        .groups = "drop"
      )
    win_boundaries <- win_summary %>%
      arrange(.data[[chr_col]], start) %>%
      group_by(.data[[chr_col]]) %>%
      slice(-1) %>%        # drop first window per chromosome
      ungroup()
    plot_obj <- ggplot(df, aes(x = .data[[pos_col]], y = .data[[z_col]])) +
      geom_point(size = 0.8, alpha = 0.6) +
      geom_segment(
        data = win_summary,
        aes(x = start, xend = end, y = mean_z, yend = mean_z, group = .data[[chr_col]]),
        inherit.aes = FALSE,
        linewidth = 0.8,
        color = "blue"
      ) +
      geom_vline(
        data = win_boundaries,
        aes(xintercept = start, group = .data[[chr_col]]),
        linetype = "dashed",
        linewidth = 0.4,
        alpha = 0.6
      ) +
      facet_wrap(~.data[[chr_col]], scales = "free_x", ncol = 2) +
      labs(x = "Genomic position", y = "z-score") +
      theme_bw() +
      theme(
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold")
      )
  }
  if (plot) return(list(windows = dat, plot = plot_obj))
  dat
}
