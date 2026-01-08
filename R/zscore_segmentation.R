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
#' @param z_col Character string giving the column name of the numeric z-score
#'   vector to segment.
#'
#' @param chr_col Character string giving the column name identifying chromosomes.
#'
#' @param pos_col Character string giving the column name of genomic positions
#'   (used only for ordering markers).
#'
#' @param method Character string specifying the changepoint search method.
#'   One of `"PELT"` (default) or `"BinSeg"`. Passed to
#'   \code{\link[changepoint]{cpt.mean}}.
#'
#' @param penalty Character string specifying the penalty for adding changepoints.
#'   One of `"MBIC"` (default), `"BIC"`, `"AIC"`, or `"Manual"`.
#'
#' @param pen_value Numeric penalty value used only when `penalty = "Manual"`.
#'
#' @param minseglen Integer giving the minimum number of markers allowed per
#'   segment. Larger values produce fewer, broader segments.
#'
#' @param window_col Character string giving the name of the output column
#'   containing window (segment) indices.
#'
#' @param verbose Logical indicating whether to print progress messages.
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
#' A `data.frame` identical to `dat` but with an additional integer column
#' specified by `window_col` indicating changepoint window membership for
#' each marker.
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
#' @importFrom changepoint cpt.mean cpts
#'
#' @export
add_changepoint_windows <- function(dat,
                                    z_col = "z",
                                    chr_col = "Chr",
                                    pos_col = "Position",
                                    method = c("PELT","BinSeg"),
                                    penalty = c("MBIC","BIC","AIC","Manual"),
                                    pen_value = NULL,
                                    minseglen = 5,
                                    window_col = ".__w__",
                                    verbose = TRUE) {

  method  <- match.arg(method)
  penalty <- match.arg(penalty)

  stopifnot(z_col %in% names(dat))
  stopifnot(chr_col %in% names(dat))
  stopifnot(pos_col %in% names(dat))

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

    fit <- cpt.mean(
      z,
      method = method,
      penalty = penalty,
      pen.value = if (penalty == "Manual") pen_value else NULL,
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
  dat
}
