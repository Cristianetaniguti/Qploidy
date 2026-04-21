#' Heterozygosity Heatmap per Sample
#'
#' Computes the proportion of markers with BAF values within a heterozygosity
#' range (\code{het_range}) for each sample and displays the result as a
#' near-square heatmap. Hovering over a cell shows the sample name,
#' heterozygosity rate, and number of markers used (requires the \pkg{plotly}
#' package; a static \pkg{ggplot2} plot is returned if \pkg{plotly} is not
#' installed).
#'
#' @param input Accepts one of:
#'   1. A data.frame with at least columns \code{MarkerName}, \code{SampleName},
#'      and the column specified by \code{col2use}.
#'   2. An object of class \code{'qploidy_standardization'} (as returned by
#'      \code{standardize()}).
#'   3. A character string: path to a \code{.tsv[.gz]} file produced by
#'      \code{write_qploidy_standardization()}.
#' @param col2use Character. Column to use for computing heterozygosity.
#'   Either \code{"baf"} (default) or \code{"ratio"}.
#' @param het_range Numeric vector of length 2 giving the lower and upper
#'   bounds of \code{col2use} that define a heterozygous marker.
#'   Defaults to \code{c(0.2, 0.8)}.
#'
#' @return If \pkg{plotly} is available, a \code{plotly} htmlwidget with hover
#'   tooltips; otherwise a \code{ggplot2} object. Values are returned
#'   invisibly so the plot is still displayed when the function is called
#'   interactively.
#'
#' @examples
#' \dontrun{
#' plot_heterozygosity(my_standardization_object)
#' plot_heterozygosity("standardization.tsv.gz", het_range = c(0.3, 0.7))
#' plot_heterozygosity(my_standardization_object, col2use = "ratio")
#' }
#'
#' @importFrom dplyr group_by summarize
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn scale_x_continuous
#'   scale_y_continuous theme_minimal theme element_blank labs
#' @importFrom scales percent_format
#'
#' @export
plot_heterozygosity <- function(
    input,
    col2use = c("baf", "ratio"),
    het_range = c(0.2, 0.8)
) {
  ## --- validate arguments --------------------------------------------------
  col2use <- match.arg(col2use)

  if (!is.numeric(het_range) || length(het_range) != 2 ||
      het_range[1] >= het_range[2] ||
      het_range[1] < 0 || het_range[2] > 1) {
    stop("'het_range' must be a numeric vector of length 2 with 0 <= het_range[1] < het_range[2] <= 1.")
  }

  ## --- parse input (mirrors depth_pca_plot) --------------------------------
  if (inherits(input, "qploidy_standardization")) {
    tdata <- as.data.frame(input$data)
  } else if (is.character(input) && length(input) == 1 && file.exists(input)) {
    tdata <- as.data.frame(read_qploidy_standardization(input)$data)
  } else if (is.data.frame(input)) {
    tdata <- input
  } else {
    stop(paste(
      "'input' must be a data.frame with columns MarkerName, SampleName, and", col2use,
      "; a 'qploidy_standardization' object; or a valid file path."
    ))
  }

  required_cols <- c("MarkerName", "SampleName", col2use)
  missing_cols <- setdiff(required_cols, colnames(tdata))
  if (length(missing_cols) > 0)
    stop(sprintf(
      "Input data is missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))

  ## --- compute per-sample heterozygosity -----------------------------------
  het_df <- tdata %>%
    group_by(SampleName) %>%
    summarize(
      het_ratio  = mean(.data[[col2use]] >= het_range[1] & .data[[col2use]] <= het_range[2], na.rm = TRUE),
      n_markers  = sum(!is.na(.data[[col2use]])),
      .groups    = "drop"
    )

  het_df <- het_df[order(het_df$het_ratio, decreasing = TRUE), ]

  ## --- lay out samples in a near-square grid -------------------------------
  n     <- nrow(het_df)
  ncols <- ceiling(sqrt(n))
  nrows <- ceiling(n / ncols)

  het_df$col_idx <- ((seq_len(n) - 1L) %% ncols) + 1L
  het_df$row_idx <- nrows - ((seq_len(n) - 1L) %/% ncols)  # row 1 = bottom (incomplete), highest row = top

  ## --- build ggplot heatmap ------------------------------------------------
  p <- ggplot(
    het_df,
    aes(
      x    = col_idx,
      y    = row_idx,
      fill = het_ratio,
      text = paste0(
        "Sample: ",          SampleName,
        "\nHeterozygosity: ", round(het_ratio * 100, 2), "%",
        "\nMarkers used: ",   n_markers
      )
    )
  ) +
    geom_tile(color = "black", linewidth = 0.4) +
    scale_fill_gradientn(
      colours = c("#2166ac", "#fdae61", "#d7191c"),
      limits  = c(0, 1),
      name    = "Heterozygosity",
      labels  = percent_format(accuracy = 1)
    ) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.title   = element_blank(),
      axis.text    = element_blank(),
      panel.grid   = element_blank()
    ) +
    labs(title = sprintf(
      "Heterozygosity per Sample  [%s in (%.2f, %.2f)]  n = %d",
      col2use, het_range[1], het_range[2], n
    ))

  ## --- add hover if plotly is available ------------------------------------
  if (requireNamespace("plotly", quietly = TRUE)) {
    out <- plotly::ggplotly(p, tooltip = "text")
    # xgap/ygap are only valid on heatmap traces; apply selectively to avoid warnings
    heatmap_idx <- which(vapply(out$x$data,
                                function(tr) isTRUE(tr$type == "heatmap"),
                                logical(1)))
    if (length(heatmap_idx) > 0)
      out <- plotly::style(out, xgap = 1.5, ygap = 1.5, traces = heatmap_idx)
  } else {
    message("Install the 'plotly' package for interactive hover tooltips.")
    out <- p
  }

  invisible(out)
}
