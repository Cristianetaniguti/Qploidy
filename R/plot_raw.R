#' Plot Raw Allele Frequency and Depth Data
#'
#' Merges a raw allele-count data frame with a genomic-position data frame and
#' then produces one or more diagnostic plots for a single sample. Supported
#' plot types mirror the ratio/depth/heterozygosity plots available in
#' [plot_qploidy_standardization()], but operate entirely on raw data (no
#' standardization object required).
#'
#' @param data A `data.frame` with at minimum columns:
#'   \describe{
#'     \item{`MarkerName`}{Character. Marker identifier (used as the join key).}
#'     \item{`SampleName`}{Character. Sample identifier.}
#'     \item{`ratio`}{Numeric in \[0, 1\]. B-allele frequency / allele ratio.}
#'     \item{`R`}{Numeric. Total read depth.}
#'   }
#' @param geno.pos A `data.frame` with columns:
#'   \describe{
#'     \item{`MarkerName`}{Character. Marker identifier (join key).}
#'     \item{`Chromosome`}{Character. Chromosome / scaffold name (mapped to
#'       `Chr` after merging).}
#'     \item{`Position`}{Numeric. Base-pair position on the chromosome.}
#'   }
#' @param sample Character string. The `SampleName` value to plot.
#' @param chr Character or numeric vector of chromosomes to include. When
#'   `NULL` (default) all chromosomes present after merging are used.
#' @param type Character vector of plot types to produce. Aliases are accepted:
#'   \describe{
#'     \item{`"all"`}{All available plot types.}
#'     \item{`"ratio"`}{Scatter plot of allele ratio vs. genomic position.}
#'     \item{`"R"`}{Smoothed read-depth curve vs. genomic position.}
#'     \item{`"het"` / `"heterozygosity"` / `"heterozygous"`}{Proportion of
#'       heterozygous loci per genomic window.}
#'     \item{`"Ratio_hist"` / `"ratio_hist"` / `"ratio_histogram"` /
#'       `"Ratio_histogram"`}{Per-chromosome ratio histogram.}
#'     \item{`"Ratio_hist_overall"` / `"ratio_hist_overall"` /
#'       `"ratio_histogram_overall"`}{Genome-wide ratio histogram.}
#'   }
#'   Default is `"all"`.
#' @param dot.size Numeric. Point size for scatter plots. Default `1`.
#' @param font_size Numeric. Base font size for all plots. Default `12`.
#' @param window_size Numeric. Window width (bp) used for the heterozygosity
#'   plot. Default `2000000`.
#' @param het_interval Numeric. Half-width around 0 and 1 used to classify a
#'   locus as homozygous. Loci with `ratio <= het_interval` or
#'   `ratio >= 1 - het_interval` are treated as homozygous. Default `0.1`.
#' @param bins Integer. Number of histogram bins. Default `100`.
#' @param ploidy Integer. Expected ploidy; used to draw expected-peak lines on
#'   ratio histograms when `add_expected_peaks = TRUE`. Default `4`.
#' @param area_single Numeric. Half-width (in ploidy-scaled units) of the
#'   coloured area around each expected histogram peak. Default `0.75`.
#' @param add_expected_peaks Logical. Draw vertical lines at expected ratio
#'   peaks on histogram plots. Default `FALSE`.
#' @param bar_colour Colour for histogram bars. Default `NA` (ggplot2 default).
#'
#' @return A `ggarrange` object (from **ggpubr**) stacking all requested plots,
#'   annotated with the sample name at the top. Returns `NULL` invisibly if no
#'   plots could be produced.
#'
#' @details
#' The function first performs a left join of `data` onto `geno.pos` by
#' `MarkerName`, renaming `Chromosome` to `Chr`. Rows without a matching
#' position are dropped. Only the requested `sample` and `chr` subset is then
#' retained before plotting.
#'
#' Plot types and their behaviour:
#' \describe{
#'   \item{`ratio`}{`ratio` vs. `Position`, one facet per chromosome.}
#'   \item{`R`}{GAM-smoothed `R` vs. `Position` with a horizontal median line,
#'     one facet per chromosome.}
#'   \item{`het`}{`ratio` is binarised into heterozygous (inside
#'     `(het_interval, 1-het_interval)`) vs. homozygous, then the proportion of
#'     heterozygous loci is computed in non-overlapping windows of
#'     `window_size` bp.}
#'   \item{`Ratio_hist`}{Histogram of `ratio` per chromosome. Expected-peak
#'     lines are added when `add_expected_peaks = TRUE` and `ploidy` is set.}
#'   \item{`Ratio_hist_overall`}{Same as `Ratio_hist` but pooling all
#'     chromosomes into a single panel.}
#' }
#'
#' @importFrom ggpubr ggarrange annotate_figure text_grob
#' @importFrom dplyr left_join filter select group_by summarise rename mutate
#' @importFrom magrittr "%>%"
#' @import ggplot2
#'
#' @export
plot_raw <- function(data,
                     geno.pos,
                     sample,
                     chr = NULL,
                     type = "all",
                     dot.size = 1,
                     font_size = 12,
                     window_size = 2000000,
                     het_interval = 0.1,
                     bins = 50,
                     bar_colour = NA,
                     ploidy = 4,
                     area_single = 0.75,
                     add_expected_peaks = FALSE) {

  # ── Input validation ──────────────────────────────────────────────────────
  if (!is.data.frame(data))    stop("'data' must be a data.frame.")
  if (!is.data.frame(geno.pos)) stop("'geno.pos' must be a data.frame.")

  required_data    <- c("MarkerName", "SampleName", "ratio", "R")
  required_geopos  <- c("MarkerName", "Chromosome", "Position")
  missing_data     <- setdiff(required_data,   names(data))
  missing_geopos   <- setdiff(required_geopos, names(geno.pos))
  if (length(missing_data)   > 0) stop("'data' is missing columns: ",    paste(missing_data,   collapse = ", "))
  if (length(missing_geopos) > 0) stop("'geno.pos' is missing columns: ", paste(missing_geopos, collapse = ", "))
  if (missing(sample) || is.null(sample) || length(sample) != 1)
    stop("'sample' must be a single character string.")

  # ── Type alias normalisation ───────────────────────────────────────────────
  type_aliases <- c(
    "all"                    = "all",
    "ratio"                  = "ratio",
    "R"                      = "R",
    "het"                    = "het",
    "heterozygosity"         = "het",
    "heterozygous"           = "het",
    "Ratio_hist"             = "Ratio_hist",
    "ratio_hist"             = "Ratio_hist",
    "ratio_histogram"        = "Ratio_hist",
    "Ratio_histogram"        = "Ratio_hist",
    "Ratio_hist_overall"     = "Ratio_hist_overall",
    "ratio_hist_overall"     = "Ratio_hist_overall",
    "ratio_histogram_overall" = "Ratio_hist_overall"
  )
  type_resolved <- type_aliases[type]
  unknown <- is.na(type_resolved)
  if (any(unknown)) warning("Unknown plot type(s) ignored: ", paste(type[unknown], collapse = ", "))
  type <- unique(unname(type_resolved[!unknown]))

  # ── Merge data + geno.pos ─────────────────────────────────────────────────
  geno.pos_clean <- geno.pos %>%
    select(MarkerName, Chr = Chromosome, Position) %>%
    mutate(Position = as.numeric(Position))

  df <- data %>%
    filter(SampleName == sample) %>%
    left_join(geno.pos_clean, by = "MarkerName") %>%
    filter(!is.na(Chr), !is.na(Position))

  if (nrow(df) == 0) stop("No data found for sample '", sample, "' after merging with geno.pos.")

  # ── Chromosome filter ──────────────────────────────────────────────────────
  if (is.null(chr)) {
    chr <- sort(unique(df$Chr))
  } else if (is.numeric(chr)) {
    chr <- sort(unique(df$Chr))[chr]
  }
  if (!all(chr %in% unique(df$Chr))) stop("One or more requested chromosomes not found in data.")
  df <- df %>% filter(Chr %in% chr)

  # ── Colour palette (Okabe-Ito) ────────────────────────────────────────────
  col_bar <- "#0072B2"

  # ── Build baf_sample-style data.frame for plot_baf_hist ───────────────────
  # plot_baf_hist expects a column named 'sample' (used as x) and 'ratio'
  hist_df       <- df
  hist_df$sample <- hist_df$ratio   # alias required by plot_baf_hist internals

  # ── Initialise plot slots ─────────────────────────────────────────────────
  p_ratio <- p_R <- p_het <- p_ratio_hist <- p_ratio_hist_overall <- NULL

  # ── ratio scatter ─────────────────────────────────────────────────────────
  if (any(type == "all" | type == "ratio")) {
    p_ratio <- df %>%
      ggplot(aes(x = Position, y = ratio)) +
      geom_point(alpha = 0.7, size = dot.size, colour = col_bar) +
      facet_grid(. ~ Chr, scales = "free_x") +
      theme_bw(base_size = font_size) +
      ylab("Ratio") +
      theme(
        axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text       = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
  }

  # ── R (depth) smooth ──────────────────────────────────────────────────────
  if (any(type == "all" | type == "R")) {
    p_R <- df %>%
      ggplot(aes(x = Position, y = R)) +
      facet_grid(. ~ Chr, scales = "free_x") +
      geom_smooth(method = "gam", colour = col_bar) +
      geom_hline(yintercept = median(df$R, na.rm = TRUE), linetype = "dashed") +
      theme_bw(base_size = font_size) +
      ylab("Read depth (R)") +
      theme(
        axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text       = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
  }

  # ── Heterozygosity rate ───────────────────────────────────────────────────
  if (any(type == "all" | type == "het")) {
    df$het_rate <- NA
    df$het_rate[df$ratio <= het_interval | df$ratio >= 1 - het_interval] <- 0
    df$het_rate[df$ratio >  het_interval & df$ratio <  1 - het_interval] <- 1

    df_lst <- split.data.frame(df, df$Chr)
    for (j in seq_along(df_lst)) {
      df_lst[[j]]$window <- 1
      if (length(unique(df_lst[[j]]$Position)) == 1) next()
      starts <- seq(1, max(df_lst[[j]]$Position), window_size)
      ends   <- c(starts[-1] - 1, max(df_lst[[j]]$Position))
      for (i in seq_along(starts)) {
        df_lst[[j]]$window[df_lst[[j]]$Position >= starts[i] &
                             df_lst[[j]]$Position <= ends[i]] <- starts[i]
      }
    }
    df_het <- do.call(rbind, df_lst)

    p_het <- df_het %>%
      group_by(Chr, window) %>%
      summarise(
        prop_het = sum(het_rate, na.rm = TRUE) /
          length(which(het_rate == 1 | het_rate == 0)),
        .groups = "drop"
      ) %>%
      ggplot(aes(x = window, y = prop_het, colour = prop_het)) +
      geom_line() +
      scale_colour_gradient(low = "black", high = "red") +
      facet_grid(. ~ Chr, scales = "free_x") +
      theme_bw(base_size = font_size) +
      ylab("Proportion of heterozygous loci") +
      xlab("Position") +
      theme(
        axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position  = "none",
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text       = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
  }

  # ── Per-chromosome ratio histogram ────────────────────────────────────────
  if (any(type == "all" | type == "Ratio_hist")) {
    p_ratio_hist <- plot_baf_hist(
      data_sample       = hist_df,
      area_single       = area_single,
      ploidy            = ploidy,
      add_estimated_peaks = FALSE,
      add_expected_peaks  = add_expected_peaks,
      BAF_hist_overall  = FALSE,
      ratio             = TRUE,
      font_size         = font_size,
      bins              = bins,
      bar_colour        = bar_colour
    )
  }

  # ── Genome-wide ratio histogram ───────────────────────────────────────────
  if (any(type == "all" | type == "Ratio_hist_overall")) {
    p_ratio_hist_overall <- plot_baf_hist(
      data_sample       = hist_df,
      area_single       = area_single,
      ploidy            = ploidy,
      add_estimated_peaks = FALSE,
      add_expected_peaks  = add_expected_peaks,
      BAF_hist_overall  = TRUE,
      ratio             = TRUE,
      font_size         = font_size,
      bins              = bins,
      bar_colour        = bar_colour
    )
  }

  # ── Assemble & return ─────────────────────────────────────────────────────
  p_all <- list(p_het, p_ratio, p_ratio_hist, p_ratio_hist_overall, p_R)
  p_all <- p_all[!sapply(p_all, is.null)]

  if (length(p_all) == 0) {
    warning("No plots produced for sample '", sample, "'.")
    return(invisible(NULL))
  }

  p_result <- ggarrange(plotlist = p_all, ncol = 1)
  p_result <- annotate_figure(
    p_result,
    top = text_grob(sample, face = "bold", size = 14)
  )
  return(p_result)
}

#' Select the Best BAF Model from Raw Allele-Count Data
#'
#' Convenience wrapper around [select_best_baf_model()] that accepts the raw
#' `data` data.frame used by [plot_raw()], filters it to a single sample, and
#' passes the `ratio` column as the `baf_vec` argument. All other arguments are
#' forwarded verbatim to [select_best_baf_model()].
#'
#' @param data A `data.frame` containing at minimum columns `SampleName` and
#'   `ratio`. Additional columns (e.g. `R`, `MarkerName`) are ignored.
#' @param sample Character string. The value of `SampleName` to select.
#' @param cn_grid Integer vector of copy-number states to test (e.g. `2:6`).
#' @param dists Character vector of distribution families to test. Passed to
#'   [select_best_baf_model()]. Default:
#'   `c("gaussian", "beta", "beta_binomial", "negative_binomial")`.
#' @param reflect Logical. Apply boundary reflection for continuous kernels.
#'   Default `TRUE`.
#' @param bw_grid Numeric vector of bandwidth values to search over. Default
#'   `c(0.02, 0.03, 0.04)`.
#' @param add_uniform_grid Logical vector. Whether to include a uniform noise
#'   component. Default `FALSE`.
#' @param uniform_weight_grid Numeric vector of uniform mixture weights to try.
#'   Default `c(0.01, 0.03, 0.05, 0.10, 0.15)`.
#' @param M Integer. Number of BAF histogram bins. Default `100`.
#' @param plot Logical. If `TRUE`, a ggplot for the best model is included in
#'   the returned object. Default `FALSE`.
#' @param param_count Optional named integer vector of free-parameter counts per
#'   distribution, used for BIC penalisation. Default `NULL` (zero for all).
#' @param count_grid_as_params Logical. Add +1 BIC penalty per tuned
#'   hyperparameter. Default `TRUE`.
#' @param min_het_frac Numeric. Minimum heterozygous fraction required to
#'   exclude CN=1 from `cn_grid`. Default `0.05`.
#' @param het_range Numeric vector of length 2. BAF interval defining
#'   heterozygous loci. Default `c(0.2, 0.8)`.
#'
#' @return An object of class `selected_BAF_model` as returned by
#'   [select_best_baf_model()].
#'
#' @seealso [select_best_baf_model()], [plot_raw()]
#'
#' @export
select_best_raw_model <- function(
    data,
    sample,
    cn_grid,
    dists               = c("gaussian", "beta", "beta_binomial", "negative_binomial"),
    reflect             = TRUE,
    bw_grid             = c(0.02, 0.03, 0.04),
    add_uniform_grid    = FALSE,
    uniform_weight_grid = c(0.01, 0.03, 0.05, 0.10, 0.15),
    M                   = 100,
    plot                = FALSE,
    param_count         = NULL,
    count_grid_as_params = TRUE,
    min_het_frac        = 0.05,
    het_range           = c(0.2, 0.8)) {

  if (!is.data.frame(data))
    stop("'data' must be a data.frame.")
  if (!all(c("SampleName", "ratio") %in% names(data)))
    stop("'data' must contain columns 'SampleName' and 'ratio'.")
  if (!is.character(sample) || length(sample) != 1)
    stop("'sample' must be a single character string.")

  df_s <- data[data$SampleName == sample, , drop = FALSE]
  if (nrow(df_s) == 0)
    stop("Sample '", sample, "' not found in 'data$SampleName'.")

  baf_vec <- df_s$ratio

  select_best_baf_model(
    baf_vec             = baf_vec,
    cn_grid             = cn_grid,
    dists               = dists,
    reflect             = reflect,
    bw_grid             = bw_grid,
    add_uniform_grid    = add_uniform_grid,
    uniform_weight_grid = uniform_weight_grid,
    M                   = M,
    plot                = plot,
    param_count         = param_count,
    count_grid_as_params = count_grid_as_params,
    min_het_frac        = min_het_frac,
    het_range           = het_range
  )
}
