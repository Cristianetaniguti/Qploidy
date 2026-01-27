# Suppress R CMD check notes for non-standard evaluation in dplyr/ggplot2
if (getRversion() >= "2.15.1") utils::globalVariables(
  c("starts_all", "starts", "outlier", "WindowID", "z_mean", "Chr", "region_id", "w_baf", "SampleName", "Position", "baf", "Start", "End", "prob_call", "CN_call")
)

#' Plot copy-number segments per window with posterior shading and BAF
#'
#' Plots a genome track for each window, showing the called copy number (CN) as a horizontal segment and coloring by the posterior probability of that call. Facets are split by chromosome. The top panel shows z-scores per window, colored by BAF weight; the bottom panel shows CN segments colored by posterior probability. The new top panel shows BAF values per SNP, colored by BAF weight, for the selected sample and chromosomes. If `qploidy_standarize_result` is NULL, marker-level data is taken from `hmm_CN$updated_data` (if available).
#'
#' @param hmm_CN An object of class \code{hmm_CN} (output from \code{hmm_estimate_CN}), containing a \code{result} data frame with one row per window and columns: \code{Sample}, \code{Chr}, \code{Start}, \code{End}, \code{CN_call}, \code{post_CN*}, etc. Optionally, may contain \code{updated_data} for marker-level fallback.
#' @param qploidy_standarize_result An object of class \code{qploidy_standardization} (output from \code{standardize}), used to extract BAF values for the BAF panel. If NULL, marker-level data is taken from \code{hmm_CN$updated_data}.
#' @param sample_id Character scalar. Which sample from \code{hmm_CN$Sample} to display. Defaults to the first unique value in \code{hmm_CN$Sample}.
#' @param cn_min,cn_max Numeric scalars. Y-axis limits for CN. Defaults span the min/max of \code{CN_call}.
#' @param show_window_lines Logical. If TRUE, show dashed vertical lines at window boundaries.
#' @param include_first_in_chr Logical. If TRUE, include the first window line in each chromosome.
#' @param line_color,line_alpha,line_width,line_linetype Appearance settings for window boundary lines.
#' @param heights Numeric vector of length 2 or 3. Relative heights of the BAF, z, and CN panels.
#' @param z_by_mean Logical. If TRUE, plot horizontal black lines for mean z per window using geom_segment; otherwise, use geom_smooth as before.
#' @param chr Character vector. If provided, only these chromosomes are displayed in the plots.
#'
#' @return A \code{ggplot} object (from ggpubr::ggarrange). Print to render, or add layers/scales as needed.
#'
#' @details
#' Posterior columns are detected by the prefix "post_CN" and matched to \code{CN_call} values, so the function is agnostic to the specific CN grid. The BAF panel shows per-SNP BAF values for the selected sample and chromosomes, colored by the BAF weight for the corresponding region. The z panel shows window z-scores, colored by BAF weight. The CN panel shows copy number segments colored by posterior probability. If marker-level data is unavailable in \code{qploidy_standarize_result}, the function will use \code{hmm_CN$updated_data} if present, otherwise an error is thrown.
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
#' # qploidy_standarize_result should be a
#' # qploidy_standardization object with $data containing BAF values
#' plot_cn_track(toy, qploidy_standarize_result, sample_id = "S1")
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @importFrom magrittr %>%
#' @importFrom ggpubr ggarrange
#' @importFrom scales squish
#' @importFrom purrr map
#'
#' @export
plot_cn_track <- function(hmm_CN,
                          qploidy_standarize_result = NULL,
                          sample_id = NULL,
                          cn_min = NULL,
                          cn_max = NULL,
                          show_window_lines = FALSE,
                          include_first_in_chr = FALSE,
                          line_color = "grey60",
                          line_alpha = 0.6,
                          line_width = 0.3,
                          line_linetype = "dashed",
                          heights = c(2, 2.5),
                          z_by_mean = TRUE,
                          chr = NULL) {

  stopifnot(inherits(hmm_CN, "hmm_CN"))
  stopifnot(inherits(qploidy_standarize_result, "qploidy_standardization") || is.null(qploidy_standarize_result))

  if (!is.null(hmm_CN$by_window)) {
    df <- hmm_CN$by_window
  } else {
    stop("hmm_CN object must have by_window or result data frame.")
  }

  if (is.null(sample_id)) sample_id <- unique(df$Sample)[1]
  x <- filter(df, Sample == sample_id)
  if (nrow(x) == 0) stop("No rows for sample_id = ", sample_id)

  # Filter chromosomes if chr argument is provided
  if (!is.null(chr)) {
    x <- x[x$Chr %in% chr, ]
    if (nrow(x) == 0) stop("No rows for selected chromosomes: ", paste(chr, collapse=", "))
  }

  # Remove NA chromosomes from all relevant data frames before plotting
  x <- x[!is.na(x$Chr), ]
  if (nrow(x) == 0) stop("No non-NA chromosomes left after filtering.")
  chr_order <- sort(unique(x$Chr), method = "radix")
  chr_order <- chr_order[order(as.numeric(gsub("[^0-9]+", "", chr_order)), chr_order)]
  x$Chr <- factor(x$Chr, levels = chr_order)

  # Marker-level data: prefer by_marker, then qploidy_standarize_result
  if (!is.null(hmm_CN$by_marker)) {
    data_sample2 <- hmm_CN$by_marker
    data_sample2 <- filter(data_sample2, SampleName == sample_id & Chr %in% unique(x$Chr))
    data_sample2 <- data_sample2[!is.na(data_sample2$Chr), ]
    data_sample2$Chr <- factor(data_sample2$Chr, levels = chr_order)
  } else if (!is.null(qploidy_standarize_result)) {
    data_sample2 <- qploidy_standarize_result$data
    data_sample2 <- filter(data_sample2, SampleName == sample_id & Chr %in% unique(x$Chr))
    data_sample2 <- data_sample2[!is.na(data_sample2$Chr), ]
    data_sample2$Chr <- factor(data_sample2$Chr, levels = chr_order)
  } else {
    stop("No marker-level data found in hmm_CN or qploidy_standarize_result.")
  }

  if (is.null(cn_min)) cn_min <- min(x$CN_call, na.rm = TRUE)
  if (is.null(cn_max)) cn_max <- max(x$CN_call, na.rm = TRUE)

  # compute P(CN_call) per window
  cn_states <- as.integer(gsub("^post_CN", "", grep("^post_CN", names(x), value = TRUE)))
  post_mat  <- as.matrix(x[, paste0("post_CN", cn_states), drop = FALSE])
  idx       <- match(x$CN_call, cn_states)
  x$prob_call <- post_mat[cbind(seq_len(nrow(x)), idx)]
  x$prob_call <- as.numeric(x$prob_call)
  if (!("prob_call" %in% names(x)) || all(is.na(x$prob_call)) || length(unique(x$prob_call)) == 1) {
    x$prob_call <- rep(1, nrow(x))
  }

  # coordinates
  x$Start <- as.numeric(x$Start)
  x$End   <- as.numeric(x$End)
  x$Mid   <- (x$Start + x$End)/2

  # Per-chromosome genomic range (based on BAF markers)
  chrom_ghost <- data_sample2 %>%
    group_by(Chr) %>%
    reframe(Position = range(Position, na.rm = TRUE), .groups = "drop")

  # One row per CN window per chromosome, ordered, with region_id
  win_tbl <- x %>%
    distinct(Chr, Start, End, w_baf) %>%
    arrange(Chr, Start) %>%
    group_by(Chr) %>%
    mutate(region_id = row_number()) %>%
    ungroup()

  # Breakpoints (window starts) per chromosome
  # Used to assign markers to region_id regardless of show_window_lines
  break_tbl <- win_tbl %>%
    group_by(Chr) %>%
    summarise(starts_all = list(sort(unique(Start))), .groups = "drop") %>%
    mutate(
      starts = map(starts_all, ~{
        s <- .x
        if (!include_first_in_chr) s <- s[s != min(s)]
        s
      })
    ) %>%
    select(Chr, starts)

  # vlines for plotting (only used when show_window_lines=TRUE)
  vlines <- win_tbl %>%
    distinct(Chr, Start) %>%
    arrange(Chr, Start)

  if (!include_first_in_chr) {
    vlines <- vlines %>%
      group_by(Chr) %>%
      filter(Start != min(Start)) %>%
      ungroup()
  }

  # Helper: tag markers with region_id based on breakpoints
  tag_markers_with_region <- function(marker_df, break_tbl) {
    marker_df %>%
      left_join(break_tbl, by = "Chr") %>%
      group_by(Chr) %>%
      mutate(
        region_id = {
          s <- starts[[1]]
          # if NA / missing starts, treat as one region
          if (length(s) == 0 || all(is.na(s))) {
            rep(1L, n())
          } else {
            brks <- c(-Inf, s, Inf)
            cut(Position, breaks = brks, labels = FALSE, right = FALSE)
          }
        }
      ) %>%
      ungroup() %>%
      select(-starts)
  }

  # Tag markers and join correct w_baf
  marker_df <- data_sample2 %>%
    tag_markers_with_region(break_tbl) %>%
    left_join(win_tbl %>% select(Chr, region_id),
                     by = c("Chr", "region_id"))

  # -------- top panel: z dots (color = w_baf, shape = outlier if present) --------
  shape_aes <- if ("outlier" %in% colnames(marker_df)) aes(shape = outlier) else NULL

  if (z_by_mean) {
    z_means <- hmm_CN$by_window %>%
      filter(Sample == sample_id & !is.na(Chr) & Chr %in% chr_order) %>%
      select(Chr, WindowID, Start, End, w_baf, z_mean = z) %>%
      arrange(factor(Chr, levels = chr_order), Start)
    z_means$Chr <- factor(z_means$Chr, levels = chr_order)
    p_z <- ggplot(marker_df, aes(x = Position, y = z, color = w_baf)) +
      geom_point(size = 1.2, alpha = 0.8, mapping = shape_aes) +
      geom_segment(
        data = z_means,
        aes(x = Start, xend = End, y = z_mean, yend = z_mean, group = WindowID),
        inherit.aes = FALSE,
        color = "black", linewidth = 0.7
      ) +
      geom_blank(
        data = chrom_ghost,
        inherit.aes = FALSE,
        aes(x = Position, y = min(marker_df$z, na.rm = TRUE))
      )
    if (show_window_lines) {
      p_z <- p_z + geom_vline(data = vlines, aes(xintercept = Start),
                   color = line_color, alpha = line_alpha,
                   linewidth = line_width, linetype = line_linetype)
    }
    if ("outlier" %in% colnames(marker_df)) {
      p_z <- p_z + scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 4), name = "Outlier")
    }
    p_z <- p_z +
      facet_wrap(~ Chr, scales = "free_x", nrow = 1) +
      scale_color_distiller(palette = "Spectral", direction = -1,
                            limits = c(0, 1),
                            oob = squish,
                            name = "BAF weight") +
      labs(x = NULL, y = "z") +
      theme_bw(base_size = 12) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        strip.background   = element_rect(fill = "grey95"),
        legend.position    = "none",
        axis.title.x       = element_blank(),
        axis.text.x        = element_text(angle = 30, vjust = 1, hjust = 1)
      )
  } else {
    p_z <- ggplot(marker_df, aes(x = Position, y = z, color = w_baf)) +
      geom_point(size = 1.2, alpha = 0.8, mapping = shape_aes) +
      geom_smooth(aes(group = Chr), method = "loess", se = FALSE,
                  color = "black", linewidth = 0.7, span = 0.2) +
      geom_blank(
        data = chrom_ghost,
        inherit.aes = FALSE,
        aes(x = Position, y = min(marker_df$z, na.rm = TRUE))
      )
    if (show_window_lines) {
      p_z <- p_z + geom_vline(data = vlines, aes(xintercept = Start),
                   color = line_color, alpha = line_alpha,
                   linewidth = line_width, linetype = line_linetype)
    }
    if ("outlier" %in% colnames(marker_df)) {
      p_z <- p_z + scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 4), name = "Outlier")
    }
    p_z <- p_z +
      facet_wrap(~ Chr, scales = "free_x", nrow = 1) +
      scale_color_distiller(palette = "Spectral", direction = -1,
                            limits = c(0, 1),
                            oob = squish,
                            name = "BAF weight") +
      labs(x = NULL, y = "z") +
      theme_bw(base_size = 12) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        strip.background   = element_rect(fill = "grey95"),
        legend.position    = "none",
        axis.title.x       = element_blank(),
        axis.text.x        = element_text(angle = 30, vjust = 1, hjust = 1)
      )
  }

  # -------- bottom panel: CN segments (color = P(CN call)) --------
  if (!("prob_call" %in% names(x))) stop("prob_call column missing in CN plot data.")
  x$prob_call <- as.numeric(x$prob_call)

  # Split x into single-marker and multi-marker windows
  x_seg <- x[x$Start != x$End, ]
  x_pt  <- x[x$Start == x$End, ]

  p_cn <- ggplot(x) +
    geom_blank(
      data = chrom_ghost,
      inherit.aes = FALSE,
      aes(x = Position, y = min(x$CN_call))
    ) +
    # Segments for multi-marker windows
    geom_segment(
      data = x_seg,
      aes(x = Start, xend = End, y = CN_call, yend = CN_call, color = prob_call),
      linewidth = 2, lineend = "butt"
    ) +
    # Points for single-marker windows
    geom_point(
      data = x_pt,
      aes(x = Start, y = CN_call, color = prob_call),
      size = 2
    )
  if (show_window_lines) {
    p_cn <- p_cn + geom_vline(data = vlines, aes(xintercept = Start),
                 color = line_color, alpha = line_alpha,
                 linewidth = line_width, linetype = line_linetype)
  }
  p_cn <- p_cn +
    facet_wrap(~ Chr, scales = "free_x", nrow = 1) +
    scale_y_continuous(breaks = seq(cn_min, cn_max, by = 1), minor_breaks = NULL) +
    scale_color_viridis_c(name = "P(CN call)", limits = c(0, 1), option = "D", direction = -1) +
    labs(x = "Genomic position (bp)", y = "Copy number") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      strip.background   = element_rect(fill = "grey95"),
      legend.position    = "bottom",
      axis.text.x        = element_text(angle = 30, vjust = 1, hjust = 1)
    )

  # -------- BAF panel (color = w_baf; uses SAME tagging/join as z panel) --------
  plot_df <- data_sample2 %>%
    tag_markers_with_region(break_tbl) %>%
    left_join(win_tbl %>% select(Chr, region_id),
                     by = c("Chr", "region_id"))

  p_baf <- ggplot(plot_df, aes(x = Position, y = baf, color = w_baf)) +
    geom_blank(
      data = chrom_ghost,
      inherit.aes = FALSE,
      aes(x = Position, y = 0)
    ) +
    geom_point(alpha = 0.7, size = 1)
  if (show_window_lines) {
    p_baf <- p_baf + geom_vline(data = vlines, aes(xintercept = Start),
               color = line_color, alpha = line_alpha,
               linewidth = line_width, linetype = line_linetype)
  }
  p_baf <- p_baf +
    facet_wrap(~Chr, scales = "free_x", nrow = 1) +
    theme_bw() +
    ylab("BAF") +
    scale_color_distiller(palette = "Spectral", direction = -1,
                          limits = c(0, 1),
                          oob = squish,
                          name = "BAF weight") +
    theme(
      axis.text.x        = element_text(angle = 30, vjust = 1, hjust = 1),
      legend.position    = "top",
      text               = element_text(size = 12),
      axis.title.x       = element_blank()
    )

  # -------- stack with ggpubr (no patchwork needed) --------
  return(list(
    baf = p_baf,
    z   = p_z,
    cn  = p_cn,
    arranged = ggarrange(
      p_baf, p_z, p_cn,
      ncol = 1, nrow = 3,
      heights = heights,
      align = "v"
    )
  ))
}


#' Compare CNV tracks across samples
#'
#' Plots CNV windows as horizontal segments for multiple samples, faceted by chromosome.
#' Segments are colored by copy number (CN_call) using a high-contrast palette:
#' baseline CN (most frequent, weighted by window length) is black; losses are blue;
#' gains are red. Transparency reflects posterior confidence (post_max).
#'
#' @param hmm_CN An object of class hmm_CN (output from hmm_estimate_CN), containing a result data frame with columns: Sample, Chr, Start, End, CN_call, post_max, etc.
#' @param samples_to_plot Character vector of sample IDs to include. If NULL, the first sample in the data is plotted.
#' @param chromosomes Optional character vector of chromosomes to include (matching result$Chr). If NULL, plots all chromosomes present after sample filtering.
#' @param facet_ncol Number of columns for facet_wrap. If NULL, will be determined automatically unless facet_nrow is set.
#' @param facet_nrow Number of rows for facet_wrap. If set, facet_ncol will be determined automatically unless both are set.
#'
#' @return A ggplot object.
#'
#' @importFrom dplyr filter group_by summarise arrange desc
#' @importFrom ggplot2 ggplot geom_segment facet_wrap
#' @importFrom ggplot2 scale_x_continuous scale_color_manual scale_alpha
#' @importFrom ggplot2 labs theme_bw theme
#' @importFrom ggplot2 element_blank element_rect element_text
#'
#' @export
compare_cn_track <- function(hmm_CN,
                             samples_to_plot = NULL,
                             chromosomes = NULL,
                             facet_ncol = NULL,
                             facet_nrow = NULL) {

  stopifnot(inherits(hmm_CN, "hmm_CN"))
  cnv_df <- hmm_CN$by_window

  req_cols <- c("Sample", "Chr", "Start", "End", "CN_call", "post_max")
  missing_cols <- setdiff(req_cols, colnames(cnv_df))
  if (length(missing_cols) > 0) {
    stop("hmm_CN$by_window is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # If samples_to_plot is NULL, use the first sample found
  if (is.null(samples_to_plot) || length(samples_to_plot) == 0) {
    samples_to_plot <- unique(cnv_df$Sample)[1]
  }

  # ---- local helper: dependency-free, dark, high-contrast palette with black baseline ----
  make_cn_palette_dark <- function(cn_levels, baseline) {

    cn_num <- sort(unique(as.integer(as.character(cn_levels))))
    cn_chr <- as.character(cn_num)

    below <- cn_num[cn_num < baseline]
    above <- cn_num[cn_num > baseline]

    vals <- setNames(rep(NA_character_, length(cn_num)), cn_chr)

    # baseline: gray
    if (as.character(baseline) %in% cn_chr) {
      vals[as.character(baseline)] <- "#888888"  # gray
    }

    # predefined dark ramps (long enough for typical CN ranges; saturates beyond)
    blues <- c(
      "#08306b", "#2171b5", "#6baed6"  # dark, medium, light blue
    )
    reds <- c(
      "#67000d", "#cb181d", "#fc9272"  # dark, medium, light red
    )
    # If more levels, interpolate with greater color jumps
    if (length(below) > length(blues)) {
      blues <- colorRampPalette(blues, space = "Lab")(length(below))
    } else {
      blues <- blues[seq_len(length(below))]
    }
    if (length(above) > length(reds)) {
      reds <- colorRampPalette(reds, space = "Lab")(length(above))
    } else {
      reds <- reds[seq_len(length(above))]
    }
    if (length(below) > 0) {
      vals[as.character(below)] <- blues
    }
    if (length(above) > 0) {
      vals[as.character(above)] <- reds
    }
    vals
  }

  # ---- filter data ----
  plot_df <- cnv_df |>
    filter(.data$Sample %in% samples_to_plot)

  if (!is.null(chromosomes)) {
    plot_df <- plot_df |>
      filter(.data$Chr %in% chromosomes)
  }

  if (nrow(plot_df) == 0) {
    stop("No rows left after filtering. Check samples_to_plot and chromosomes.")
  }

  # order samples top-to-bottom in the user-provided order (first listed = top)
  samples_present <- intersect(samples_to_plot, unique(plot_df$Sample))
  plot_df$Sample <- factor(plot_df$Sample, levels = rev(samples_present))

  # CN_call as integer -> ordered factor
  plot_df$CN_call <- as.integer(as.character(plot_df$CN_call))

  # baseline CN = most frequent CN, weighted by window length (End-Start)
  plot_df$w <- pmax(0, plot_df$End - plot_df$Start)
  cn_totals <- plot_df |>
    group_by(.data$CN_call) |>
    summarise(total_bp = sum(.data$w, na.rm = TRUE), .groups = "drop") |>
    arrange(desc(.data$total_bp))

  baseline_cn <- cn_totals$CN_call[1]

  cn_levels <- sort(unique(plot_df$CN_call))
  plot_df$CN_call <- factor(plot_df$CN_call, levels = cn_levels)

  # palette
  cn_pal <- make_cn_palette_dark(levels(plot_df$CN_call), baseline_cn)

  # Order chromosomes by extracting the numeric part and sorting accordingly
  chr_levels <- unique(plot_df$Chr)
  chr_nums <- suppressWarnings(as.numeric(gsub("[^0-9]+", "", chr_levels)))
  chr_order <- chr_levels[order(chr_nums, chr_levels)]
  plot_df$Chr <- factor(plot_df$Chr, levels = chr_order)

  # ---- build plot ----
  n_facets <- length(unique(plot_df$Chr))
  ncol_final <- facet_ncol
  nrow_final <- facet_nrow
  if (!is.null(facet_nrow) && is.null(facet_ncol)) {
    ncol_final <- ceiling(n_facets / facet_nrow)
  }
  if (!is.null(facet_ncol) && is.null(facet_nrow)) {
    nrow_final <- NULL
  }
  ggplot(plot_df) +
    geom_segment(
      aes(
        x = .data$Start, xend = .data$End,
        y = .data$Sample, yend = .data$Sample,
        color = .data$CN_call, alpha = .data$post_max
      ),
      linewidth = 4,
      lineend = "butt"
    ) +
    facet_wrap(~ Chr, scales = "free_x", ncol = ncol_final, nrow = nrow_final) +
    scale_x_continuous(
      labels = function(x) format(x / 1e6, trim = TRUE),
      name = "Position (Mb)"
    ) +
    scale_color_manual(
      values = cn_pal,
      drop = FALSE,
      name = "CN"
    ) +
    scale_alpha(range = c(0.35, 1), guide = "none") +
    labs(y = NULL) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = NA),
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 9)
    )
}
