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
#' @param z_by_mean Logical. If TRUE, plot horizontal black lines for mean z per window using geom_segment; otherwise, use geom_smooth as before.
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
                          heights = c(2, 2.5),
                          z_by_mean = TRUE) {

  stopifnot(inherits(hmm_CN, "hmm_CN"))
  stopifnot(inherits(qploidy_standarize_result, "qploidy_standardization"))

  # packages used
  # (keep as imports in your pkg / DESCRIPTION; here assumed available)
  # dplyr, tidyr, purrr, ggplot2, ggpubr, scales

  df <- hmm_CN$result

  # defaults BEFORE filtering
  if (is.null(sample_id)) sample_id <- unique(df$Sample)[1]
  x <- filter(df, Sample == sample_id)
  if (nrow(x) == 0) stop("No rows for sample_id = ", sample_id)

  # Sort chromosomes in natural order (e.g., chr1e, chr2e, chr10e) without extra dependencies
  chr_order <- sort(unique(x$Chr), method = "radix")
  chr_order <- chr_order[order(as.numeric(gsub("[^0-9]+", "", chr_order)), chr_order)]
  x$Chr <- factor(x$Chr, levels = chr_order)

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

  # Pull marker-level data for the sample
  data_sample2 <- qploidy_standarize_result$data %>%
    filter(SampleName == sample_id & Chr %in% unique(x$Chr))
  data_sample2$Chr <- factor(data_sample2$Chr, levels = chr_order)

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
    left_join(win_tbl %>% select(Chr, region_id, w_baf),
                     by = c("Chr", "region_id"))

  # OPTIONAL: if any w_baf are slightly out of [0,1], squish for color scaling
  # marker_df$w_baf <- squish(marker_df$w_baf, range = c(0, 1))

  # -------- top panel: z dots (color = w_baf) --------
  if (z_by_mean) {
    # Calculate mean z by WindowID using hmm_CN$result
    z_means <- hmm_CN$result %>%
      filter(Sample == sample_id) %>%
      select(Chr, WindowID, Start, End, w_baf, z_mean = z) %>%
      arrange(factor(Chr, levels = chr_order), Start)
    z_means$Chr <- factor(z_means$Chr, levels = chr_order)
    p_z <- ggplot(marker_df, aes(x = Position, y = z, color = w_baf)) +
      geom_point(size = 1.2, alpha = 0.8) +
      geom_segment(
        data = z_means,
        aes(x = Start, xend = End, y = z_mean, yend = z_mean, group = WindowID, color = NULL),
        inherit.aes = FALSE,
        color = "black", linewidth = 0.7
      ) +
      geom_blank(
        data = chrom_ghost,
        inherit.aes = FALSE,
        aes(x = Position, y = min(marker_df$z, na.rm = TRUE))
      ) +
      { if (show_window_lines)
        geom_vline(data = vlines, aes(xintercept = Start),
                   color = line_color, alpha = line_alpha,
                   linewidth = line_width, linetype = line_linetype)
        else NULL } +
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
      geom_point(size = 1.2, alpha = 0.8) +
      geom_smooth(aes(group = Chr), method = "loess", se = FALSE,
                  color = "black", linewidth = 0.7, span = 0.2) +
      geom_blank(
        data = chrom_ghost,
        inherit.aes = FALSE,
        aes(x = Position, y = min(marker_df$z, na.rm = TRUE))
      ) +
      { if (show_window_lines)
        geom_vline(data = vlines, aes(xintercept = Start),
                   color = line_color, alpha = line_alpha,
                   linewidth = line_width, linetype = line_linetype)
        else NULL } +
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
    ) +
    { if (show_window_lines)
      geom_vline(data = vlines, aes(xintercept = Start),
                          color = line_color, alpha = line_alpha,
                          linewidth = line_width, linetype = line_linetype)
      else NULL } +
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
    left_join(win_tbl %>% select(Chr, region_id, w_baf),
                     by = c("Chr", "region_id"))

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
                          linewidth = line_width, linetype = line_linetype)
      else NULL } +
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
  ggarrange(
    p_baf, p_z, p_cn,
    ncol = 1, nrow = 3,
    heights = heights,
    align = "v"
  )
}
