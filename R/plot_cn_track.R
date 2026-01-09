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
                          heights = c(2, 2.5)) {

  stopifnot(inherits(hmm_CN, "hmm_CN"))
  stopifnot(inherits(qploidy_standarize_result, "qploidy_standardization"))

  df <- hmm_CN$result

  # defaults BEFORE filtering
  if (is.null(sample_id)) sample_id <- unique(df$Sample)[1]
  x <- filter(df, Sample == sample_id)
  if (nrow(x) == 0) stop("No rows for sample_id = ", sample_id)

  if (is.null(cn_min)) cn_min <- min(x$CN_call, na.rm = TRUE)
  if (is.null(cn_max)) cn_max <- max(x$CN_call, na.rm = TRUE)

  # compute P(CN_call) per window
  cn_states <- as.integer(gsub("^post_CN", "", grep("^post_CN", names(x), value = TRUE)))
  post_mat  <- as.matrix(x[, paste0("post_CN", cn_states), drop = FALSE])
  idx       <- match(x$CN_call, cn_states)
  x$prob_call <- post_mat[cbind(seq_len(nrow(x)), idx)]
  # Ensure prob_call is numeric and not all NA/constant
  x$prob_call <- as.numeric(x$prob_call)
  if (!("prob_call" %in% names(x)) || all(is.na(x$prob_call)) || length(unique(x$prob_call)) == 1) {
    x$prob_call <- rep(1, nrow(x))
  }

  # coordinates
  x$Start <- as.numeric(x$Start)
  x$End   <- as.numeric(x$End)
  x$Mid   <- (x$Start + x$End)/2

  # dashed verticals (one per unique Start, per Chr)
  vlines <- NULL
  if (show_window_lines) {
    vlines <- x |>
      distinct(Chr, Start) |>
      arrange(Chr, Start)
    if (!include_first_in_chr) {
      vlines <- vlines |>
        group_by(Chr) |>
        filter(Start != min(Start)) |>
        ungroup()
    }
  }

  data_sample2 <- qploidy_standarize_result$data %>% filter(SampleName == sample_id & Chr %in% unique(x$Chr))

  # Per-chromosome genomic range (based on BAF; you can also combine with Start/End if you prefer)
  chrom_ghost <- data_sample2 %>%
    group_by(Chr) %>%
    reframe(Position = range(Position, na.rm = TRUE), .groups = "drop")
  # This gives 2 rows per Chr: min and max Position

  # -------- top panel: z dots (fill = w_baf) --------
  # Join window w_baf to each marker by region
  marker_df <- data_sample2
  # Assign region_id to each marker as in BAF panel
  if(nrow(x) == 1 || all(table(x$Chr) == 1)) {
    marker_df <- marker_df %>% group_by(Chr) %>% mutate(region_id = 1) %>% ungroup()
    w_baf_tbl <- x %>% select(Chr, w_baf) %>% mutate(region_id = 1)
  } else {
    marker_df <- marker_df %>%
      group_by(Chr) %>%
      group_modify(function(df, key) {
        chr_starts <- vlines %>%
          filter(Chr == key$Chr[1]) %>%
          arrange(Start) %>%
          pull(Start)
        brks <- c(-Inf, chr_starts, Inf)
        df$region_id <- cut(df$Position, breaks = brks, labels = FALSE, right = FALSE)
        df
      }) %>%
      ungroup()
    # Define regions_tbl and w_baf_tbl for multi-window chromosomes
    regions_tbl <- vlines %>%
      arrange(Chr, Start) %>%
      count(Chr, name = "n_vlines") %>%
      mutate(n_regions = n_vlines + 1L) %>%
      select(Chr, n_regions) %>%
      rowwise() %>%
      mutate(region_id = list(seq_len(n_regions))) %>%
      unnest(region_id) %>%
      ungroup()
    w_baf_tbl <- regions_tbl %>%
      group_by(Chr) %>%
      mutate(w_baf = x$w_baf[region_id]) %>%
      ungroup()
  }
  marker_df <- marker_df %>% left_join(w_baf_tbl, by = c("Chr", "region_id"))

  p_z <- ggplot(marker_df, aes(x = Position, y = z, color = w_baf)) +
    geom_point(size = 1.2, alpha = 0.8) +
    geom_smooth(aes(group = Chr), method = "loess", se = FALSE, color = "black", linewidth = 0.7, span = 0.2) +
    { if (show_window_lines)
      geom_vline(data = vlines, aes(xintercept = Start),
                 color = line_color, alpha = line_alpha,
                 linewidth = line_width, linetype = line_linetype)
      else NULL } +
    facet_wrap(~ Chr, scales = "free_x", nrow = 1) +
    scale_color_distiller(palette = "RdBu", direction = -1, limits = c(0, 1), name = "BAF weight") +
    labs(x = NULL, y = "z", title = sample_id) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      strip.background   = element_rect(fill = "grey95"),
      legend.position    = "none",
      axis.title.x       = element_blank(),
      axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)
    )

  # -------- bottom panel: CN segments (color = P(CN call)) --------
  # Ensure prob_call is numeric and present
  if (!("prob_call" %in% names(x))) stop("prob_call column missing in CN plot data.")
  x$prob_call <- as.numeric(x$prob_call)

  p_cn <- ggplot(x) +
    geom_blank(
      data = chrom_ghost,
      inherit.aes = FALSE,
      aes(x = Position, y = min(x$CN_call))
    ) +
    geom_segment(aes(x = Start, xend = End, y = CN_call, yend = CN_call, color = prob_call),
                 linewidth = 2, lineend = "butt") +
    { if (show_window_lines)
      geom_vline(data = vlines, aes(xintercept = Start),
                 color = line_color, alpha = line_alpha,
                 linewidth = line_width, linetype = line_linetype) else NULL } +
    facet_wrap(~ Chr, scales = "free_x", nrow = 1) +
    scale_y_continuous(breaks = seq(cn_min, cn_max, by = 1), minor_breaks = NULL) +
    scale_color_viridis_c(name = "P(CN call)", limits = c(0, 1), option = "D") +
    labs(x = "Genomic position (bp)", y = "Copy number") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      strip.background   = element_rect(fill = "grey95"),
      legend.position    = "bottom",
      axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)
    )

  # -------- BAF -----

  # Defensive: if only one window per chromosome, region_id should be 1 for all per chromosome
  if(nrow(x) == 1 || all(table(x$Chr) == 1)) {
    # Assign region_id = 1 for all points per chromosome
    data_tagged <- data_sample2 %>% group_by(Chr) %>% mutate(region_id = 1) %>% ungroup()
    # Build w_baf_tbl for each chromosome
    w_baf_tbl <- x %>% select(Chr, w_baf) %>% mutate(region_id = 1)
  } else {
    data_tagged <- data_sample2 %>%
      group_by(Chr) %>%
      group_modify(function(df, key) {
        chr_starts <- vlines %>%
          filter(Chr == key$Chr[1]) %>%
          arrange(Start) %>%
          pull(Start)
        brks <- c(-Inf, chr_starts, Inf)
        df$region_id <- cut(df$Position, breaks = brks, labels = FALSE, right = FALSE)
        df
      }) %>%
      ungroup()

    regions_tbl <- vlines %>%
      arrange(Chr, Start) %>%
      count(Chr, name = "n_vlines") %>%
      mutate(n_regions = n_vlines + 1L) %>%
      select(Chr, n_regions) %>%
      rowwise() %>%
      mutate(region_id = list(seq_len(n_regions))) %>%
      unnest(region_id) %>%
      ungroup()

    w_baf_tbl <- regions_tbl %>%
      group_by(Chr) %>%
      mutate(w_baf = x$w_baf[region_id]) %>%
      ungroup()
  }

  plot_df <- data_tagged %>%
    left_join(w_baf_tbl, by = c("Chr", "region_id"))

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
                 linewidth = line_width, linetype = line_linetype) else NULL } +
    facet_wrap(~Chr, scales = "free_x", nrow=1) +
    theme_bw() +
    ylab("BAF") +
    scale_color_distiller(palette = "RdBu", direction = -1, limits = c(0, 1), name = "BAF weight")+
    theme(
      axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
      legend.position = "top",
      text = element_text(size = 12),
      axis.title.x       = element_blank()
    )

  # -------- stack with ggpubr (no patchwork needed) --------
  # align = "v" keeps x-axes aligned; widths are matched automatically
  ggarrange(
    p_baf, p_z, p_cn,
    ncol = 1, nrow = 3,
    heights = heights,
    align = "v"
  )
}