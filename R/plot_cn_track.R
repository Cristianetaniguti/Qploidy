# Suppress R CMD check notes for non-standard evaluation in dplyr/ggplot2
if (getRversion() >= "2.15.1") utils::globalVariables(
  c("starts_all", "starts", "outlier", "WindowID", "z_mean", "Chr", "region_id", "w_baf", "SampleName", "Position", "baf", "Start", "End", "prob_call", "CN_call")
)

##' Plot copy-number segments per window with posterior shading and BAF
##'
##' Produces a multi-panel genome track visualization for a selected sample, including:
##' - A BAF distribution density plot (all markers, all chromosomes/windows)
##' - Stacked legends for BAF weight and P(CN call)
##' - A summary panel with CN grid, distribution, mu, sigma, sample name, and log-likelihood
##' - BAF panel (per-SNP BAF values, colored by BAF weight, with a gray density background)
##' - Z-score panel (window z-scores, colored by BAF weight)
##' - CN panel (copy number segments, colored by posterior probability)
##'
##' The panels are arranged as: first row (BAF distribution, legends, summary), second row (BAF), third row (z), fourth row (CN).
##' If `qploidy_standarize_result` is NULL, marker-level data is taken from `hmm_CN$updated_data` (if available).
##'
##' @param hmm_CN An object of class \code{hmm_CN} (output from \code{hmm_estimate_CN}), containing a \code{result} data frame with one row per window and columns: \code{Sample}, \code{Chr}, \code{Start}, \code{End}, \code{CN_call}, \code{post_CN*}, etc. Optionally, may contain \code{updated_data} for marker-level fallback.
##' @param qploidy_standarize_result An object of class \code{qploidy_standardization} (output from \code{standardize}), used to extract BAF values for the BAF panel. If NULL, marker-level data is taken from \code{hmm_CN$updated_data}.
##' @param data Optional. A \code{data.frame} with columns \code{MarkerName}, \code{SampleName}, \code{X}, \code{Y}, \code{R}, and \code{ratio}. Required when \code{qploidy_standarize_result} is \code{NULL}.
##' @param geno.pos Optional. A \code{data.frame} with columns \code{MarkerName}, \code{Chromosome}, and \code{Position}. Required when \code{qploidy_standarize_result} is \code{NULL}.
##' @param sample_id Character scalar. Which sample from \code{hmm_CN$Sample} to display. Defaults to the first unique value in \code{hmm_CN$Sample}.
##' @param cn_min,cn_max Numeric scalars. Y-axis limits for CN. Defaults span the min/max of \code{CN_call}.
##' @param show_window_lines Logical. If TRUE, show dashed vertical lines at window boundaries.
##' @param include_first_in_chr Logical. If TRUE, include the first window line in each chromosome.
##' @param line_color,line_alpha,line_width,line_linetype Appearance settings for window boundary lines.
##' @param heights Numeric vector of length 2 or 3. Relative heights of the BAF, z, and CN panels.
##' @param z_by_mean Logical. If TRUE, plot horizontal black lines for mean z per window using geom_segment; otherwise, use geom_smooth as before.
##' @param chr Character vector. If provided, only these chromosomes are displayed in the plots.
##' @param depth_zero_as_x Logical. If TRUE, markers with R == 0 are shown as 'x' in the z-score panel; otherwise, they are shown as normal dots. Default is FALSE.
##' @param summarized Logical. If TRUE, only summary lines/polygons are plotted (no individual points). Default is FALSE.
##'
##' @return A \code{ggplot} object (from ggpubr::ggarrange). Print to render, or add layers/scales as needed.
##'
##' @details
##' Posterior columns are detected by the prefix "post_CN" and matched to \code{CN_call} values, so the function is agnostic to the specific CN grid.
##' The BAF panel shows per-SNP BAF values for the selected sample and chromosomes, colored by the BAF weight for the corresponding region, and includes a gray density background for BAF distribution.
##' The z panel shows window z-scores, colored by BAF weight. Outliers (if present) are marked as 'x'. If `depth_zero_as_x = TRUE`, markers with R == 0 are also shown as 'x' and distinguished in the legend.
##' The CN panel shows copy number segments colored by posterior probability. The top summary panel displays the CN grid, emission distribution, estimated mu and sigma, sample name, and final log-likelihood.
##' If marker-level data is unavailable in \code{qploidy_standarize_result}, the function will use \code{hmm_CN$updated_data} if present, otherwise an error is thrown.
##'
##' @section Expected columns:
##' The function assumes posterior columns named exactly \code{post_CN<k>} for each copy-number state \code{k}. If your column naming differs, rename them before calling this function.
##'
##' @examples
##' \dontrun{
##' library(dplyr)
##' library(ggplot2)
##' toy <- data.frame(
##'   Sample   = "S1",
##'   Chr      = c("chr1","chr1","chr1"),
##'   Start    = c(1, 1e6, 2e6),
##'   End      = c(1e6-1, 2e6-1, 3e6-1),
##'   CN_call  = c(2,3,2),
##'   post_CN2 = c(0.95, 0.05, 0.9),
##'   post_CN3 = c(0.05, 0.94, 0.1)
##' )
##' # qploidy_standarize_result should be a
##' # qploidy_standardization object with $data containing BAF values
##' plot_cn_track(toy, qploidy_standarize_result, sample_id = "S1", depth_zero_as_x = TRUE)
##' }
##'
##' @import ggplot2
##' @import dplyr
##' @import tidyr
##' @importFrom magrittr %>%
##' @importFrom ggpubr ggarrange get_legend
##' @importFrom scales squish
##' @importFrom stats density
##' @importFrom purrr map
##'
##' @export
plot_cn_track <- function(hmm_CN,
                          qploidy_standarize_result = NULL,
                          data = NULL,
                          geno.pos = NULL,
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
                          chr = NULL,
                          summarized = FALSE,
                          depth_zero_as_x = FALSE) {

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

  # Marker-level data: prefer by_marker, then qploidy_standarize_result, then data+geno.pos
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
  } else if (!is.null(data) && !is.null(geno.pos)) {
    req_data <- c("MarkerName", "SampleName")
    req_gp   <- c("MarkerName", "Chromosome", "Position")
    miss_d  <- setdiff(req_data, names(data))
    miss_gp <- setdiff(req_gp,   names(geno.pos))
    if (length(miss_d))  stop(paste("'data' is missing columns:",   paste(miss_d,  collapse = ", ")))
    if (length(miss_gp)) stop(paste("'geno.pos' is missing columns:", paste(miss_gp, collapse = ", ")))
    gp <- geno.pos[, req_gp]
    colnames(gp)[colnames(gp) == "Chromosome"] <- "Chr"
    data_sample2 <- merge(
      data[data$SampleName == sample_id, , drop = FALSE],
      gp, by = "MarkerName", all.x = FALSE
    )
    data_sample2 <- data_sample2[data_sample2$Chr %in% unique(x$Chr) & !is.na(data_sample2$Chr), , drop = FALSE]
    data_sample2$Chr      <- factor(data_sample2$Chr, levels = chr_order)
    data_sample2$Position <- as.numeric(data_sample2$Position)
    # Alias ratio -> baf and R -> z so downstream plotting code works
    if (!"baf" %in% names(data_sample2) && "ratio" %in% names(data_sample2))
      data_sample2$baf <- data_sample2$ratio
    if (!"z"   %in% names(data_sample2) && "R"     %in% names(data_sample2))
      data_sample2$z   <- data_sample2$R
  } else {
    stop("No marker-level data found. Provide qploidy_standarize_result, hmm_CN$by_marker, or both 'data' and 'geno.pos'.")
  }

  data_sample2$Position <- as.numeric(data_sample2$Position)

  # Detect y-axis labels based on original column names in data_sample2
  baf_label <- if ("baf" %in% names(data_sample2)) "BAF" else if ("ratio" %in% names(data_sample2)) "Ratio" else "BAF"
  z_label   <- if ("z"   %in% names(data_sample2)) "z"   else if ("R"     %in% names(data_sample2)) "R"     else "z"

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
            brks <- c(-Inf, as.numeric(s), Inf)
            cut(as.numeric(Position), breaks = brks, labels = FALSE, right = FALSE)
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
  # Ensure baf/z columns exist (needed when data+geno.pos path was used)
  if (!"baf" %in% names(marker_df) && "ratio" %in% names(marker_df))
    marker_df$baf <- as.numeric(marker_df$ratio)
  if (!"z" %in% names(marker_df) && "R" %in% names(marker_df))
    marker_df$z <- as.numeric(marker_df$R)
  marker_df$Position <- as.numeric(marker_df$Position)
  # Add a column to marker_df to indicate R == 0
  marker_df$R0 <- ifelse(!is.null(marker_df$R) & marker_df$R == 0, TRUE, FALSE)

  # Determine shape aesthetic: priority is outlier+R0, then outlier, then R0, else none
  shape_aes <- NULL
  shape_scale <- NULL
  if (depth_zero_as_x) {
    if ("outlier" %in% colnames(marker_df)) {
      # Combine outlier and R0 into a single shape variable
      marker_df$shape_flag <- ifelse(marker_df$R0, "R0", ifelse(marker_df$outlier, "Outlier", "Normal"))
      shape_aes <- aes(shape = shape_flag)
      shape_scale <- scale_shape_manual(
        values = c("Normal" = 16, "Outlier" = 4, "R0" = 4),
        name = "Marker flag",
        breaks = c("Normal", "Outlier", "R0"),
        labels = c("Normal", "Outlier", "R = 0")
      )
    } else if ("R0" %in% colnames(marker_df)) {
      marker_df$shape_flag <- ifelse(marker_df$R0, "R0", "Normal")
      shape_aes <- aes(shape = shape_flag)
      shape_scale <- scale_shape_manual(
        values = c("Normal" = 16, "R0" = 4),
        name = "Marker flag",
        breaks = c("Normal", "R0"),
        labels = c("Normal", "R = 0")
      )
    }
  } else {
    # No special shape for R0
    if ("outlier" %in% colnames(marker_df)) {
      marker_df$shape_flag <- ifelse(marker_df$outlier, "Outlier", "Normal")
      shape_aes <- aes(shape = shape_flag)
      shape_scale <- scale_shape_manual(
        values = c("Normal" = 16, "Outlier" = 4),
        name = "Marker flag",
        breaks = c("Normal", "Outlier"),
        labels = c("Normal", "Outlier")
      )
    } else {
      shape_aes <- NULL
      shape_scale <- NULL
    }
  }

  if (summarized) {
    # Only plot z means (no points)
    z_means <- hmm_CN$by_window %>%
      filter(Sample == sample_id & !is.na(Chr) & Chr %in% chr_order) %>%
      select(Chr, WindowID, Start, End, w_baf, z_mean = z) %>%
      arrange(factor(Chr, levels = chr_order), Start)
    z_means$Chr   <- factor(z_means$Chr, levels = chr_order)
    z_means$Start <- as.numeric(z_means$Start)
    z_means$End   <- as.numeric(z_means$End)
    p_z <- ggplot() +
      geom_segment(
        data = z_means,
        aes(x = Start, xend = End, y = z_mean, yend = z_mean, group = WindowID),
        color = "black", linewidth = 0.7
      ) +
      geom_blank(
        data = chrom_ghost,
        inherit.aes = FALSE,
        aes(x = Position, y = min(z_means$z_mean, na.rm = TRUE))
      ) +
      facet_wrap(~ Chr, scales = "free_x", nrow = 1) +
      labs(x = NULL, y = z_label) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        strip.background   = element_rect(fill = "grey95"),
        legend.position    = "none",
        axis.title.x       = element_blank(),
        axis.text.x        = element_text(angle = 30, vjust = 1, hjust = 1, size = 8)
      ) +
      scale_x_continuous(labels = scales::scientific_format())
  } else if (z_by_mean) {
    z_means <- hmm_CN$by_window %>%
      filter(Sample == sample_id & !is.na(Chr) & Chr %in% chr_order) %>%
      select(Chr, WindowID, Start, End, w_baf, z_mean = z) %>%
      arrange(factor(Chr, levels = chr_order), Start)
    z_means$Chr   <- factor(z_means$Chr, levels = chr_order)
    z_means$Start <- as.numeric(z_means$Start)
    z_means$End   <- as.numeric(z_means$End)
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
    if (!is.null(shape_scale)) {
      p_z <- p_z + shape_scale
    }
    p_z <- p_z +
      facet_wrap(~ Chr, scales = "free_x", nrow = 1) +
      scale_color_viridis_c(option = "plasma", direction = -1,
                            limits = c(0, 1),
                            oob = squish,
                            name = "BAF weight") +
      labs(x = NULL, y = z_label) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        strip.background   = element_rect(fill = "grey95"),
        legend.position    = "none",
        axis.title.x       = element_blank(),
        axis.text.x        = element_text(angle = 30, vjust = 1, hjust = 1, size = 8)
      ) +
      scale_x_continuous(labels = scales::scientific_format())
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
    if (!is.null(shape_scale)) {
      p_z <- p_z + shape_scale
    }
    p_z <- p_z +
      facet_wrap(~ Chr, scales = "free_x", nrow = 1) +
      scale_color_viridis_c(option = "plasma", direction = -1,
                            limits = c(0, 1),
                            oob = squish,
                            name = "BAF weight") +
      labs(x = NULL, y = z_label) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        strip.background   = element_rect(fill = "grey95"),
        legend.position    = "none",
        axis.title.x       = element_blank(),
        axis.text.x        = element_text(angle = 30, vjust = 1, hjust = 1, size = 8)
      ) +
      scale_x_continuous(labels = scales::scientific_format())
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
      axis.text.x        = element_text(angle = 30, vjust = 1, hjust = 1, size = 8)
    ) +
    scale_x_continuous(labels = scales::scientific_format())

  # -------- BAF panel (color = w_baf; uses SAME tagging/join as z panel) --------
  plot_df <- data_sample2 %>%
    tag_markers_with_region(break_tbl) %>%
    left_join(win_tbl %>% select(Chr, region_id),
                     by = c("Chr", "region_id"))
  # Ensure baf/z columns exist (needed when data+geno.pos path was used)
  if (!"baf" %in% names(plot_df) && "ratio" %in% names(plot_df))
    plot_df$baf <- as.numeric(plot_df$ratio)
  if (!"z" %in% names(plot_df) && "R" %in% names(plot_df))
    plot_df$z <- as.numeric(plot_df$R)
  plot_df$Position <- as.numeric(plot_df$Position)
  dens_list <- lapply(split(plot_df, plot_df$Chr), function(df) {
    if (nrow(df) < 2) return(NULL)
    dens <- stats::density(df$baf, na.rm = TRUE)
    data.frame(
      Chr = unique(df$Chr),
      y = dens$x,
      density = dens$y
    )
  })
  dens_df <- do.call(rbind, dens_list)

  # Add per-chromosome Position range and per-chromosome density max to dens_df
  chrom_pos_range <- plot_df %>%
    group_by(Chr) %>%
    summarise(pos_min = min(Position, na.rm = TRUE), pos_max = max(Position, na.rm = TRUE), .groups = "drop")
  dens_df <- dens_df %>%
    left_join(chrom_pos_range, by = "Chr") %>%
    group_by(Chr) %>%
    mutate(max_dens = max(density, na.rm = TRUE)) %>%
    ungroup()

  if (summarized) {
    # Only plot the BAF polygon (no points)
    p_baf <- ggplot() +
      geom_blank(
        data = chrom_ghost,
        inherit.aes = FALSE,
        aes(x = Position, y = 0)
      ) +
      geom_polygon(
        data = dens_df,
        aes(x = pos_min + (pos_max - pos_min) * density / max_dens,
            y = y, group = Chr),
        fill = "black", alpha = 0.6, color = NA,
        inherit.aes = FALSE
      )
    if (show_window_lines) {
      p_baf <- p_baf + geom_vline(data = vlines, aes(xintercept = Start),
                 color = line_color, alpha = line_alpha,
                 linewidth = line_width, linetype = line_linetype)
    }
    p_baf <- p_baf +
      facet_wrap(~Chr, scales = "free_x", nrow = 1) +
      theme_bw() +
      ylab(baf_label) +
      theme(
        axis.text.x        = element_text(angle = 30, vjust = 1, hjust = 1, size = 8),
        legend.position    = "top",
        text               = element_text(size = 12),
        axis.title.x       = element_blank()
      ) +
    scale_x_continuous(labels = scales::scientific_format())
  } else {
    p_baf <- ggplot(plot_df, aes(x = Position, y = baf, color = w_baf)) +
      geom_blank(
        data = chrom_ghost,
        inherit.aes = FALSE,
        aes(x = Position, y = 0)
      ) +
      # Add density polygon on y-axis (background)
      geom_polygon(
        data = dens_df,
        aes(x = pos_min + (pos_max - pos_min) * density / max_dens,
            y = y, group = Chr),
        fill = "#848d98", alpha = 0.3, color = NA,
        inherit.aes = FALSE
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
      ylab(baf_label) +
      scale_color_viridis_c(option = "plasma", direction = -1,
                            limits = c(0, 1),
                            oob = squish,
                            name = "BAF weight") +
      theme(
        axis.text.x        = element_text(angle = 30, vjust = 1, hjust = 1, size = 8),
        legend.position    = "top",
        text               = element_text(size = 12),
        axis.title.x       = element_blank()
      ) +
    scale_x_continuous(labels = scales::scientific_format())
  }

  # -------- stack with ggpubr (no patchwork needed) --------
  # BAF distribution plot (all markers, all chromosomes/windows)
  p_baf_dist <- ggplot(plot_df, aes(x = baf)) +
    geom_density(fill = "#848d98", alpha = 0.5, color = NA) +
    theme_bw() +
    labs(title = paste(baf_label, "distribution"), x = baf_label, y = "Density") +
    scale_x_continuous(limits = range(plot_df$baf, na.rm = TRUE)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "grey60", fill = NA),
      legend.position = "top"
    )


  # Extract legends from BAF and CN panels (call only on ggplot objects)
  legend_baf <- get_legend(p_baf)
  legend_cn <- get_legend(p_cn)

  # Stack legends vertically
  legends_v <- ggarrange(legend_baf, legend_cn, ncol = 1, nrow = 2, heights = c(1, 1))

  # Prepare text panel for HMM params
  if(!any(grepl("params$",names(hmm_CN)))) {
    params <- hmm_CN$params_samples[[sample_id]]
  } else {
    params <- hmm_CN$params
  }
  cn_grid <- if (!is.null(params$cn_grid)) paste(params$cn_grid, collapse = ", ") else "NA"
  dist <- if (!is.null(params$distribution)) as.character(params$distribution) else "NA"
  sample_name <- if (!is.null(sample_id)) as.character(sample_id) else "NA"
  min_snps <- if (!is.null(params$min_snps_per_window)) round(params$min_snps_per_window, 3) else "NA"
  param_text <- paste0(
    "Sample Name: ", sample_name, "\n",
    "Est. overall ploidy: ", if (!is.null(params$exp_ploidy)) params$exp_ploidy else "NA", "\n",
    "CN grid: ", cn_grid, "\n",
    "Min SNPs p/window: ", min_snps, "\n",
    "Distribution: ", dist, "\n"
  )
  param_panel <- ggplot() +
    theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = param_text, hjust = 0.5, vjust = 0.5, size = 3.5, family = "mono")

  # First row: three columns (BAF distribution, stacked legends, HMM params)
  first_row <- ggarrange(
    p_baf_dist,
    legends_v,
    param_panel,
    ncol = 3, nrow = 1, widths = c(3, 1, 2)
  )

  # Stack all panels vertically
  arranged <- ggarrange(
    first_row,
    p_baf + theme(legend.position = "none"),
    p_z + theme(legend.position = "none"),
    p_cn + theme(legend.position = "none"),
    ncol = 1, nrow = 4,
    heights = c(2, heights[1], heights[2], heights[2]),
    align = "v"
  )

  return(list(
    baf = p_baf,
    z   = p_z,
    cn  = p_cn,
    baf_dist = p_baf_dist,
    arranged = arranged
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
##' @param facet_nrow Number of rows for facet_wrap. If set, facet_ncol will be determined automatically unless both are set.
##' @param gray_CN Integer or NULL. If provided, this copy-number value is colored gray (used as the baseline color). All CNs below it are colored blue and above red. If NULL (default), the baseline is auto-detected as the most frequent CN weighted by window length.
#' @param add_het Logical. If TRUE, a heterozygosity column is added to the right side of the plot, with one square per sample colored from blue (low) to red (high). Requires \code{hmm_dosage_calls}. Default is TRUE.
#' @param hmm_dosage_calls An object of class \code{hmm_dosage_calls} (inherits \code{data.frame}) with columns \code{SampleName}, \code{Chr}, \code{dosage}, and optionally others. Used to compute per-sample heterozygosity as the fraction of markers with dosage not equal to 0 or 1 (among non-NA dosage values). Required when \code{add_het = TRUE}.
#' @param interactive Logical. If TRUE, returns an interactive \code{plotly} figure instead of a static ggplot. Hovering over segments shows: Sample, Chr, Start, End, CN, post_max, w_baf, and Heterozygosity (if \code{add_het = TRUE} and \code{hmm_dosage_calls} is provided). Requires the \pkg{plotly} package. Default is FALSE.
#'
#' @return A ggplot/ggarrange object (static), or a plotly figure when \code{interactive = TRUE}.
#'
#' @importFrom dplyr filter group_by summarise arrange desc mutate left_join
#' @importFrom ggplot2 ggplot geom_segment geom_tile geom_point facet_wrap
#' @importFrom ggplot2 scale_x_continuous scale_color_manual scale_alpha scale_fill_gradient
#' @importFrom ggplot2 labs theme_bw theme aes
#' @importFrom ggplot2 element_blank element_rect element_text
#' @importFrom ggpubr ggarrange get_legend
#'
#' @export
compare_cn_track <- function(hmm_CN,
                             samples_to_plot = NULL,
                             chromosomes = NULL,
                             facet_ncol = NULL,
                             facet_nrow = NULL,
                             gray_CN = NULL,
                             add_het = TRUE,
                             hmm_dosage_calls = NULL,
                             interactive = FALSE) {

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

  # baseline CN = user-supplied gray_CN, or most frequent CN weighted by window length
  if (!is.null(gray_CN)) {
    baseline_cn <- as.integer(gray_CN)
  } else {
    plot_df$w <- pmax(0, plot_df$End - plot_df$Start)
    cn_totals <- plot_df |>
      group_by(.data$CN_call) |>
      summarise(total_bp = sum(.data$w, na.rm = TRUE), .groups = "drop") |>
      arrange(desc(.data$total_bp))
    baseline_cn <- cn_totals$CN_call[1]
  }

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
  p_main <- ggplot(plot_df) +
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
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  # ---- compute het_df (shared by both interactive and static paths) ----
  het_df <- NULL
  if (add_het && !is.null(hmm_dosage_calls)) {
    if (!inherits(hmm_dosage_calls, "hmm_dosage_calls")) {
      stop("hmm_dosage_calls must be an object of class 'hmm_dosage_calls'")
    }
    het_raw <- hmm_dosage_calls |>
      filter(.data$SampleName %in% samples_present)
    if (!is.null(chromosomes)) {
      het_raw <- het_raw |> filter(.data$Chr %in% chromosomes)
    }
    het_df <- het_raw |>
      group_by(.data$SampleName) |>
      summarise(
        het = sum(!is.na(.data$dosage) & .data$dosage != 0 & .data$dosage != 1,
                  na.rm = TRUE) / max(sum(!is.na(.data$dosage)), 1L),
        .groups = "drop"
      )
    het_df$SampleName <- factor(het_df$SampleName, levels = levels(plot_df$Sample))
    het_df$label <- "Het."
  }

  # ---- interactive (plotly) output ----
  if (interactive) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop("Package 'plotly' is required for interactive = TRUE. Install with: install.packages('plotly')")
    }

    # Build per-segment tooltip lookup for het
    het_lookup <- if (!is.null(het_df)) {
      setNames(round(het_df$het, 3), as.character(het_df$SampleName))
    } else NULL

    plot_df$Mid <- (as.numeric(plot_df$Start) + as.numeric(plot_df$End)) / 2
    plot_df$tooltip_text <- paste0(
      "Sample: ", plot_df$Sample,
      "<br>Chr: ", plot_df$Chr,
      "<br>Start: ", format(as.numeric(plot_df$Start), big.mark = ",", scientific = FALSE),
      "<br>End: ", format(as.numeric(plot_df$End), big.mark = ",", scientific = FALSE),
      "<br>CN: ", plot_df$CN_call,
      "<br>post_max: ", round(plot_df$post_max, 3),
      "<br>w_baf: ", if ("w_baf" %in% names(plot_df)) round(plot_df$w_baf, 3) else "NA",
      "<br>w_het: ", round(plot_df$n_het/plot_df$n_snps, 3)
    )

    # map CN colors to hex per row
    plot_df$seg_color <- cn_pal[as.character(plot_df$CN_call)]

    chrs <- levels(plot_df$Chr)[levels(plot_df$Chr) %in% unique(as.character(plot_df$Chr))]
    n_chr <- length(chrs)

    # one subplot column per chromosome.
    # Draw one add_segments trace per CN level per chromosome so every CN level
    # gets a proper legend entry (shown only on the first chr where that CN appears).
    # Use only observed CN levels to avoid NA in showlegend.
    cn_levels_all <- as.character(sort(unique(plot_df$CN_call)))

    # For each CN level, find the first chromosome where it has data
    first_chr_for_cn <- vapply(cn_levels_all, function(cn) {
      found <- Filter(function(ch) {
        any(as.character(plot_df$CN_call[plot_df$Chr == ch]) == cn)
      }, chrs)
      if (length(found) == 0L) chrs[1] else found[1]
    }, character(1))

    chr_figs <- lapply(seq_along(chrs), function(i) {
      ch <- chrs[i]
      d <- plot_df[plot_df$Chr == ch, ]
      fig <- plotly::plot_ly()

      # One trace per CN level — guarantees a legend swatch for every level
      for (cn in cn_levels_all) {
        dsub <- d[as.character(d$CN_call) == as.character(cn), ]
        if (nrow(dsub) == 0) next
        fig <- fig |>
          plotly::add_segments(
            x    = as.numeric(dsub$Start),
            xend = as.numeric(dsub$End),
            y    = as.character(dsub$Sample),
            yend = as.character(dsub$Sample),
            line = list(color = cn_pal[[as.character(cn)]], width = 10),
            name = as.character(cn),
            legendgroup = as.character(cn),
            showlegend = (ch == first_chr_for_cn[[as.character(cn)]]),
            opacity = 1,
            hoverinfo = "none"
          )
      }

      # Invisible tooltip anchors at segment midpoints
      sample_order <- rev(levels(plot_df$Sample))
      fig <- fig |>
        plotly::add_markers(
          x    = d$Mid,
          y    = as.character(d$Sample),
          text = d$tooltip_text,
          hoverinfo = "text",
          marker = list(color = "rgba(0,0,0,0)", size = 6, line = list(width = 0)),
          showlegend = FALSE
        ) |>
        plotly::layout(
          xaxis = list(title = ch, tickformat = ".2s", tickangle = 270),
          yaxis = list(
            title = "",
            categoryarray = sample_order,
            categoryorder = "array",
            showticklabels = (i == 1)
          )
        )
      fig
    })

    fig_main <- plotly::subplot(chr_figs,
                                nrows = 1,
                                shareY = TRUE,
                                margin = 0.002)

    if (!is.null(het_df)) {
      het_df$tooltip_text <- paste0(
        "Sample: ", het_df$SampleName,
        "<br>Heterozygosity: ", round(het_df$het, 3)
      )

      sample_order <- rev(levels(plot_df$Sample))
      fig_het <- plotly::plot_ly() |>
        plotly::add_markers(
          x    = het_df$label,
          y    = as.character(het_df$SampleName),
          text = het_df$tooltip_text,
          hoverinfo = "text",
          marker = list(
            symbol    = "square",
            size      = 20,
            color     = het_df$het,
            colorscale = list(list(0, "#2166ac"), list(1, "#d6604d")),
            cmin      = 0,
            cmax      = 1,
            showscale = TRUE,
            colorbar  = list(title = "Het.", len = 0.4, y = 0.5),
            line      = list(width = 0)
          ),
          showlegend = FALSE
        ) |>
        plotly::layout(
          xaxis = list(title = ""),
          yaxis = list(
            title = "",
            categoryarray = sample_order,
            categoryorder = "array",
            showticklabels = FALSE
          )
        )
      return(plotly::subplot(fig_main, fig_het,
                             nrows = 1, widths = c(0.92, 0.08),
                             shareY = FALSE, margin = 0.0001))
    }
    return(fig_main)
  }

  # ---- static (ggplot) output ----
  if (!is.null(het_df)) {
    het_df$facet_dummy <- " "  # blank strip to match p_main chromosome strip height

    p_het <- ggplot(het_df, aes(x = .data$label, y = .data$SampleName, fill = .data$het)) +
      geom_tile(color = "white") +
      scale_fill_gradient(
        low = "#2166ac", high = "#d6604d",
        name = "Heterozygosity",
        limits = c(0, 1)
      ) +
      facet_wrap(~ facet_dummy) +
      theme_bw() +
      labs(x = "", y = NULL) +
      theme(
        axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        panel.grid       = element_blank(),
        legend.position  = "top",
        axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.background = element_rect(fill = "grey95"),
        strip.text       = element_text(colour = "grey95")  # invisible text, keeps height
      )

    p_main <- p_main + theme(legend.position = "top")

    leg_main <- get_legend(p_main)
    leg_het  <- get_legend(p_het)

    plots_row <- ggarrange(
      p_main + theme(legend.position = "none"),
      p_het  + theme(legend.position = "none"),
      ncol = 2, widths = c(15, 1), align = "hv"
    )
    legends_row <- ggarrange(leg_main, leg_het, ncol = 2, widths = c(7, 4))
    return(ggarrange(legends_row, plots_row, nrow = 2, heights = c(2, 20)))
  }

  return(p_main)
}
