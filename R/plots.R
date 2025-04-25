globalVariables(c("Chr", "baf", "z", "ratio", "median", "window", "prop_het"))

#' Plot BAF
#'
#' This function generates a BAF (B-allele frequency) plot for visualizing
#' genomic data. It allows customization of dot size, expected and estimated
#' peaks, centromere positions, and area colors.
#'
#' @param data_sample A data.frame containing BAF and genomic position
#' information. Must include columns `Chr`, `Position`, and `sample`.
#' @param area_single Numeric value defining the area around the expected peak
#' to be considered.
#' @param ploidy Integer or vector specifying the expected ploidy. If a vector,
#' it must match the number of chromosomes in `data_sample`.
#' @param dot.size Numeric value for the size of the dots in the plot. Default
#' is 1.
#' @param add_estimated_peaks Logical. If TRUE, adds lines for estimated peaks.
#' Default is FALSE.
#' @param add_expected_peaks Logical. If TRUE, adds lines for expected peaks.
#' Default is FALSE.
#' @param colors Logical. If TRUE, adds area colors to the plot. Default is
#' FALSE.
#' @param centromeres Named vector defining centromere positions for each
#' chromosome. Names must match chromosome IDs in `data_sample`.
#' @param add_centromeres Logical. If TRUE, adds vertical lines at centromere
#' positions. Default is FALSE.
#' @param font_size Numeric value for the font size of plot labels. Default is
#' 12.
#' @return A ggplot object representing the BAF plot.
#'
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom dplyr %>%
#' @export
plot_baf <- function(data_sample,
                     area_single,
                     ploidy,
                     dot.size = 1,
                     add_estimated_peaks = FALSE,
                     add_expected_peaks = FALSE,
                     centromeres = NULL,
                     add_centromeres = FALSE,
                     colors = FALSE,
                     font_size = 12) {
  if (add_centromeres) {
    # centromeres <- c("1 = 49130338, 5 = 49834357")
    if (length(centromeres) == 1 && any(grepl(",", centromeres))) {
      centromeres <- gsub("\ ", "", centromeres)
      centromeres <- unlist(strsplit(centromeres, ","))
      centromeres <- sapply(centromeres, function(x) strsplit(x, "="))
      centromeres_df <- data.frame(
        Chr = sapply(centromeres, "[[", 1),
        value = sapply(centromeres, "[[", 2)
      )
    } else {
      centromeres_df <- data.frame(Chr = names(centromeres), value = centromeres)
    }
  }

  if (add_estimated_peaks || add_expected_peaks || colors) {
    if (length(ploidy) == 1) {
      ploidy <- rep(ploidy, length(unique(data_sample$Chr)))
    } else if (length(ploidy) > 1 &&
      length(ploidy) < length(unique(data_sample$Chr))) {
      stop("Provide a ploidy for each chromosome.")
    }

    data_sample2 <- rets2 <- data.frame()
    for (z in seq_along(unique(data_sample$Chr))) {
      ymin <- seq(0, 1, 1 / ploidy[z]) - (area_single / (ploidy[z] * 2))
      ymax <- seq(0, 1, 1 / ploidy[z]) + (area_single / (ploidy[z] * 2))

      ymin[which(ymin < 0)] <- 0
      ymax[which(ymax > 1)] <- 1
      rets <- data.frame(ymin, ymax,
        xmax = Inf, xmin = -Inf,
        Chr = sort(unique(data_sample$Chr))[z]
      )

      data_chr <- data_sample[which(data_sample$Chr ==
        sort(unique(data_sample$Chr))[z]), ]
      idx_tot <- FALSE
      idx <- list()
      for (i in seq_len(nrow(rets))) {
        idx <- data_chr$sample >= rets$ymin[i] & data_chr$sample <= rets$ymax[i]
        idx_tot <- idx_tot | idx
      }

      data_chr$color <- NA
      data_chr$color[which(!idx_tot)] <- "red"
      data_chr$color[which(idx_tot)] <- "black"

      rets2 <- rbind(rets2, rets)
      data_sample2 <- rbind(data_sample2, data_chr)
    }
  } else {
    data_sample2 <- data_sample
  }

  p_baf <- data_sample2 %>% ggplot(aes(x = Position, y = sample)) +
    {
      if (colors) {
        geom_point(aes(color = color), alpha = 0.7, size = dot.size)
      } else {
        geom_point(alpha = 0.7, size = dot.size)
      }
    } +
    scale_color_manual(values = c("red", "black")) +
    facet_grid(~Chr, scales = "free_x") +
    theme_bw() +
    {
      if (add_expected_peaks) {
        geom_rect(
          data = rets2, inherit.aes = FALSE,
          aes(
            ymin = ymin, ymax = ymax,
            xmax = xmax, xmin = xmin,
            alpha = 0.001
          ),
          color = "red"
        )
      }
    } +
    ylab("BAF") +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "none", text = element_text(size = font_size)
    ) +
    {
      if (add_centromeres) {
        geom_vline(
          data = centromeres_df,
          aes(xintercept = as.numeric(value)),
          color = "blue",
          linewidth = 0.8
        )
      }
    }
  return(p_baf)
}

#' Plot BAF Histogram
#'
#' This function generates a histogram of BAF (B-allele frequency) values. It
#' supports options for adding estimated and expected peaks, area colors, and
#' filtering homozygous calls.
#'
#' @param data_sample A data.frame containing BAF and genomic position
#' information. Must include columns `Chr`, `Position`, and `sample`.
#' @param area_single Numeric value defining the area around the expected peak
#' to be considered.
#' @param ploidy Integer or vector specifying the expected ploidy. If a vector,
#' it must match the number of chromosomes in `data_sample`.
#' @param colors Logical. If TRUE, adds area colors to the histogram. Default
#' is FALSE.
#' @param add_estimated_peaks Logical. If TRUE, adds lines for estimated peaks.
#' Default is TRUE.
#' @param add_expected_peaks Logical. If TRUE, adds lines for expected peaks.
#' Default is FALSE.
#' @param BAF_hist_overall Logical. If TRUE, plots the BAF histogram for the
#' entire genome. Default is FALSE.
#' @param ratio Logical. If TRUE, plots the raw ratio instead of BAF. Default
#' is FALSE.
#' @param rm_homozygous Logical. If TRUE, removes homozygous calls from the
#' histogram. Default is FALSE.
#' @param font_size Numeric value for the font size of plot labels. Default is
#' 12.
#' @return A ggplot object representing the BAF histogram.
#'
#' @import ggplot2
#' @importFrom dplyr filter
#' @export
plot_baf_hist <- function(data_sample,
                          area_single,
                          ploidy,
                          colors = FALSE,
                          add_estimated_peaks = TRUE,
                          add_expected_peaks = FALSE,
                          BAF_hist_overall = FALSE,
                          ratio = FALSE,
                          rm_homozygous = FALSE,
                          font_size = 12) {
  if ((add_estimated_peaks || add_expected_peaks || colors) &&
    !all(is.na(ploidy))) {
    if (length(ploidy) == 1) {
      ploidy <- rep(ploidy, length(unique(data_sample$Chr)))
    } else if (length(ploidy) != 1 &&
      length(ploidy) != length(unique(data_sample$Chr))) {
      stop("Provide a ploidy for each chromosome.")
    }

    if (BAF_hist_overall) {
      Chrs <- "all"
    } else {
      Chrs <- sort(unique(data_sample$Chr))
    }
    data_sample2 <- modes.df3 <- data.frame()
    for (z in seq_along(Chrs)) {
      ymin <- seq(0, 1, 1 / ploidy[z]) - (area_single / (ploidy[z] * 2))
      ymax <- seq(0, 1, 1 / ploidy[z]) + (area_single / (ploidy[z] * 2))

      ymin[which(ymin < 0)] <- 0
      ymax[which(ymax > 1)] <- 1
      rets <- data.frame(ymin, ymax, xmax = Inf, xmin = -Inf, Chr = Chrs[z])

      if (all(Chrs == "all")) {
        split_data <- data_sample
      } else {
        split_data <- data_sample[which(data_sample$Chr ==
          sort(unique(data_sample$Chr))[z]), ]
      }

      if (ratio) split_data$sample <- split_data$ratio

      idx_tot <- FALSE
      modes.df2 <- data.frame()
      for (i in seq_len(nrow(rets))) {
        idx <- split_data$sample >= rets$ymin[i] & split_data$sample <= rets$ymax[i]
        idx_tot <- idx_tot | idx

        estimated <- mode(split_data$sample[which(split_data$sample >
          rets$ymin[i] &
          split_data$sample <
            rets$ymax[i])])
        expected <- seq(0, 1, 1 / ploidy[z])[i]
        modes.df <- cbind(
          Chr = Chrs[z],
          pivot_longer(data.frame(estimated, expected),
            cols = 1:2
          )
        )
        modes.df2 <- rbind(modes.df2, modes.df)
      }

      modes.df2$alpha <- modes.df2$name
      modes.df2$alpha <- gsub("expected", 1, modes.df2$alpha)
      modes.df2$alpha <- as.numeric(gsub("estimated", 0.5, modes.df2$alpha))

      split_data$color <- NA
      split_data$color[which(!idx_tot)] <- "outside area"
      split_data$color[which(idx_tot)] <- "inside area"

      modes.df3 <- rbind(modes.df2, modes.df3)
      data_sample2 <- rbind(data_sample2, split_data)
    }
  } else {
    data_sample2 <- data_sample
    add_estimated_peaks <- add_expected_peaks <- colors <- FALSE
  }

  if (rm_homozygous) data_sample2 <- data_sample2 %>% filter(sample != 0 & sample != 1)
  if (ratio) data_sample2$sample <- data_sample2$ratio

  if (!BAF_hist_overall) {
    p_hist <- data_sample2 %>% ggplot(aes(x = sample)) +
      {
        if (colors) geom_histogram(aes(fill = color)) else geom_histogram()
      } +
      # scale_x_continuous(breaks = round(seq(0, 1, 1/ploidy),2)) +
      {
        if (add_estimated_peaks) {
          geom_vline(
            data = modes.df3[which(modes.df3$name == "estimated"), ],
            aes(
              xintercept = value,
              color = name,
              linetype = name,
              alpha = alpha
            ),
            linewidth = 0.8
          )
        }
      } +
      {
        if (add_expected_peaks) {
          geom_vline(
            data = modes.df3[which(modes.df3$name == "expected"), ],
            aes(
              xintercept = value,
              color = name,
              linetype = name,
              alpha = alpha
            ),
            linewidth = 0.8
          )
        }
      } +
      {
        if (colors) scale_fill_manual(values = c("red", "black"))
      } +
      {
        if (add_expected_peaks || add_estimated_peaks) {
          scale_color_manual(values = c("blue", "purple"))
        }
      } +
      scale_linetype_manual(values = c("dashed", "solid"), guide = "none") +
      scale_alpha(range = c(0.7, 1), guide = "none") +
      facet_grid(~Chr, scales = "free_x") +
      theme_bw() +
      xlab("BAF") +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom", text = element_text(size = font_size)
      ) +
      {
        if (add_expected_peaks || add_estimated_peaks) labs(color = "Peaks")
      } +
      {
        if (colors) labs(fill = "Area")
      }
  } else {
    if (rm_homozygous) data_sample <- data_sample %>% filter(sample != 0 & sample != 1)
    if (ratio) data_sample$sample <- data_sample$ratio

    p_hist_all <- data_sample %>% ggplot(aes(x = sample)) +
      geom_histogram() +
      {
        if (colors) geom_histogram(aes(fill = color)) else geom_histogram()
      } +
      # scale_x_continuous(breaks = round(seq(0, 1, 1/ploidy),2)) +
      {
        if (add_estimated_peaks) {
          geom_vline(
            data = modes.df3[which(modes.df3$name == "estimated"), ],
            aes(
              xintercept = value,
              color = name,
              linetype = name,
              alpha = alpha
            ),
            linewidth = 0.8
          )
        }
      } +
      {
        if (add_expected_peaks) {
          geom_vline(
            data = modes.df3[which(modes.df3$name == "expected"), ],
            aes(
              xintercept = value,
              color = name,
              linetype = name,
              alpha = alpha
            ),
            linewidth = 0.8
          )
        }
      } +
      {
        if (colors) scale_fill_manual(values = c("red", "black"))
      } +
      {
        if (add_expected_peaks || add_estimated_peaks) {
          scale_color_manual(values = c("blue", "purple"))
        }
      } +
      scale_linetype_manual(values = c("dashed", "solid"), guide = "none") +
      scale_alpha(range = c(0.7, 1), guide = "none") +
      theme_bw() +
      {
        if (ratio) xlab("ratio") else xlab("BAF")
      } +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom", text = element_text(size = font_size)
      ) +
      {
        if (add_expected_peaks || add_estimated_peaks) labs(color = "Peaks")
      } +
      {
        if (colors) labs(fill = "Area")
      }

    p_hist <- p_hist_all
  }

  return(p_hist)
}

#' Plot Method for Qploidy Standardization
#'
#' This function generates various plots for visualizing the results of Qploidy
#' standardization. It supports multiple plot types, including BAF, z-score,
#' and histograms.
#'
#' @param x An object of class `qploidy_standardization`.
#' @param sample Character string indicating the sample ID to plot.
#' @param chr Character or numeric vector specifying the chromosomes to plot.
#' Default is NULL (plots all chromosomes).
#' @param type Character vector defining the plot types. Options include:
#'   - "all": Generates all available plot types.
#'   - "het": Plots heterozygous locus counts across genomic windows.
#'   - "BAF": Plots B-allele frequency (BAF) for each chromosome.
#'   - "zscore": Plots z-scores for each chromosome.
#'   - "BAF_hist": Plots BAF histograms for each chromosome.
#'   - "BAF_hist_overall": Plots a BAF histogram for the entire genome.
#'   - "Ratio_hist_overall": Plots a histogram of raw ratios for the entire
#'   genome.
#'   - "ratio": Plots raw ratios for each chromosome.
#'   Default is "all".
#' @param area_single Numeric value defining the area around the expected peak
#' to be considered. Default is 0.75.
#' @param ploidy Integer specifying the expected ploidy. Default is 4.
#' @param dot.size Numeric value for the size of the dots in the plots. Default
#' is 1.
#' @param font_size Numeric value for the font size of plot labels. Default is
#' 12.
#' @param add_estimated_peaks Logical. If TRUE, adds lines for estimated peaks.
#' Default is FALSE.
#' @param add_expected_peaks Logical. If TRUE, adds lines for expected peaks.
#' Default is FALSE.
#' @param centromeres Named vector defining centromere positions for each
#' chromosome. Names must match chromosome IDs in `x`.
#' @param add_centromeres Logical. If TRUE, adds vertical lines at centromere
#' positions. Default is FALSE.
#' @param colors Logical. If TRUE, adds area colors to the plots. Default is
#' FALSE.
#' @param window_size Numeric value defining the genomic position window for
#' heterozygous locus counts. Default is 2000000.
#' @param het_interval Numeric value defining the interval to consider as
#' heterozygous. Default is 0.1.
#' @param rm_homozygous Logical. If TRUE, removes homozygous calls from BAF
#' histogram plots. Default is FALSE.
#' @param ... Additional plot parameters.
#' @return A ggarrange object containing the requested plots.
#'
#' @details The function supports the following plot types:
#'
#' - **all**: Generates all available plot types.
#' - **het**: Plots the proportion of heterozygous loci across genomic windows,
#' useful for identifying regions with high or low heterozygosity.
#' - **BAF**: Plots the B-allele frequency (BAF) for each chromosome, showing
#' the distribution of allele frequencies.
#' - **zscore**: Plots z-scores for each chromosome, which can help identify
#' outliers or regions with unusual data distributions.
#' - **BAF_hist**: Plots histograms of BAF values for each chromosome,
#' providing a summary of allele frequency distributions.
#' - **BAF_hist_overall**: Plots a single histogram of BAF values for the
#' entire genome, summarizing allele frequency distributions genome-wide.
#' - **Ratio_hist_overall**: Plots a histogram of raw ratios for the entire
#' genome, useful for visualizing overall ratio distributions.
#' - **ratio**: Plots raw ratios for each chromosome, showing the distribution
#' of observed ratios.
#'
#'
#' @importFrom ggpubr ggarrange annotate_figure text_grob
#' @import ggplot2
#' @importFrom dplyr filter group_by summarise
#' @importFrom tidyr pivot_wider
#' @export
plot_qploidy_standardization <- function(x,
                                         sample = NULL,
                                         chr = NULL,
                                         type = c(
                                           "all", "het", "BAF", "zscore",
                                           "BAF_hist", "ratio",
                                           "BAF_hist_overall",
                                           "Ratio_hist_overall"
                                         ),
                                         area_single = 0.75,
                                         ploidy = 4,
                                         dot.size = 1,
                                         font_size = 12,
                                         add_estimated_peaks = FALSE,
                                         add_expected_peaks = FALSE,
                                         centromeres = NULL,
                                         add_centromeres = FALSE,
                                         colors = FALSE,
                                         window_size = 2000000,
                                         het_interval = 0.1,
                                         rm_homozygous = FALSE, ...) {
  if (!inherits(x, "qploidy_standardization")) {
    stop("Object is not of class qploidy_standardization")
  }

  if (is.null(sample)) stop("Define sample ID")

  data_sample <- x$data[which(x$data$SampleName == sample), ]

  if (is.null(chr)) {
    chr <- sort(unique(data_sample$Chr))
  } else if (is.numeric(chr)) {
    chr <- sort(unique(data_sample$Chr))[chr]
  }

  if(!all(chr %in% unique(data_sample$Chr))) {
    stop("Chromosome not found in data.")
  }

  data_sample <- data_sample %>%
    filter(Chr %in% chr) %>%
    select(MarkerName, SampleName, Chr, Position, baf, z, ratio)

  data_sample$Position <- as.numeric(data_sample$Position)

  if (nrow(data_sample) == 0) stop("Sample or chromosome not found.")

  baf_sample <- data_sample %>% pivot_wider(
    names_from = SampleName,
    values_from = baf
  )
  zscore_sample <- data_sample %>% pivot_wider(
    names_from = SampleName,
    values_from = z
  )

  colnames(baf_sample)[ncol(baf_sample)] <- "sample"

  baf_point <- baf_hist <- p_z <- raw_ratio <- het_rate <- baf_hist_overall <-
    ratio_hist_overall <- NULL

  if (any(type == "all" | type == "BAF")) {
    baf_point <- plot_baf(baf_sample,
      area_single,
      ploidy,
      dot.size = dot.size,
      add_estimated_peaks,
      add_expected_peaks,
      centromeres,
      add_centromeres,
      colors,
      font_size = font_size
    )
  }

  if (add_centromeres) {
    # centromeres <- c("1 = 49130338, 5 = 49834357")
    if (length(centromeres) == 1 && any(grepl(",", centromeres))) {
      centromeres <- gsub("\ ", "", centromeres)
      centromeres <- unlist(strsplit(centromeres, ","))
      centromeres <- sapply(centromeres, function(x) strsplit(x, "="))
      centromeres_df <- data.frame(
        Chr = sapply(centromeres, "[[", 1),
        value = sapply(centromeres, "[[", 2)
      )
    } else {
      centromeres_df <- data.frame(Chr = names(centromeres), value = centromeres)
    }

    for (i in seq_along(centromeres)) {
      idx <- which(baf_sample$Chr == centromeres_df$Chr[i] &
        baf_sample$Position < centromeres_df$value[i])
      baf_sample$Chr[idx] <- paste0(baf_sample$Chr[idx], ".1")
      idx <- which(baf_sample$Chr == centromeres_df$Chr[i] &
        baf_sample$Position >= centromeres_df$value[i])
      baf_sample$Chr[idx] <- paste0(baf_sample$Chr[idx], ".2")
      if (any(type == "zscore")) {
        idx <- which(zscore_sample$Chr == centromeres_df$Chr[i] &
          zscore_sample$Position < centromeres_df$value[i])
        zscore_sample$Chr[idx] <- paste0(zscore_sample$Chr[idx], ".1")
        idx <- which(zscore_sample$Chr == centromeres_df$Chr[i] &
          zscore_sample$Position >= centromeres_df$value[i])
        zscore_sample$Chr[idx] <- paste0(zscore_sample$Chr[idx], ".2")
      }
    }
  }

  if (any(type == "all" | type == "BAF_hist")) {
    baf_hist <- plot_baf_hist(
      data_sample = baf_sample,
      area_single,
      ploidy,
      colors,
      add_estimated_peaks,
      add_expected_peaks,
      BAF_hist_overall = FALSE,
      rm_homozygous = rm_homozygous,
      font_size = font_size
    )
  }

  if (any(type == "all" | type == "BAF_hist_overall")) {
    baf_hist_overall <- plot_baf_hist(
      data_sample = baf_sample,
      area_single,
      ploidy,
      colors,
      add_estimated_peaks,
      add_expected_peaks,
      BAF_hist_overall = TRUE,
      rm_homozygous = rm_homozygous,
      font_size = font_size
    )
  }

  if (any(type == "all" | type == "Ratio_hist_overall")) {
    ratio_hist_overall <- plot_baf_hist(
      data_sample = baf_sample,
      area_single,
      ploidy,
      colors,
      add_estimated_peaks,
      add_expected_peaks,
      BAF_hist_overall = TRUE,
      ratio = TRUE,
      rm_homozygous = rm_homozygous,
      font_size = font_size
    )
  }

  if (any(type == "zscore")) {
    colnames(zscore_sample)[ncol(zscore_sample)] <- "z"
    p_z <- zscore_sample %>%
      ggplot(aes(x = Position, y = z)) +
      facet_grid(. ~ Chr, scales = "free") +
      geom_smooth(method = "gam") +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        text = element_text(size = font_size)
      ) +
      geom_hline(yintercept = median(zscore_sample$z), linetype = "dashed")
  }

  if (any(type == "ratio")) {
    raw_ratio <- data_sample %>% ggplot(aes(x = Position, y = ratio)) +
      geom_point(alpha = 0.7, size = dot.size) +
      facet_grid(. ~ Chr, scales = "free") +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        text = element_text(size = font_size)
      )
  }

  if (any(type == "het")) {
    data_sample$het_rate <- NA
    data_sample$het_rate[which(data_sample$ratio <= 0 + het_interval |
      data_sample$ratio >= 1 - het_interval)] <- 0
    data_sample$het_rate[which(data_sample$ratio > 0 + het_interval &
      data_sample$ratio < 1 - het_interval)] <- 1

    data_sample_lst <- split.data.frame(data_sample, data_sample$Chr)

    for (j in seq_along(data_sample_lst)) {
      data_sample_lst[[j]]$window <- 1
      if (length(unique(data_sample_lst[[j]]$Position)) == 1) next()
      intervals_start <- c(seq(
        1, max(data_sample_lst[[j]]$Position),
        window_size
      ))
      intervals_end <- c(seq(
        1, max(data_sample_lst[[j]]$Position),
        window_size
      ) - 1, max(data_sample_lst[[j]]$Position))
      intervals_end <- intervals_end[-1]
      for (i in seq_along(intervals_start)) {
        data_sample_lst[[j]]$window[which(data_sample_lst[[j]]$Position >=
          intervals_start[i] &
          data_sample_lst[[j]]$Position <=
            intervals_end[i])] <- intervals_start[i]
      }
    }

    data_sample <- do.call(rbind, data_sample_lst)

    het_rate <- data_sample %>%
      group_by(Chr, window) %>%
      summarise(prop_het = sum(het_rate, na.rm = TRUE) /
        length(which(het_rate == 1 | het_rate == 0))) %>%
      ggplot(aes(x = window, y = prop_het, color = prop_het)) +
      geom_line() +
      scale_color_gradient(low = "black", high = "red") +
      facet_grid(. ~ Chr, scales = "free") +
      theme_bw() +
      ylab("proportion of heterozygous loci") +
      xlab("Position") +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none", text = element_text(size = font_size)
      )
  }

  p_all <- list(
    het_rate, raw_ratio, baf_point, baf_hist, baf_hist_overall,
    ratio_hist_overall, p_z
  )

  rm <- which(sapply(p_all, is.null))
  if (length(rm) != 0) p_all <- p_all[-rm]

  p_result <- ggarrange(plotlist = p_all, ncol = 1)

  p_result <- annotate_figure(p_result, top = text_grob(sample,
    face = "bold",
    size = 14
  ))

  return(p_result)
}

#' Plot graphics for ploidy visual inspection for each resolution
#'
#' This function generates and saves plots for visual inspection of ploidy at
#' different resolutions: chromosome, chromosome-arm, and sample levels. It is
#' designed for parallelization purposes and supports customization of
#' centromere positions and chromosome selection.
#'
#' @param data_standardized An object of class `qploidy_standardization`
#' containing standardized data for ploidy analysis.
#' @param sample A character string specifying the sample name to be analyzed.
#' @param ploidy A numeric value indicating the expected ploidy of the sample.
#' This parameter is required.
#' @param centromeres A named vector with centromere positions (in base pairs)
#' for each chromosome. The names must match the chromosome IDs in the dataset.
#' This is used for chromosome-arm level resolution.
#' @param file_name A character string defining the output file path and name
#' prefix for the saved plots. The function appends resolution-specific
#' suffixes to this prefix. If NULL, plots are not saved to files.
#' @param chr A vector specifying the chromosomes to include in the analysis.
#' If NULL, all chromosomes are included.
#' @param types_chromosome A character vector defining the plot types for
#' chromosome-level resolution. Options include:
#'   - "het": Plots heterozygous locus counts.
#'   - "BAF": Plots B-allele frequency (BAF).
#'   - "zscore": Plots z-scores.
#'   - "BAF_hist": Plots BAF histograms for each chromosome.
#'   - "ratio": Plots raw ratios for each chromosome.
#'   Default is c("Ratio_hist", "BAF_hist", "zscore").
#' @param types_chromosome_arm A character vector defining the plot types for
#' chromosome-arm level resolution. Options include:
#'   - "het": Plots heterozygous locus counts.
#'   - "BAF": Plots B-allele frequency (BAF).
#'   - "zscore": Plots z-scores.
#'   - "BAF_hist": Plots BAF histograms for each chromosome arm.
#'   - "ratio": Plots raw ratios for each chromosome arm.
#'   Default is c("Ratio_hist", "BAF_hist", "zscore").
#' @param types_sample A character vector defining the plot types for
#' sample-level resolution. Options include:
#'   - "Ratio_hist_overall": Plots a histogram of raw ratios for the entire
#'   genome.
#'   - "BAF_hist_overall": Plots a BAF histogram for the entire genome.
#'   Default is c("Ratio_hist_overall", "BAF_hist_overall").
#' @return A list containing the generated plots for each resolution:
#' - `chromosome`: Plot for chromosome-level resolution.
#' - `chromosome_arm`: Plot for chromosome-arm level resolution (if centromeres
#' are provided).
#' - `sample`: Plot for sample-level resolution.
#' @details The function generates three types of plots:
#'
#' - **Chromosome-level resolution**: Plots raw ratio, BAF histograms,
#' z-scores, heterozygous locus counts, and BAF for each chromosome.
#' - **Chromosome-arm level resolution**: Similar to chromosome-level but
#' splits data by chromosome arms using centromere positions.
#' - **Sample-level resolution**: Combines all markers in the sample to
#' generate overall raw ratio and BAF histograms.
#'
#' The plots are saved as PNG files with the following suffixes:
#' - `_res:chromosome.png`
#' - `_res:chromosome_arm.png`
#' - `_res:sample.png`
#'
#' If `file_name` is NULL, the plots are not saved to files but are returned in
#' the output list.
#'
#' @import ggplot2
#' @export
all_resolutions_plots <- function(
    data_standardized,
    sample,
    ploidy,
    centromeres,
    types_chromosome = c("Ratio_hist", "BAF_hist", "zscore"),
    types_chromosome_arm = c("Ratio_hist", "BAF_hist", "zscore"),
    types_sample = c("Ratio_hist_overall", "BAF_hist_overall"),
    file_name = NULL,
    chr = NULL) {
  # Input checks
  if (!inherits(data_standardized, "qploidy_standardization")) {
    stop("The 'data_standardized' parameter must be of class 'qploidy_standardization'.")
  }

  if (!is.character(sample) || length(sample) != 1) {
    stop("The 'sample' parameter must be a single character string.")
  }

  if (!is.numeric(ploidy)) {
    stop("The 'ploidy' parameter must be a numeric value.")
  }

  if (!is.null(chr) && !is.vector(chr)) {
    stop("The 'chr' parameter must be a vector or NULL.")
  }

  p_list <- list()
  # Raw ratio and BAF histograms (chromosome level resolution)
  p <- plot_qploidy_standardization(
    x = data_standardized,
    sample = sample,
    type = types_chromosome,
    chr = chr,
    add_expected_peaks = TRUE,
    ploidy = ploidy
  )
  p_list[[1]] <- p
  if (!is.null(file_name)) ggsave(p, filename = paste0(file_name, "_res:chromosome.png"))

  # Raw ratio and BAF histograms combining all markers in the sample (chromosome-arm level resolution)
  if (!is.null(centromeres)) {
    if (!is.vector(centromeres) || is.null(names(centromeres))) {
      stop("The 'centromeres' parameter must be a named vector with chromosome IDs as names.")
    }

    ploidy_cent <- ploidy
    if(length(ploidy) > 1 && length(ploidy) < length(chr)) ploidy_cent <- rep(ploidy, 2)

    p <- plot_qploidy_standardization(
      x = data_standardized,
      sample = sample,
      type = types_chromosome_arm,
      chr = chr,
      ploidy = ploidy_cent,
      add_expected_peaks = TRUE,
      add_centromeres = TRUE,
      centromeres = centromeres
    )

    p_list[[2]] <- p
    if (!is.null(file_name)) ggsave(p, filename = paste0(file_name, "_res:chromosome_arm.png"))
  }

  # Raw ratio and BAF histograms combining all markers in the sample (sample level resolution)
  p <- plot_qploidy_standardization(
    x = data_standardized,
    sample = sample,
    type = types_sample,
    chr = chr,
    ploidy = ploidy,
    add_expected_peaks = TRUE
  )

  p_list[[3]] <- p
  if (!is.null(file_name)) ggsave(p, filename = paste0(file_name, "_res:sample.png"))

  names(p_list) <- c("chromosome", "chromosome_arm", "sample")
  return(p_list)
}
