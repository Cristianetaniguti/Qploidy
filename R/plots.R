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
#' @importFrom magrittr "%>%"
#'
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
#' @importFrom magrittr "%>%"
#'
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
#' @importFrom magrittr "%>%"
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
    chr <- unique(data_sample$Chr)[chr] # BUGfix - in case the values are not sorted at the data
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

#' Plot X–Y scatter with ploidy dosage guide lines and optional sample highlighting
#'
#' @description
#' Creates a dot plot of allele counts `X` (A-allele) vs `Y` (B-allele) and overlays
#' expected dosage guide lines for a given ploidy. You can either color all samples
#' (`sample = "all"`) or highlight a single sample while rendering the others in gray.
#' Samples that have no plotted points (e.g., due to non-finite `X`/`Y`) are
#' automatically omitted from the legend.
#'
#' @details
#' The function:
#' * Drops rows where `X` or `Y` are non-finite before plotting and builds the legend
#'   from the remaining points.
#' * When `sample = "all"`, colors points by `SampleName` (requires a `SampleName` column).
#' * When `sample` is a specific name, only that sample is colored with
#'   `highlight_color`; all others use `other_color`. If the requested sample has
#'   no plotted points, it is omitted from the legend.
#' * Draws dosage guide lines for dosages `d = 0, …, ploidy`:
#'   `d = 0` → horizontal line `Y = 0`; `d = ploidy` → vertical line `X = 0`;
#'   intermediate dosages are lines through the origin with slope
#'   \eqn{(d/ploidy) / (1 - d/ploidy)}.
#' * Uses a fixed aspect ratio (`coord_fixed`) so that one unit on the x- and y-axes
#'   has the same length.
#'
#' @param df A `data.frame` containing at least columns `X` and `Y`. If
#'   `sample = "all"` or a specific sample is to be highlighted, `df` should also
#'   contain a `SampleName` column.
#' @param ploidy Integer (≥ 2). Ploidy used to compute and draw dosage guide lines.
#' @param sample Character. Either `"all"` to color all samples by `SampleName`,
#'   or the name of a single sample to highlight. Default: `"all"`.
#' @param highlight_color Color used for the highlighted sample when
#'   `sample != "all"`. Default: `"tomato"`.
#' @param other_color Color used for non-highlighted samples when
#'   `sample != "all"`. Default: `"grey75"`.
#'
#' @return A **ggplot** object.
#'
#'
#' @seealso [plot_baf_with_ploidy_guides()] for plotting standardized BAF-based coordinates.
#'
#' @import ggplot2
#' @export
plot_xy_with_ploidy_guides <- function(df, ploidy = 2,
                                       sample = NULL,
                                       highlight_color = "tomato",
                                       other_color = "grey75") {
  stopifnot(all(c("X","Y") %in% names(df)))
  if (!"SampleName" %in% names(df) && identical(sample, "all")) {
    stop("Column 'SampleName' is required when sample='all'.")
  }

  # Keep only rows that will actually plot
  df$.__ok__ <- is.finite(df$X) & is.finite(df$Y)
  df_plot <- df[df$.__ok__, , drop = FALSE]
  if (nrow(df_plot) == 0) stop("No points to plot after removing NAs in coordinates.")

  # Ensure SampleName is a factor with only observed levels
  if ("SampleName" %in% names(df_plot)) {
    if (is.factor(df_plot$SampleName)) {
      df_plot$SampleName <- droplevels(df_plot$SampleName)
    } else {
      df_plot$SampleName <- factor(df_plot$SampleName)
    }
  }

  max_x <- max(df_plot$X, na.rm = TRUE)
  max_y <- max(df_plot$Y, na.rm = TRUE)

  # Build guide definitions
  dosage <- 0:ploidy
  ratio  <- dosage / ploidy
  slope  <- ratio / (1 - ratio)                 # Inf when ratio == 1
  type   <- ifelse(dosage == 0, "h",
                   ifelse(dosage == ploidy, "v", "abline"))
  guides <- data.frame(dosage, ratio, slope, type)

  # Base plot (legend reflects only samples present in df_plot)
  if(is.null(sample)){
    p <- ggplot(df_plot, aes(X, Y)) +
      geom_point(alpha = 0.85, size = 2, na.rm = TRUE) +
      guides(color = guide_legend(override.aes = list(alpha = 1)))
  } else if (identical(sample, "all")) {
    p <- ggplot(df_plot, aes(X, Y)) +
      geom_point(aes(color = SampleName), alpha = 0.85, size = 2, na.rm = TRUE) +
      scale_color_discrete(drop = TRUE) +
      guides(color = guide_legend(override.aes = list(alpha = 1)))
  } else {
    if (!"SampleName" %in% names(df_plot)) stop("Column 'SampleName' not found.")
    has_highlight <- any(df_plot$SampleName == sample, na.rm = TRUE)
    if (has_highlight) {
      df_plot$.__hl__ <- factor(ifelse(df_plot$SampleName == sample, "highlight", "other"),
                                levels = c("highlight","other"))
      brks <- "highlight"; labs <- sample
    } else {
      df_plot$.__hl__ <- factor("other", levels = c("highlight","other"))
      brks <- NULL; labs <- NULL
      warning(sprintf("Sample '%s' has no plotted points. Legend will omit it.", sample))
    }

    p <- ggplot(df_plot, aes(X, Y)) +
      geom_point(aes(color = .__hl__), alpha = 0.85, size = 2, na.rm = TRUE) +
      scale_color_manual(values = c(highlight = highlight_color, other = other_color),
                         breaks = brks, labels = labs, drop = TRUE) +
      guides(color = guide_legend(title = NULL, override.aes = list(alpha = 1)))
  }

  lim <- max(max_x, max_y)            # same range on both axes
  pad <- 0.05 * lim                   # small headroom for corner labels

  p <- p +
    scale_x_continuous(limits = c(0, lim + pad), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, lim + pad), expand = c(0, 0)) +
    coord_equal(clip = "off") +       # same as coord_fixed(ratio = 1)
    labs(x = "X (reads for allele A)",
         y = "Y (reads for allele B)") +
    theme_minimal(base_size = 12) +
    theme(aspect.ratio = 1,           # make the panel square
          plot.margin = margin(12, 24, 12, 24))

  # Intermediate dosages (lines through origin)
  ab <- subset(guides, type == "abline")
  if (nrow(ab)) {
    p <- p + geom_abline(data = ab, aes(slope = slope, intercept = 0),
                         linetype = "dashed", alpha = 0.6)

    lab_x <- pmin(max_x * 0.85, ifelse(ab$slope > 0, max_y * 0.85 / ab$slope, max_x * 0.85))
    lab_y <- ab$slope * lab_x
    p <- p + annotate(
      "text",
      x = lab_x, y = lab_y,
      label = paste0("d=", ab$dosage, " (", round(ab$ratio, 2), ")"),
      hjust = 1, vjust = -0.2, size = 3
    )
  }

  # d = 0 (Y=0) and d = ploidy (X=0)
  if (any(guides$type == "h")) {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.6) +
      annotate("text", x = max_x * 0.95, y = 0, label = "d=0 (0)", vjust = -0.6, size = 3)
  }
  if (any(guides$type == "v")) {
    p <- p + geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.6) +
      annotate("text", x = 0, y = max_y * 0.95, label = paste0("d=", ploidy, " (1)"),
               hjust = -0.1, vjust = 1.1, angle = 90, size = 3)
  }

  p
}

#' Plot BAF-derived X–Y scatter with ploidy dosage guides (optional sample highlighting)
#'
#' @description
#' Reconstructs Cartesian coordinates from standardized B-allele frequency (BAF) and
#' read depth `R` as \eqn{X = (1 - \mathrm{BAF}) \cdot R} and \eqn{Y = \mathrm{BAF} \cdot R},
#' then draws a scatter plot with expected dosage guide lines for a given ploidy.
#' You can color all samples by `SampleName` or highlight a single sample and render
#' the rest in gray. Samples with no plottable points (non-finite coordinates) are
#' automatically omitted from the legend.
#'
#' @details
#' * Coordinates are computed from BAF and depth: \eqn{X=(1-\mathrm{BAF})R}, \eqn{Y=\mathrm{BAF}R}.
#' * If `fallback_to_ratio = TRUE` and `baf` is `NA`, values from `ratio` are used.
#'   The effective BAF is clamped into \[0, 1\].
#' * When `normalize_depth = TRUE`, all points are projected to the same radius
#'   (depth) given by `radius` (or `stats::median(R)` if `radius` is `NULL`), which
#'   emphasizes dosage bands rather than depth variation. When `FALSE`, each point
#'   uses its own `R`.
#' * Dosage guide lines are drawn for \eqn{d \in \{0,\dots,\mathrm{ploidy}\}}:
#'   `d = 0` → horizontal line `Y = 0`; `d = ploidy` → vertical line `X = 0`;
#'   intermediate dosages are lines through the origin with slope
#'   \eqn{(d/\mathrm{ploidy}) / (1 - d/\mathrm{ploidy})}.
#' * The legend is built from actually plotted rows only; if a requested `sample`
#'   has no plottable points, it is omitted from the legend.
#' * Uses a fixed aspect ratio (`coord_fixed`) so x and y units are comparable.
#'
#' @param df A `data.frame` with required columns:
#'   - `baf` (numeric in \[0,1\]): standardized B-allele frequency.
#'   - `R` (numeric): total read depth.
#'   - `SampleName` (character/factor): sample label used for coloring.
#'   Optional column `ratio` may be present and used when `fallback_to_ratio = TRUE`.
#' @param ploidy Integer (≥ 2). Ploidy used to compute dosage guide lines.
#' @param fallback_to_ratio Logical. If `TRUE`, fill `NA` values in `baf` with
#'   corresponding values from `ratio` (when available). Default: `FALSE`.
#' @param normalize_depth Logical. If `TRUE`, place all points on a common radius
#'   (see `radius`); if `FALSE`, use each point's `R`. Default: `TRUE`.
#' @param radius Numeric scalar radius to use when `normalize_depth = TRUE`.
#'   If `NULL`, uses `stats::median(df$R, na.rm = TRUE)`. Ignored when
#'   `normalize_depth = FALSE`. Default: `NULL`.
#' @param sample Character. Either `"all"` to color all samples by `SampleName`,
#'   or the name of a single sample to highlight. Default: `"all"`.
#' @param highlight_color Color for the highlighted sample when `sample != "all"`.
#'   Default: `"tomato"`.
#' @param other_color Color for non-highlighted samples when `sample != "all"`.
#'   Default: `"grey75"`.
#'
#' @return A **ggplot** object.
#'
#'
#' @seealso [plot_xy_with_ploidy_guides()] for plotting raw `X`/`Y` counts with guides.
#'
#' @import ggplot2
#' @importFrom stats median
#' @export
plot_baf_with_ploidy_guides <- function(df,
                                        ploidy = 2,
                                        fallback_to_ratio = FALSE,
                                        normalize_depth = TRUE,
                                        radius = NULL,
                                        sample = NULL,
                                        highlight_color = "tomato",
                                        other_color = "grey75") {
  stopifnot("R" %in% names(df))
  if (!"baf" %in% names(df)) stop("Column 'baf' not found.")
  if (!"SampleName" %in% names(df)) stop("Column 'SampleName' not found (needed for coloring).")

  # Effective BAF
  baf_eff <- df$baf
  if (fallback_to_ratio && "ratio" %in% names(df)) {
    baf_eff[is.na(baf_eff)] <- df$ratio[is.na(baf_eff)]
  }
  baf_eff <- pmin(pmax(baf_eff, 0), 1)

  # Depth to use
  if (normalize_depth) {
    if (is.null(radius)) radius <- stats::median(df$R, na.rm = TRUE)
    Ruse <- rep(radius, nrow(df))
  } else {
    Ruse <- df$R
  }

  # Coordinates from BAF and depth
  df$Xb <- (1 - baf_eff) * Ruse
  df$Yb <- baf_eff * Ruse

  # Keep only points that will actually plot
  df$.__ok__ <- is.finite(df$Xb) & is.finite(df$Yb)
  df_plot <- df[df$.__ok__, , drop = FALSE]
  if (nrow(df_plot) == 0) stop("No points to plot after removing NAs in coordinates.")

  # Guides
  dosage <- 0:ploidy
  ratio  <- dosage / ploidy
  slope  <- ratio / (1 - ratio)
  type   <- ifelse(dosage == 0, "h", ifelse(dosage == ploidy, "v", "abline"))
  guides <- data.frame(dosage, ratio, slope, type)

  max_x <- max(df_plot$Xb, na.rm = TRUE)
  max_y <- max(df_plot$Yb, na.rm = TRUE)

  # Coloring
  if(is.null(sample)){
    p <- ggplot(df_plot, aes(Xb, Yb)) +
      geom_point(alpha = 0.85, size = 2, na.rm = TRUE) +
      guides(color = guide_legend(override.aes = list(alpha = 1)))
  } else   if (identical(sample, "all")) {
    p <- ggplot(df_plot, aes(Xb, Yb)) +
      geom_point(aes(color = SampleName), alpha = 0.85, size = 2, na.rm = TRUE) +
      scale_color_discrete(drop = TRUE) +
      guides(color = guide_legend(override.aes = list(alpha = 1)))
  } else {
    has_highlight <- any(df_plot$SampleName == sample, na.rm = TRUE)
    if (has_highlight) {
      df_plot$.__hl__ <- factor(ifelse(df_plot$SampleName == sample, "highlight", "other"),
                                levels = c("highlight","other"))
      brks <- "highlight"; labs <- sample
    } else {
      df_plot$.__hl__ <- factor("other", levels = c("highlight","other"))
      brks <- NULL; labs <- NULL
    }
    p <- ggplot(df_plot, aes(Xb, Yb)) +
      geom_point(aes(color = .__hl__), alpha = 0.85, size = 2, na.rm = TRUE) +
      scale_color_manual(values = c(highlight = highlight_color, other = other_color),
                         breaks = brks, labels = labs, drop = TRUE) +
      guides(color = guide_legend(title = NULL, override.aes = list(alpha = 1)))
  }


  lim <- max(max_x, max_y)            # same range on both axes
  pad <- 0.05 * lim                   # small headroom for corner labels
  # Rest of the plot
  p <- p +
    scale_x_continuous(limits = c(0, lim + pad), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, lim + pad), expand = c(0, 0)) +
    coord_equal(clip = "off") +       # same as coord_fixed(ratio = 1)
    labs(
      x = "X = (1 - BAF) * R",
      y = "Y = BAF * R"
    ) +
    theme_minimal(base_size = 12) +
    theme(aspect.ratio = 1,           # make the panel square
          plot.margin = margin(12, 24, 12, 24))

  # Intermediate dosage guides
  ab <- subset(guides, type == "abline" & is.finite(slope))
  if (nrow(ab)) {
    p <- p + geom_abline(data = ab, aes(slope = slope, intercept = 0),
                         linetype = "dashed", alpha = 0.6)
    lab_x <- pmin(max_x * 0.85, ifelse(ab$slope > 0, max_y * 0.85 / ab$slope, max_x * 0.85))
    lab_y <- ab$slope * lab_x
    p <- p + annotate("text", x = lab_x, y = lab_y,
                      label = paste0("d=", ab$dosage, " (", round(ab$ratio, 2), ")"),
                      hjust = 1, vjust = -0.2, size = 3)
  }

  # d = 0 and d = ploidy
  if (any(guides$type == "h")) {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.6) +
      annotate("text", x = max_x * 0.95, y = 0, label = "d=0 (0)", vjust = -0.6, size = 3)
  }
  if (any(guides$type == "v")) {
    p <- p + geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.6) +
      annotate("text", x = 0, y = max_y * 0.95, label = paste0("d=", ploidy, " (1)"),
               hjust = -0.1, vjust = 1.1, angle = 90, size = 3)
  }

  p
}


