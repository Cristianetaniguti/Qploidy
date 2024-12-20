#' Plot BAF
#'
#' @param data_sample data.frame with BAF and genomic position information
#' @param area_single area around the expected peak to be considered
#' @param ploidy expected ploidy overall (single integer value) or for each chromosome
#' @param dot.size graphic dot size
#' @param add_estimated_peaks add expected peaks lines
#' @param add_expected_peaks add estimated peaks
#' @param colors add area colors
#' @param centromeres vector defining centromeres positions
#' @param add_centromeres logical defining if centromeres positions will be displayed
#' @param font_size graphic labels font size
#'
#' @import ggplot2
#' @import tidyr
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
                     font_size = 12){

  if(add_centromeres){
    #centromeres <- c("1 = 49130338, 5 = 49834357")
    if(length(centromeres) == 1 & any(grepl(",", centromeres))){
      centromeres <- gsub("\ ", "", centromeres)
      centromeres <- unlist(strsplit(centromeres, ","))
      centromeres <- sapply(centromeres, function(x) strsplit(x, "="))
      centromeres_df <- data.frame(Chr = sapply(centromeres, "[[", 1), value = sapply(centromeres, "[[", 2))
    } else {
      centromeres_df <- data.frame(Chr = names(centromeres), value = centromeres)
    }
  }

  if(add_estimated_peaks | add_expected_peaks | colors){

    if(length(ploidy) == 1) {
      ploidy <- rep(ploidy, length(unique(data_sample$Chr)))
    } else if(length(ploidy) > 1 & length(ploidy) < length(unique(data_sample$Chr)))
      stop("Provide a ploidy for each chromosome.")

    data_sample2 <- rets2 <- data.frame()
    for(z in 1:length(unique(data_sample$Chr))){
      ymin <- seq(0, 1, 1/ploidy[z]) - (area_single/(ploidy[z]*2))
      ymax <- seq(0, 1, 1/ploidy[z]) + (area_single/(ploidy[z]*2))

      ymin[which(ymin < 0)] <- 0
      ymax[which(ymax > 1)] <- 1
      rets <- data.frame(ymin, ymax, xmax = Inf, xmin = -Inf, Chr = sort(unique(data_sample$Chr))[z])

      data_chr <- data_sample[which(data_sample$Chr == sort(unique(data_sample$Chr))[z]),]
      idx_tot <- FALSE
      idx <- list()
      for(i in 1:nrow(rets)){
        idx <- data_chr$sample >= rets$ymin[i] & data_chr$sample <= rets$ymax[i]
        idx_tot <- idx_tot | idx
      }

      data_chr$color <- NA
      data_chr$color[which(!idx_tot)] <- "red"
      data_chr$color[which(idx_tot)] <- "black"

      rets2 <- rbind(rets2, rets)
      data_sample2 <- rbind(data_sample2, data_chr)
    }
  } else data_sample2 <- data_sample

  p_baf <- data_sample2 %>% ggplot(aes(x=Position, y=sample)) +
    {if(colors) geom_point(aes(color = color), alpha=0.7, size=dot.size) else geom_point(alpha=0.7, size=dot.size)} +
    scale_color_manual(values = c("red", "black")) +
    facet_grid(~ Chr, scales = "free_x") + theme_bw() +
    {if(add_expected_peaks) geom_rect(data = rets2, inherit.aes=FALSE,
                                      aes(ymin=ymin,ymax=ymax,
                                          xmax = xmax, xmin = xmin,
                                          alpha=0.001),
                                      color = "red")} +
    ylab("BAF") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position="none", text = element_text(size = font_size)) +
    {if(add_centromeres) geom_vline(data = centromeres_df,
                                    aes(xintercept= as.numeric(value)),
                                    color = "blue",
                                    linewidth = 0.8)}
  return(p_baf)
}

#' Plot BAF histogram
#'
#' @param data_sample data.frame with BAF and genomic position information
#' @param area_single area around the expected peak to be considered
#' @param ploidy expected ploidy
#' @param add_estimated_peaks add expected peaks lines
#' @param add_expected_peaks add estimated peaks
#' @param colors add area colors
#' @param BAF_hist_overall if TRUE it plots the BAF histogram for entire genome
#' @param ratio if TRUE plot the raw ratio
#' @param rm_homozygous if TRUE removes the homozygous calls
#' @param font_size graphic labels font size
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
                          font_size = 12){

  if((add_estimated_peaks | add_expected_peaks | colors) & !all(is.na(ploidy))){
    if(length(ploidy) == 1) {
      ploidy <- rep(ploidy, length(unique(data_sample$Chr)))
    } else if(length(ploidy) != 1 & length(ploidy) != length(unique(data_sample$Chr)))
      stop("Provide a ploidy for each chromosome.")

    if(BAF_hist_overall) Chrs <- "all" else Chrs = sort(unique(data_sample$Chr))
    data_sample2 <- modes.df3 <- data.frame()
    for(z in 1:length(Chrs)){
      ymin <- seq(0, 1, 1/ploidy[z]) - (area_single/(ploidy[z]*2))
      ymax <- seq(0, 1, 1/ploidy[z]) + (area_single/(ploidy[z]*2))

      ymin[which(ymin < 0)] <- 0
      ymax[which(ymax > 1)] <- 1
      rets <- data.frame(ymin, ymax, xmax = Inf, xmin = -Inf, Chr = Chrs[z])

      if(all(Chrs == "all")) {
        split_data <- data_sample
      } else {
        split_data <- data_sample[which(data_sample$Chr == sort(unique(data_sample$Chr))[z]),]
      }

      if(ratio) split_data$sample <- split_data$ratio

      idx_tot <- FALSE
      modes.df2 <- data.frame()
      for(i in 1:nrow(rets)) {
        idx <- split_data$sample >= rets$ymin[i] & split_data$sample <= rets$ymax[i]
        idx_tot <- idx_tot | idx

        estimated <- mode(split_data$sample[which(split_data$sample > rets$ymin[i] & split_data$sample < rets$ymax[i])])
        expected <- seq(0, 1, 1/ploidy[z])[i]
        modes.df <- cbind(Chr = Chrs[z], pivot_longer(data.frame(estimated, expected), cols = 1:2))
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

  if(rm_homozygous) data_sample2 <- data_sample2 %>% filter(sample != 0 & sample != 1)
  if(ratio) data_sample2$sample <- data_sample2$ratio

  if(!BAF_hist_overall){
    p_hist <- data_sample2 %>% ggplot(aes(x=sample)) +
      {if(colors) geom_histogram(aes(fill = color)) else  geom_histogram()} +
      #scale_x_continuous(breaks = round(seq(0, 1, 1/ploidy),2)) +
      {if(add_estimated_peaks) geom_vline(data = modes.df3[which(modes.df3$name == "estimated"),],
                                          aes(xintercept= value,
                                              color = name,
                                              linetype = name,
                                              alpha = alpha),
                                          linewidth = 0.8)} +
      {if(add_expected_peaks) geom_vline(data = modes.df3[which(modes.df3$name == "expected"),],
                                         aes(xintercept= value,
                                             color = name,
                                             linetype = name,
                                             alpha = alpha),
                                         linewidth = 0.8)} +
      {if(colors) scale_fill_manual(values = c("red", "black"))} +
      {if(add_expected_peaks | add_estimated_peaks) scale_color_manual(values = c("blue", "purple"))}+
      scale_linetype_manual(values = c("dashed", "solid"), guide="none") +
      scale_alpha(range = c(0.7, 1), guide="none") +
      facet_grid(~ Chr, scales = "free_x") + theme_bw() +  xlab("BAF") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.position="bottom", text = element_text(size = font_size)) +
      {if(add_expected_peaks | add_estimated_peaks) labs(color= "Peaks")} +
      {if(colors) labs(fill= "Area")}

  } else {
    if(rm_homozygous) data_sample <- data_sample %>% filter(sample != 0 & sample != 1)
    if(ratio) data_sample$sample <- data_sample$ratio

    p_hist_all <- data_sample %>% ggplot(aes(x=sample)) + geom_histogram() +
      {if(colors) geom_histogram(aes(fill = color)) else  geom_histogram()} +
      #scale_x_continuous(breaks = round(seq(0, 1, 1/ploidy),2)) +
      {if(add_estimated_peaks) geom_vline(data = modes.df3[which(modes.df3$name == "estimated"),],
                                          aes(xintercept= value,
                                              color = name,
                                              linetype = name,
                                              alpha = alpha),
                                          linewidth = 0.8)} +
      {if(add_expected_peaks) geom_vline(data = modes.df3[which(modes.df3$name == "expected"),],
                                         aes(xintercept= value,
                                             color = name,
                                             linetype = name,
                                             alpha = alpha),
                                         linewidth = 0.8)} +
      {if(colors) scale_fill_manual(values = c("red", "black"))} +
      {if(add_expected_peaks | add_estimated_peaks) scale_color_manual(values = c("blue", "purple"))}+
      scale_linetype_manual(values = c("dashed", "solid"), guide="none") +
      scale_alpha(range = c(0.7, 1), guide="none") +
      theme_bw() +  {if(ratio) xlab("ratio") else xlab("BAF")} +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.position="bottom", text = element_text(size = font_size)) +
      {if(add_expected_peaks | add_estimated_peaks) labs(color= "Peaks")} +
      {if(colors) labs(fill= "Area")}

    p_hist <- p_hist_all
  }

  return(p_hist)
}


#' 3D plot to check raw data allele intensities (if array) or depth (if sequencing) data properties
#' Qploidy can be used if that is not much variation among the individuals for the same marker.
#' So we expect that the overall structure of the resulting graphic is a ramp.
#'
#' @param data_qploidy data.frame input for Qploidy
#' @param R_lim limit the sum of intensities/depth axis to adjust the dimentions in case there are outliers
#' @param n subset size
#'
#'
#' @importFrom plotly plot_ly layout
#' @import tidyr
#'
#' @export
plot_check <- function(data_qploidy, R_lim = NULL, n = 3000){

  plot_df <- pivot_wider(data_qploidy[,c(1,2,5)], names_from = SampleName, values_from = R)
  mk_names <- plot_df$MarkerName
  plot_df <- as.matrix(plot_df[,-1])
  rownames(plot_df) <- mk_names

  plot_df <- plot_df[sample(1:nrow(plot_df), n, replace =FALSE),]
  plot_df <- plot_df[order(apply(plot_df, 1, sum)),]

  if(!is.null(R_lim)) {
    plot_df[which(plot_df > R_lim)] <- NA
    rm.mk <- apply(plot_df, 1, function(x) sum(is.na(x))/length(x))
    rm.mk <- which(rm.mk > 0.5)
    if(length(rm.mk)>0) plot_df <- plot_df[-rm.mk,]
  }

  axx <- list(title = "Individuals")
  axy <- list(title = "Markers")
  axz <- list(title = "Sum of intensities/depth")

  cat("Values summary:\n")
  print(summary(as.numeric(plot_df)))

  p <- plot_ly(z=plot_df, type="surface") %>%
    layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  return(p)
}

#' Plot method for object of class 'qploidy_standardization'
#'
#' @param x object of class 'qploidy_standardization'
#' @param sample character indicating sample ID
#' @param area_single area around the expected peak to be considered
#' @param ploidy expected ploidy
#' @param add_estimated_peaks add expected peaks lines
#' @param add_expected_peaks add estimated peaks
#' @param colors add area colors
#' @param type character defining the graphic type. "all" plots all graphics,
#'             it is equivalent to "het","ratio","BAF","zscore","BAF_hist", "BAF_hist_overall".
#'             "het" plots the heterozygous locus counts inside genomic positions window defined in 'window_size'.
#'             "ratio" plots raw Y/(X+Y) or alternative count/(alternative counts + reference counts).
#'             "BAF" plots standardized ratios. "zscore" is the smoothed conditional means curve of standardized sum of intensities/counts.
#'             "BAF_hist" is the histogram of standardized ratios. "BAF_hist_overall" add the histogram including all markers.
#' @param window_size genomic position window to calculate the number of heterozygous locus
#' @param het_interval interval to be considered as heterozygous (heterozygous ratio > 1 - het_interval and heterozygous ratio < 0 + het_interval)
#' @param chr chromosome index to be plotted
#' @param dot.size dot size
#' @param centromeres centromeres position for each chromosome
#' @param add_centromeres if TRUE add a vertical line at the centromere position
#' @param rm_homozygous if TRUE remove the homozygous for the BAF histogram plots
#' @param font_size graphic labels font size
#'
#' @param ... plot parameters
#'
#' @importFrom ggpubr ggarrange annotate_figure text_grob
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#'
#' @return printed information about Qploidy standardization process
#'
#'
#' @export
plot_qploidy_standardization <- function(x,
                                         sample = NULL,
                                         chr = NULL,
                                         type = c("all", "het", "BAF","zscore","BAF_hist", "BAF_hist_overall", "Ratio_hist_overall"),
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
                                         rm_homozygous = FALSE, ...){

  if(!inherits(x, "qploidy_standardization")) stop("Object is not of class qploidy_standardization")

  if(is.null(sample)) stop("Define sample ID")

  data_sample <- x$data[which(x$data$SampleName == sample),]

  if(is.numeric(chr)) chr <- sort(unique(data_sample$Chr))[chr] else if(is.null(chr)) chr <- sort(unique(data_sample$Chr))

  data_sample <- data_sample %>% filter(Chr %in% chr) %>% select(MarkerName, SampleName, Chr, Position, baf, z, ratio)

  if(nrow(data_sample) == 0) stop("Sample or chromosome not found.")

  baf_sample <- data_sample %>% pivot_wider(names_from = SampleName, values_from = baf)
  zscore_sample <- data_sample %>% pivot_wider(names_from = SampleName, values_from = z)

  colnames(baf_sample)[ncol(baf_sample)] <-  "sample"

  baf_point <- baf_hist <- p_z <- raw_ratio <- het_rate <- baf_hist_overall <- ratio_hist_overall <- NULL

  if(any(type == "all" | type == "BAF")){
    baf_point <- plot_baf(baf_sample,
                          area_single,
                          ploidy,
                          dot.size = dot.size,
                          add_estimated_peaks,
                          add_expected_peaks,
                          centromeres,
                          add_centromeres,
                          colors,
                          font_size = font_size)
  }

  if(add_centromeres){
    #centromeres <- c("1 = 49130338, 5 = 49834357")
    if(length(centromeres) == 1 & any(grepl(",", centromeres))){
      centromeres <- gsub("\ ", "", centromeres)
      centromeres <- unlist(strsplit(centromeres, ","))
      centromeres <- sapply(centromeres, function(x) strsplit(x, "="))
      centromeres_df <- data.frame(Chr = sapply(centromeres, "[[", 1), value = sapply(centromeres, "[[", 2))
    } else {
      centromeres_df <- data.frame(Chr = names(centromeres), value = centromeres)
    }

    for(i in 1:length(centromeres)){
      idx <- which(baf_sample$Chr == centromeres_df$Chr[i] & baf_sample$Position < centromeres_df$value[i])
      baf_sample$Chr[idx] <- paste0(baf_sample$Chr[idx], ".1")
      idx <- which(baf_sample$Chr == centromeres_df$Chr[i] & baf_sample$Position >= centromeres_df$value[i])
      baf_sample$Chr[idx] <- paste0(baf_sample$Chr[idx], ".2")
      if(any(type == "zscore")){
        idx <- which(zscore_sample$Chr == centromeres_df$Chr[i] & zscore_sample$Position < centromeres_df$value[i])
        zscore_sample$Chr[idx] <- paste0(zscore_sample$Chr[idx], ".1")
        idx <- which(zscore_sample$Chr == centromeres_df$Chr[i] & zscore_sample$Position >= centromeres_df$value[i])
        zscore_sample$Chr[idx] <- paste0(zscore_sample$Chr[idx], ".2")
      }
    }
  }

  if(any(type == "all" | type == "BAF_hist")){
    baf_hist <- plot_baf_hist(data_sample = baf_sample,
                              area_single,
                              ploidy,
                              colors,
                              add_estimated_peaks,
                              add_expected_peaks,
                              BAF_hist_overall = FALSE,
                              rm_homozygous = rm_homozygous,
                              font_size = font_size)
  }

  if(any(type == "all" | type == "BAF_hist_overall")){
    baf_hist_overall <- plot_baf_hist(data_sample = baf_sample,
                                      area_single,
                                      ploidy,
                                      colors,
                                      add_estimated_peaks,
                                      add_expected_peaks,
                                      BAF_hist_overall = TRUE,
                                      rm_homozygous = rm_homozygous,
                                      font_size = font_size)
  }

  if(any(type == "all" | type == "Ratio_hist_overall")){
    ratio_hist_overall <- plot_baf_hist(data_sample = baf_sample,
                                        area_single,
                                        ploidy,
                                        colors,
                                        add_estimated_peaks,
                                        add_expected_peaks,
                                        BAF_hist_overall = TRUE,
                                        ratio = TRUE,
                                        rm_homozygous = rm_homozygous,
                                        font_size = font_size)
  }

  if(any(type == "zscore")){
    colnames(zscore_sample)[ncol(zscore_sample)] <- "z"
    p_z <- zscore_sample  %>%
      ggplot(aes(x = Position , y = z)) +
      facet_grid(.~Chr, scales = "free") +
      geom_smooth(method = "gam") +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                         text = element_text(size = font_size)) +
      geom_hline(yintercept=median(zscore_sample$z), linetype="dashed")

  }

  if(any(type == "ratio")){
    raw_ratio <- data_sample %>%  ggplot(aes(x = Position, y = ratio)) + geom_point(alpha =0.7, size=dot.size) +
      facet_grid(.~Chr, scales = "free") + theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            text = element_text(size = font_size))
  }

  if(any(type == "het")){
    data_sample$het_rate <- NA
    data_sample$het_rate[which(data_sample$ratio <= 0+het_interval | data_sample$ratio >= 1 - het_interval)] <- 0
    data_sample$het_rate[which(data_sample$ratio > 0+het_interval & data_sample$ratio < 1 - het_interval)] <- 1

    data_sample_lst <- split.data.frame(data_sample, data_sample$Chr)

    for(j in 1:length(data_sample_lst)){
      data_sample_lst[[j]]$window <- 1
      if(length(unique(data_sample_lst[[j]]$Position)) == 1) next()
      intervals_start <- c(seq(1,max(data_sample_lst[[j]]$Position), window_size))
      intervals_end <- c(seq(1,max(data_sample_lst[[j]]$Position), window_size) - 1, max(data_sample_lst[[j]]$Position))
      intervals_end <- intervals_end[-1]
      for(i in 1:length(intervals_start)){
        data_sample_lst[[j]]$window[which(data_sample_lst[[j]]$Position >= intervals_start[i] & data_sample_lst[[j]]$Position <= intervals_end[i])] <- intervals_start[i]
      }
    }

    data_sample <- do.call(rbind, data_sample_lst)

    het_rate <- data_sample %>%
      group_by(Chr, window) %>%
      summarise(prop_het = sum(het_rate, na.rm = TRUE)/length(which(het_rate == 1 | het_rate == 0))) %>%
      ggplot(aes(x = window, y = prop_het, color = prop_het)) + geom_line() +
      scale_color_gradient(low = "black", high = "red") +
      facet_grid(.~Chr, scales = "free") + theme_bw() + ylab("proportion of heterozygous loci") + xlab("Position")+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.position = "none", text = element_text(size = font_size))
  }

  p_all <- list(het_rate, raw_ratio, baf_point, baf_hist,baf_hist_overall, ratio_hist_overall, p_z)

  rm <- which(sapply(p_all, is.null))
  if(length(rm) != 0) p_all <- p_all[-rm]

  p_result <- ggarrange(plotlist = p_all, ncol = 1)

  p_result <- annotate_figure(p_result, top = text_grob(sample, face = "bold", size = 14))

  return(p_result)
}

#' Plot graphics for ploidy visual inspection for each resolution
#' It was made for parallelization purpose
#'
#' @param data_standardized result from standardization function
#' @param sample sample name
#' @param ploidy sample ploidy if known
#' @param centromeres vector with centromeres positions
#' @param file_name character defining the output file path and name
#' @param chr vector containing which chromosomes should be included
#'
#' @importFrom ggplot2 ggsave
#'
#' @export
all_resolutions_plots <- function(data_standardized, sample, ploidy, centromeres, file_name, chr){
  # Raw ratio and BAF histograms (chromosome level resolution)
  p <- plot_qploidy_standardization(x = data_standardized,
                                    sample = sample,
                                    type = c("Ratio_hist", "BAF_hist", "zscore"),
                                    chr = chr,
                                    add_expected_peaks = TRUE,
                                    ploidy = ploidy)

  ggsave(p, filename = paste0(file_name, "_res:chromosome.png"))

  # Raw ratio and BAF histograms combining all markers in the sample (chromosome-arm level resolution)
  p <- plot_qploidy_standardization(x = data_standardized,
                                    sample = sample,
                                    type = c("Ratio_hist", "BAF_hist", "zscore"),
                                    chr = chr,
                                    ploidy = ploidy,
                                    add_expected_peaks = TRUE,
                                    add_centromeres = TRUE,
                                    centromeres = centromeres)

  ggsave(p, filename = paste0(file_name, "_res:chromosome_arm.png"))

  # Raw ratio and BAF histograms combining all markers in the sample (sample level resolution)
  p <- plot_qploidy_standardization(x = data_standardized,
                                    sample = sample,
                                    type = c("Ratio_hist_overall", "BAF_hist_overall"),
                                    chr = chr,
                                    ploidy = ploidy,
                                    add_expected_peaks = TRUE)

  ggsave(p, filename = paste0(file_name, "_res:sample.png"))
}
