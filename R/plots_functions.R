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
                     colors = FALSE){

  if(add_centromeres){
    #centromeres <- c("1 = 49130338, 5 = 49834357")
    centromeres <- gsub("\ ", "", centromeres)
    centromeres <- unlist(strsplit(centromeres, ","))
    centromeres <- sapply(centromeres, function(x) strsplit(x, "="))
    centromeres_df <- data.frame(Chr = sapply(centromeres, "[[", 1), value = sapply(centromeres, "[[", 2))
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
          legend.position="none") +
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
#'
#' @export
plot_baf_hist <- function(data_sample,
                          area_single,
                          ploidy,
                          colors = FALSE,
                          add_estimated_peaks = TRUE,
                          add_expected_peaks = FALSE,
                          BAF_hist_overall = TRUE){


  if(add_estimated_peaks | add_expected_peaks | colors){
    if(length(ploidy) == 1) {
      ploidy <- rep(ploidy, length(unique(data_sample$Chr)))
    } else if(length(ploidy) != 1 & length(ploidy) != length(unique(data_sample$Chr)))
      stop("Provide a ploidy for each chromosome.")

    data_sample2 <- modes.df3 <- data.frame()
    for(z in 1:length(unique(data_sample$Chr))){
      ymin <- seq(0, 1, 1/ploidy[z]) - (area_single/(ploidy[z]*2))
      ymax <- seq(0, 1, 1/ploidy[z]) + (area_single/(ploidy[z]*2))

      ymin[which(ymin < 0)] <- 0
      ymax[which(ymax > 1)] <- 1
      rets <- data.frame(ymin, ymax, xmax = Inf, xmin = -Inf, Chr = sort(unique(data_sample$Chr))[z])

      split_data <- data_sample[which(data_sample$Chr == sort(unique(data_sample$Chr))[z]),]
      idx_tot <- FALSE
      modes.df2 <- data.frame()
      for(i in 1:nrow(rets)) {
        idx <- split_data$sample >= rets$ymin[i] & split_data$sample <= rets$ymax[i]
        idx_tot <- idx_tot | idx

        estimated <- mode(split_data$sample[which(split_data$sample > rets$ymin[i] & split_data$sample < rets$ymax[i])])
        expected <- seq(0, 1, 1/ploidy[z])[i]
        modes.df <- cbind(Chr = unique(split_data$Chr), pivot_longer(data.frame(estimated, expected), cols = 1:2))
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
  }

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
          legend.position="bottom") +
    {if(add_expected_peaks | add_estimated_peaks) labs(color= "Peaks")} +
    {if(colors) labs(fill= "Area")}

  if(BAF_hist_overall){
    p_hist_all <- data_sample %>% ggplot(aes(x=sample)) + geom_histogram() +
      scale_color_manual(values = c("blue", "purple")) +
      scale_linetype_manual(values = c("dashed", "solid"), guide="none") +
      scale_alpha(range = c(0.7, 1), guide="none") +
      theme_bw() +  xlab("BAF") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.position="bottom")

    p_hist <- ggarrange(p_hist_all, p_hist, widths = c(1,length(unique(data_sample$Chr))))
  }

  return(p_hist)
}
