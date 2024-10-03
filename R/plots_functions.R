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
#' @param ratio if TRUE plot the raw ratio
#' @param rm_homozygous if TRUE removes the homozygous calls
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
                          rm_homozygous = FALSE){

  if(add_estimated_peaks | add_expected_peaks | colors){
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
            legend.position="bottom") +
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
            legend.position="bottom") +
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
