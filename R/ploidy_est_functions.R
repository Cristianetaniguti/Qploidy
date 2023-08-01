
#' Estimate ploidy using area method
#'
#' @param data_sample data.frame with BAF values and genomic position information
#' @param ploidys range of ploidies to by tested
#' @param area area around the expected peak to be considered
#'
#' @export
area_estimate_ploidy_by_chr <- function(data_sample, ploidys, area){

  by_chr <- split(data_sample, data_sample$Chr)

  freq <- pascalTriangle(ploidys[2])
  freq <- freq[-1]
  ploidys <- ploidys[1]:ploidys[2]
  p.values_tot <- vector()
  dots.int_tot <- corr_tot <- max_sd_tot <- modes_paste_tot <- as.list(rep(NA, length(by_chr)))
  means <- list()
  for(j in 1:length(ploidys)){
    # Area method
    ymin <- seq(0, 1, 1/ploidys[j]) - (area/(ploidys[j]*2))
    ymax <- seq(0, 1, 1/ploidys[j]) + (area/(ploidys[j]*2))

    ymin[which(ymin < 0)] <- 0
    ymax[which(ymax > 1)] <- 1
    rets <- data.frame(ymin, ymax)

    prop_tot <- modes_tot  <- sd_areas_tot <- as.list(rep(NA, length(by_chr)))
    for(z in 1:length(by_chr)){
      for(i in 1:nrow(rets)){
        prop <- apply(by_chr[[z]][,-c(1,2)], 2, function(x) sum(x >= rets$ymin[i] & x <= rets$ymax[i]))
        modes <- apply(by_chr[[z]][,-c(1,2)], 2, function(x) median(x[x >= rets$ymin[i] & x <= rets$ymax[i]]))
        sd_areas <- apply(by_chr[[z]][,-c(1,2)], 2, function(x) sd(x[x >= rets$ymin[i] & x <= rets$ymax[i]]))

        prop_tot[[z]] <- rbind(prop_tot[[z]], prop)
        modes_tot[[z]] <- rbind(modes_tot[[z]], modes)
        sd_areas_tot[[z]] <- rbind(sd_areas_tot[[z]], sd_areas)
      }
    }

    # remove NAs
    prop_tot <- lapply(prop_tot, function(x) x[-1,])
    for(z in 1:length(prop_tot)){
      if(is.null(ncol(prop_tot[[z]]))) dots.int <- sum(prop_tot[[z]])/dim(by_chr[[z]])[1] else
        dots.int <- apply(prop_tot[[z]], 2, function(x) sum(x)/dim(by_chr[[z]])[1])
      dots.int_tot[[z]] <- rbind(dots.int_tot[[z]], dots.int)

      modes_paste <- apply(modes_tot[[z]], 2, function(x) paste0(round(x[-1],3), collapse = "/"))
      modes_paste_tot[[z]] <- rbind(modes_paste_tot[[z]], modes_paste)

      corr <- apply(modes_tot[[z]], 2, function(x) cor(x = x[-1], y = seq(0,1,1/ploidys[j])))
      corr_tot[[z]] <- rbind(corr_tot[[z]], corr)

      max_sd <- apply(sd_areas_tot[[z]], 2, function(x) max(x[-1]))
      max_sd_tot[[z]] <- rbind(max_sd_tot[[z]], max_sd)
    }
  }

  dots.int_tot <- lapply(dots.int_tot, function(x) x[-1,])
  modes_paste_tot <- lapply(modes_paste_tot, function(x) x[-1,])
  corr_tot <- lapply(corr_tot, function(x) x[-1,])
  max_sd_tot <- lapply(max_sd_tot, function(x) x[-1,])

  # Area method
  result.ploidy  <- diff.count <- second <- diff.second <- diff.first.second <- list()
  for(z in 1:length(dots.int_tot)){
    if(is.null(rownames(dots.int_tot[[z]]))) {
      names(dots.int_tot[[z]]) <- ploidys
      result.ploidy[[z]] <- ploidys[order(dots.int_tot[[z]], decreasing =T)][1]

      diff.count[[z]] <- dots.int_tot[[z]][order(dots.int_tot[[z]], decreasing =T)][1]
      second[[z]] <- ploidys[order(dots.int_tot[[z]], decreasing =T)][2]
      diff.second[[z]] <- dots.int_tot[[z]][order(dots.int_tot[[z]], decreasing =T)][2]
    } else {
      rownames(dots.int_tot[[z]]) <- ploidys
      result.ploidy[[z]] <- apply(dots.int_tot[[z]], 2, function(x) ploidys[order(x, decreasing =T)][1])
      diff.count[[z]] <- apply(dots.int_tot[[z]], 2, function(x) x[order(x, decreasing =T)][1])
      second[[z]] <- apply(dots.int_tot[[z]], 2, function(x) ploidys[order(x, decreasing =T)][2])
      diff.second[[z]] <- apply(dots.int_tot[[z]], 2, function(x) x[order(x, decreasing =T)][2])
    }
    diff.first.second[[z]] <- diff.count[[z]] - diff.second[[z]]
    # filt <- which((diff.count[[z]] - diff.second[[z]]) < filter)
    # result.ploidy[[z]][filt] <- NA
    # second[[z]][filt] <- NA
  }

  result.ploidy <- do.call(rbind, result.ploidy)
  diff.first.second <- do.call(rbind, diff.first.second)

  if(is.vector(max_sd_tot[[1]])){ # only one individual selected
    max_sd_tot <- lapply(max_sd_tot, as.matrix)
    corr_tot <- lapply(corr_tot, as.matrix)
    modes_paste_tot <- lapply(modes_paste_tot, as.matrix)
    dots.int_tot <- lapply(dots.int_tot, as.matrix)
  } else if(dim(max_sd_tot[[1]])[1] == 1){
    max_sd_tot <- lapply(max_sd_tot, t)
    corr_tot <- lapply(corr_tot, t)
    modes_paste_tot <- lapply(modes_paste_tot, t)
    dots.int_tot <- lapply(dots.int_tot, t)
  }

  sd_tot_mt <- corr_tot_mt <- dots.int_tot_mt <- modes_paste_tot_mt <- matrix(NA, nrow = dim(result.ploidy)[1], ncol = dim(result.ploidy)[2])
  for(i in 1:dim(result.ploidy)[2]){ # Ind
    for(j in 1:dim(result.ploidy)[1]){ # Chr
      sd_tot_mt[j, i] <- max_sd_tot[[j]][which(ploidys == result.ploidy[j,i]),i]
      corr_tot_mt[j,i] <- corr_tot[[j]][which(ploidys == result.ploidy[j,i]),i]
      modes_paste_tot_mt[j,i] <- modes_paste_tot[[j]][which(ploidys == result.ploidy[j,i]),i]
      dots.int_tot_mt[j,i] <- dots.int_tot[[j]][which(ploidys == result.ploidy[j,i]),i]
    }
  }

  result.ploidy <- t(result.ploidy)
  diff.first.second <- t(diff.first.second)
  sd_tot_mt <- t(sd_tot_mt)
  corr_tot_mt <- t(corr_tot_mt)
  modes_paste_tot_mt <- t(modes_paste_tot_mt)
  dots.int_tot_mt <- t(dots.int_tot_mt)

  colnames(result.ploidy) <- colnames(diff.first.second) <- paste0("Chr",names(by_chr))
  colnames(sd_tot_mt) <- colnames(corr_tot_mt) <- colnames(modes_paste_tot_mt) <- colnames(dots.int_tot_mt) <- paste0("Chr",names(by_chr))

  rownames(diff.first.second) <- rownames(result.ploidy)
  rownames(sd_tot_mt) <- rownames(corr_tot_mt) <- rownames(modes_paste_tot_mt) <- rownames(dots.int_tot_mt) <- rownames(result.ploidy)

  est.ploidy.chr_df <- list(result.ploidy,         # Estimated ploidy by area method
                            dots.int_tot_mt,       # Proportion of dots inside selected area
                            diff.first.second,     # Difference between first and second place in area method
                            sd_tot_mt,             # standard deviation inside area
                            corr_tot_mt,           # Highest correlation
                            modes_paste_tot_mt)    # Modes inside areas
  return(est.ploidy.chr_df)
}

#' Generate plots for overall analysis
#'
#' @param est.ploidy.chr_df list of data.frames with ploidy estimation information
#' @param filter_diff filter by difference on area proportion between first and second place
#' @param filter_corr filter by correlation of observed and expected peak position
#'
#' @import tidyr
#' @import dplyr
#' @importFrom reshape2 melt
#' @import ggplot2
#' @importFrom ggpubr ggarrange
#'
#' @export
plots_overall <- function(est.ploidy.chr_df, filter_diff=0, filter_corr=0){
  idx <- which(apply(est.ploidy.chr_df[[1]], 1, function(x) length(unique(x)) > 1))
  if(length(idx) > 0){
    ploidies <- est.ploidy.chr_df[[1]][-idx,]
    if(is.null(dim(ploidies))) { # If only one individual
      ploidies <- t(ploidies)
      rownames(ploidies) <- rownames(est.ploidy.chr_df[[1]])[-idx]
    }
    ploidies.prop <- est.ploidy.chr_df[[2]][-idx,]
    if(is.null(dim(ploidies.prop))) { # If only one individual
      ploidies.prop <- t(ploidies.prop)
      rownames(ploidies.prop) <- rownames(est.ploidy.chr_df[[1]])[-idx]
    }
    ploidies.sd <- est.ploidy.chr_df[[4]][-idx,]
    if(is.null(dim(ploidies.sd))) { # If only one individual
      ploidies.sd <- t(ploidies.sd)
      rownames(ploidies.sd) <- rownames(est.ploidy.chr_df[[1]])[-idx]
    }
    ploidies.corr <- est.ploidy.chr_df[[5]][-idx,]
    if(is.null(dim(ploidies.corr))) { # If only one individual
      ploidies.corr <- t(ploidies.corr)
      rownames(ploidies.corr) <- rownames(est.ploidy.chr_df[[1]])[-idx]
    }
    ploidies.diff <- est.ploidy.chr_df[[3]][-idx,]
    if(is.null(dim(ploidies.diff))) { # If only one individual
      ploidies.diff <- t(ploidies.diff)
      rownames(ploidies.diff) <- rownames(est.ploidy.chr_df[[1]])[-idx]
    }
    # Graphics using only area method and aneuploidy individuals
    aneuploids <- est.ploidy.chr_df[[1]][idx,]
    aneuploids.prop <- est.ploidy.chr_df[[2]][idx,]
    aneuploids.sd <- est.ploidy.chr_df[[4]][idx,]
    aneuploids.corr <- est.ploidy.chr_df[[5]][idx,]
    aneuploids.diff <- est.ploidy.chr_df[[3]][idx,]

    if(is.vector(aneuploids)){
      aneuploids <- data.frame(Var1 = rownames(est.ploidy.chr_df[[1]])[idx], Var2=names(aneuploids),melt(aneuploids, value.name = "ploidy"))
      aneuploids.prop <- data.frame(Var1 = rownames(est.ploidy.chr_df[[1]])[idx], Var2=names(aneuploids.prop),melt(aneuploids.prop, value.name = "prop"))
      aneuploids.sd <- data.frame(Var1 = rownames(est.ploidy.chr_df[[1]])[idx], Var2=names(aneuploids.sd),melt(aneuploids.sd, value.name = "sd"))
      aneuploids.corr <- data.frame(Var1 = rownames(est.ploidy.chr_df[[1]])[idx], Var2=names(aneuploids.corr),melt(aneuploids.corr, value.name = "corr"))
      aneuploids.diff <- data.frame(Var1 = rownames(est.ploidy.chr_df[[1]])[idx], Var2=names(aneuploids.diff),melt(aneuploids.diff, value.name = "diff"))
    } else {
      aneuploids <- melt(aneuploids, value.name = "ploidy")
      aneuploids.prop <- melt(aneuploids.prop, value.name = "prop")
      aneuploids.sd <- melt(aneuploids.sd, value.name = "sd")
      aneuploids.corr <- melt(aneuploids.corr, value.name = "corr")
      aneuploids.diff <- melt(aneuploids.diff, value.name = "diff")
    }

    # Aneuploid table and samples list
    main_ploidy <- aneuploids %>% group_by(Var1) %>% mutate(main_ploidy = mode(ploidy))

    aneuploidy_df <- main_ploidy[-which(main_ploidy$ploidy == main_ploidy$main_ploidy),]
    colnames(aneuploidy_df) <- c("sample", "chromosome", "aneuploidy_ploidy", "main_ploidy")
    aneuploidy_list <- aneuploidy_df %>% group_by(chromosome, aneuploidy_ploidy, main_ploidy) %>% summarise(samples = paste0(sample, collapse = " "))
    aneuploidy_list <- aneuploidy_list[order(aneuploidy_list$main_ploidy, aneuploidy_list$aneuploidy_ploidy, aneuploidy_list$chromosome),]
    aneuploid_text <- paste0("<br/>Main ploidy - ", aneuploidy_list$main_ploidy, "; ploidy - ", aneuploidy_list$aneuploidy_ploidy, "; in chromosome - ", aneuploidy_list$chromosome,":<br/>", aneuploidy_list$samples)

    aneuploidy_all_weird <- unique(c(as.character(aneuploids.corr$Var1)[which(aneuploids.corr$corr < filter_corr)],
                                     as.character(aneuploids.diff$Var1)[which(aneuploids.diff$diff < filter_diff)]))

    if(length(aneuploidy_all_weird) == 0) aneuploidy_all_weird <- "Set filters"

    # Graphics
    df_n <- melt(table(aneuploids$ploidy), varnames = "ploidy", value.name = "#individuals*#chrom")
    df_n$ploidy <- as.factor(df_n$ploidy)
    p1 <- ggplot(df_n, aes(x = ploidy, y=`#individuals*#chrom`, fill =ploidy)) +
      geom_bar(stat="identity") +
      theme_bw()

    df <- merge(aneuploids, aneuploids.prop, by = c(1,2))
    df <- merge(df, aneuploids.sd, by = c(1,2))
    df <- merge(df, aneuploids.corr, by = c(1,2))
    df <- merge(df, aneuploids.diff, by = c(1,2))

    df <- melt(df, id.vars = c("ploidy", "Var1", "Var2"))

    df$variable <- gsub("sd", "standard deviation", df$variable)
    df$variable <- gsub("prop", "proportion in area", df$variable)
    df$variable <- gsub("corr", "correlation (Pearson)", df$variable)
    df$variable <- gsub("diff", "difference in proportion", df$variable)

    df$ploidy <- as.factor(df$ploidy)
    p2 <- ggplot(df, aes(x = ploidy, y=value, group=ploidy, fill=ploidy)) +
      geom_boxplot() +
      facet_grid(variable~., scales = "free") +
      theme_bw()

    p_aneuploids <- ggarrange(p1, p2, common.legend = T)

    df <- as.data.frame(table(data.frame(chr=aneuploids$Var2, ploidy = aneuploids$ploidy)))
    df <- df %>% group_by(chr, ploidy) %>% summarise(Freq_all = sum(Freq))
    p_aneuploids2 <- ggplot(df, aes(chr,ploidy)) +
      geom_tile(aes(fill = Freq_all)) +
      geom_text(aes(label=Freq_all)) +
      scale_fill_gradient(low="white", high = "red")

  } else {
    ploidies <- est.ploidy.chr_df[[1]]
    if(is.null(dim(ploidies))) { # If only one individual
      ploidies <- t(ploidies)
      rownames(ploidies) <- rownames(est.ploidy.chr_df[[1]])[-idx]
    }
    ploidies.prop <- est.ploidy.chr_df[[2]]
    if(is.null(dim(ploidies.prop))) { # If only one individual
      ploidies.prop <- t(ploidies.prop)
      rownames(ploidies.prop) <- rownames(est.ploidy.chr_df[[1]])[-idx]
    }
    ploidies.sd <- est.ploidy.chr_df[[4]]
    if(is.null(dim(ploidies.sd))) { # If only one individual
      ploidies.sd <- t(ploidies.sd)
      rownames(ploidies.sd) <- rownames(est.ploidy.chr_df[[1]])[-idx]
    }
    ploidies.corr <- est.ploidy.chr_df[[5]]
    if(is.null(dim(ploidies.corr))) { # If only one individual
      ploidies.corr <- t(ploidies.corr)
      rownames(ploidies.corr) <- rownames(est.ploidy.chr_df[[1]])[-idx]
    }
    ploidies.diff <- est.ploidy.chr_df[[3]]
    if(is.null(dim(ploidies.diff))) { # If only one individual
      ploidies.diff <- t(ploidies.diff)
      rownames(ploidies.diff) <- rownames(est.ploidy.chr_df[[1]])[-idx]
    }

    p_aneuploids <- paste0("No aneuploid samples")
    p_aneuploids2 <- paste0("No aneuploid samples")
    aneuploid_text <- paste0("No aneuploid samples")
    aneuploidy_all_weird <- paste0("No aneuploid samples")
    aneuploidy_df <-  paste0("No aneuploid samples")
  }

  # Graphics using only area method and non-aneuploidy individuals
  if(nrow(ploidies) < 2) {
    ploidy_all_tex <- list(rownames(ploidies))
    names(ploidy_all_tex) <- ploidies[,1]
  } else {
  ploidy_all_tex <- split(names(ploidies[,1]), ploidies[,1])
  }
  text_ploidy <- vector()
  for(i in 1:length(ploidy_all_tex)){
    line1 <- paste0("<br/>ploidy ",names(ploidy_all_tex)[i], ":<br/>")
    line2 <- paste0(ploidy_all_tex[[i]], collapse = " ")
    text_ploidy <- paste0(text_ploidy, line1, line2)
  }

  df_n <- melt(table(ploidies[,1]), varnames = "ploidy", value.name = "#individuals")
  df_n$ploidy <- as.factor(df_n$ploidy)
  p1 <- ggplot(df_n, aes(x = ploidy, y=`#individuals`, fill =ploidy)) +
    geom_bar(stat="identity") +
    theme_bw()

  if(is.vector(ploidies)){
    ploidies <- data.frame(Var1 = rownames(est.ploidy.chr_df[[1]])[idx], Var2=names(ploidies),melt(ploidies, value.name = "ploidy"))
    ploidies.prop <- data.frame(Var1 = rownames(est.ploidy.chr_df[[1]])[idx], Var2=names(ploidies.prop),melt(ploidies.prop, value.name = "prop"))
    ploidies.sd <- data.frame(Var1 = rownames(est.ploidy.chr_df[[1]])[idx], Var2=names(ploidies.sd),melt(ploidies.sd, value.name = "sd"))
    ploidies.corr <- data.frame(Var1 = rownames(est.ploidy.chr_df[[1]])[idx], Var2=names(ploidies.corr),melt(ploidies.corr, value.name = "corr"))
    ploidies.diff <- data.frame(Var1 = rownames(est.ploidy.chr_df[[1]])[idx], Var2=names(ploidies.diff),melt(ploidies.diff, value.name = "diff"))
  } else {
    ploidies <- melt(ploidies, value.name = "ploidy")
    ploidies.prop <- melt(ploidies.prop, value.name = "prop")
    ploidies.sd <- melt(ploidies.sd, value.name = "sd")
    ploidies.corr <- melt(ploidies.corr, value.name = "corr")
    ploidies.diff <- melt(ploidies.diff, value.name = "diff")
  }

  df <- merge(ploidies, ploidies.prop, by = c(1,2))
  df <- merge(df, ploidies.sd, by = c(1,2))
  df <- merge(df, ploidies.corr, by = c(1,2))
  df <- merge(df, ploidies.diff, by = c(1,2))

  df <- melt(df, id.vars = c("ploidy", "Var1", "Var2"))

  df$variable <- gsub("sd", "standard deviation", df$variable)
  df$variable <- gsub("prop", "proportion in area", df$variable)
  df$variable <- gsub("corr", "correlation (Pearson)", df$variable)
  df$variable <- gsub("diff", "difference in proportion", df$variable)

  df$ploidy <- as.factor(df$ploidy)
  p2 <- ggplot(df, aes(x = ploidy, y=value, group=ploidy, fill=ploidy)) +
    geom_boxplot() +
    facet_grid(variable~., scales = "free") +
    theme_bw()

  p <- ggarrange(p1, p2, common.legend = T)

  # Weird individuals
  # set filters
  # All
  ploidy_all_weird <- unique(c(as.character(ploidies.corr$Var1)[which(ploidies.corr$corr < filter_corr)],
                               as.character(ploidies.diff$Var1)[which(ploidies.diff$diff < filter_diff)]))

  if(length(ploidy_all_weird) == 0) ploidy_all_weird <- "Set filters"


  list(p, p_aneuploids, p_aneuploids2, text_ploidy, aneuploid_text,
       aneuploidy_df, ploidy_all_weird, aneuploidy_all_weird)
}

#' Plot BAF
#'
#' @param data_sample data.frame with BAF and genomic position information
#' @param area_single area around the expected peak to be considered
#' @param ploidy expected ploidy
#' @param dot.size graphic dot size
#' @param add_lines add expected peaks lines
#' @param colors add area colors
#'
#' @export
plot_baf <- function(data_sample, area_single, ploidy, dot.size, add_lines, colors){
  ymin <- seq(0, 1, 1/ploidy) - (area_single/(ploidy*2))
  ymax <- seq(0, 1, 1/ploidy) + (area_single/(ploidy*2))

  ymin[which(ymin < 0)] <- 0
  ymax[which(ymax > 1)] <- 1
  rets <- data.frame(ymin, ymax, xmax = Inf, xmin = -Inf)

  idx_tot <- FALSE
  idx <- list()
  for(i in 1:nrow(rets)){
    idx <- data_sample$sample >= rets$ymin[i] & data_sample$sample <= rets$ymax[i]
    idx_tot <- idx_tot | idx
  }

  data_sample$color <- NA
  data_sample$color[which(!idx_tot)] <- "red"
  data_sample$color[which(idx_tot)] <- "black"

  p_baf <- data_sample %>% ggplot(aes(x=Position, y=sample)) +
    {if(colors) geom_point(aes(color = color), alpha=0.7, size=dot.size) else geom_point(alpha=0.7, size=dot.size)} +
    scale_color_manual(values = c("red", "black")) +
    {if(add_lines) geom_rect(data = rets, inherit.aes=FALSE,
                             aes(ymin=ymin,ymax=ymax,
                                 xmax = xmax, xmin = xmin,
                                 alpha=0.001, color = "red"))} +
    ylab("BAF") +
    facet_grid(~ Chr, scales = "free_x") + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position="none")
  return(p_baf)
}


plot_baf_hist <- function(data_sample, area_single, ploidy, colors, add_lines){
  ymin <- seq(0, 1, 1/ploidy) - (area_single/(ploidy*2))
  ymax <- seq(0, 1, 1/ploidy) + (area_single/(ploidy*2))

  ymin[which(ymin < 0)] <- 0
  ymax[which(ymax > 1)] <- 1
  rets <- data.frame(ymin, ymax, xmax = Inf, xmin = -Inf)

  data_sample$color <- NA
  data_sample$color[which(data_sample$color == "black")] <- "inside area"
  data_sample$color[which(data_sample$color == "red")] <- "outside area"

  modes.df2 <- list()
  for(i in 1:nrow(rets)) {
    split_data <- split(data_sample$sample, data_sample$Chr)
    x <- split_data[[1]]
    modes.df2[[i]] <- matrix(unlist(sapply(split_data, function(x){
      estimated <- median(x[which(x > rets$ymin[i] & x < rets$ymax[i])])
      expected <- seq(0, 1, 1/ploidy)[i]
      return(data.frame(estimated, expected))
    })), nrow = 2)
    colnames(modes.df2[[i]]) <- names(split_data)
    rownames(modes.df2[[i]]) <- c("modes estimated", "modes expected")
  }

  modes.df2 <- lapply(modes.df2, function(x) {
    y <- t(x)
    y <- data.frame(Chr=rownames(y), y)
    pivot_longer(y, cols = 2:3)
  })

  modes.df2 <- do.call(rbind, modes.df2)
  modes.df2$alpha <- modes.df2$name
  modes.df2$alpha <- gsub("modes.expected", 1, modes.df2$alpha)
  modes.df2$alpha <- as.numeric(gsub("modes.estimated", 0.5, modes.df2$alpha))

  p_hist <- data_sample %>% ggplot(aes(x=sample)) +
    {if(colors) geom_histogram(aes(fill = color)) else  geom_histogram()} +
    scale_x_continuous(breaks = round(seq(0, 1, 1/ploidy),2)) +
    {if(add_lines) geom_vline(data = modes.df2, aes(xintercept= value,
                                                          color = name,
                                                          linetype = name,
                                                          alpha = alpha),
                                    linewidth = 0.8)}+
    {if(colors) scale_fill_manual(values = c("red", "black"))} +
    scale_color_manual(values = c("blue", "purple")) +
    scale_linetype_manual(values = c("dashed", "solid"), guide="none") +
    scale_alpha(range = c(0.7, 1), guide="none") +
    facet_grid(~ Chr, scales = "free_x") + theme_bw() +  xlab("BAF") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position="bottom") +
    labs(color= "Peaks", fill = "Area")

  p_hist_all <- data_sample %>% ggplot(aes(x=sample)) +
    {if(colors) geom_histogram(aes(fill = color)) else  geom_histogram()} +
    scale_x_continuous(breaks = round(seq(0, 1, 1/ploidy),2)) +
    {if(colors) scale_fill_manual(values = c("red", "black"))} +
    scale_color_manual(values = c("blue", "purple")) +
    scale_linetype_manual(values = c("dashed", "solid"), guide="none") +
    scale_alpha(range = c(0.7, 1), guide="none") +
    theme_bw() +  xlab("BAF") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position="bottom") +
    labs(color= "Peaks", fill = "Area")

  p_hist <- ggarrange(p_hist_all, p_hist, widths = c(1,max(data_sample$Chr)-2))
  return(p_hist)
}
