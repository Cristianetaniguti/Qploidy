globalVariables(c("color", "xmax", "xmin", "value", "name",
                  "Var1", "ploidy", "chromosome", "aneuploidy_ploidy",
                  "#individuals*#chrom", "Freq", "Freq_all", "#individuals",
                  "LG", "MarkerName", "R", "SampleName", "baf", "sd"))

#' Estimate ploidy using area method
#'
#' @param qploidy_standardization object of class qploidy_standardization
#' @param samples if "all" all samples contained in the qploidy_standardization object will be evaluate. If vector with sample names is provided, only those will be evaluated.
#' @param level character identifying the level of the analysis. If `chromosome` the number of copies will be estimated by chromosome, if `sample` it will be by sample; if `chromosome-arm` it will be peformed by chromosome arm (only if `centromeres` argument is defined.)
#' @param ploidies range of ploidies to by tested
#' @param area area around the expected peak to be considered
#' @param centromeres vector with centromeres genomic position in bp. The vector should be named with the chromosomes ID. The information will only be used if `chromosome-arm` level is defined.
#'
#' @import tidyr
#' @import dplyr
#'
#' @export
area_estimate_ploidy <- function(qploidy_standardization = NULL,
                                 samples = "all",
                                 level = "chromosome",
                                 ploidies = NULL,
                                 area = 0.75,
                                 centromeres= NULL){

  if(all(samples == "all")) {
    data_sample <- qploidy_standardization$data[,c(9,10,2,8)] %>% filter(!is.na(baf))
  } else {
    data_sample <- qploidy_standardization$data[,c(9,10,2,8)] %>% filter(SampleName %in% samples & !is.na(baf))
  }

  one_sample <- data_sample$Position[which(data_sample$SampleName %in% data_sample$SampleName[1])]
  if(any(duplicated(one_sample))) {
    warning("There are duplicated marker positions. Only the first one will be kept.")
    data_sample <- data_sample %>% group_by(SampleName) %>% filter(!duplicated(Position))
  }

  data_sample <- pivot_wider(data_sample, names_from = SampleName, values_from = baf)

  if(is.null(data_sample) | is.null(ploidies)) stop("Define data_sample and ploidies arguments")

  if(level == "chromosome" | level == "chromosome-arm"){
    by_chr <- split(data_sample, data_sample$Chr)

    if(!is.null(centromeres) & level == "chromosome-arm"){
      chrs <- match(names(centromeres), names(by_chr))
      if(length(chrs) == 0) stop("Chromosome names in centromeres vector do not match the ones in dataset.")
      for(i in 1:length(chrs)){
        idx <- which(by_chr[[chrs[i]]]$Position <= centromeres[i])
        by_chr[[chrs[i]]]$Chr[idx] <- paste0(by_chr[[chrs[i]]]$Chr[idx],".1")
        idx <- which(by_chr[[chrs[i]]]$Position > centromeres[i])
        by_chr[[chrs[i]]]$Chr[idx] <- paste0(by_chr[[chrs[i]]]$Chr[idx],".2")
        by_chr[[chrs[i]]] <- split(by_chr[[chrs[i]]], by_chr[[chrs[i]]]$Chr)
      }
      names(by_chr) <- NULL
      skip <- which(sapply(by_chr, length) != 2)
      if(length(skip) != 0){
        for(i in 1:length(skip)){
          by_chr[[i]] <- list(by_chr[[i]])
          names(by_chr[[i]]) <- unique(by_chr[[i]][[1]]$Chr)
        }
      }

      by_chr <- unlist(by_chr, recursive = F)
    }
  } else if(level == "sample"){
    by_chr <- list(data_sample)
  }

  freq <- pascalTriangle(ploidies[2])
  freq <- freq[-1]
  ploidies <- ploidies[1]:ploidies[2]
  ploidies <- unique(c(1, ploidies))
  dots.int_tot <- corr_tot <- max_sd_tot <- modes_paste_tot <- as.list(rep(NA, length(by_chr)))
  for(j in 1:length(ploidies)){
    # Area method
    ymin <- seq(0, 1, 1/ploidies[j]) - (area/(ploidies[j]*2))
    ymax <- seq(0, 1, 1/ploidies[j]) + (area/(ploidies[j]*2))

    ymin[which(ymin < 0)] <- 0
    ymax[which(ymax > 1)] <- 1
    rets <- data.frame(ymin, ymax)

    prop_tot <- modes_tot  <- sd_areas_tot <- as.list(rep(NA, length(by_chr)))
    for(z in 1:length(by_chr)){
      for(i in 1:nrow(rets)){
        prop <- apply(by_chr[[z]][,-c(1,2)], 2, function(x) sum(x >= rets$ymin[i] & x <= rets$ymax[i], na.rm = TRUE))
        modes <- apply(by_chr[[z]][,-c(1,2)], 2, function(x) mode(x[x >= rets$ymin[i] & x <= rets$ymax[i]]))
        sd_areas <- apply(by_chr[[z]][,-c(1,2)], 2, function(x) sd(x[x >= rets$ymin[i] & x <= rets$ymax[i]], na.rm = TRUE))

        prop_tot[[z]] <- rbind(prop_tot[[z]], prop)
        modes_tot[[z]] <- rbind(modes_tot[[z]], modes)
        sd_areas_tot[[z]] <- rbind(sd_areas_tot[[z]], sd_areas)
      }
    }

    # remove NAs
    prop_tot <- lapply(prop_tot, function(x) x[-1,])
    for(z in 1:length(prop_tot)){
      if(is.null(ncol(prop_tot[[z]]))) dots.int <- sum(prop_tot[[z]], na.rm = TRUE)/dim(by_chr[[z]])[1] else
        dots.int <- apply(prop_tot[[z]], 2, function(x) sum(x, na.rm = TRUE)/dim(by_chr[[z]])[1])
      dots.int_tot[[z]] <- rbind(dots.int_tot[[z]], dots.int)

      modes_paste <- apply(modes_tot[[z]], 2, function(x) paste0(round(x[-1],3), collapse = "/"))
      modes_paste_tot[[z]] <- rbind(modes_paste_tot[[z]], modes_paste)

      corr <- apply(modes_tot[[z]], 2, function(x) cor(x = x[-1], y = seq(0,1,1/ploidies[j])))
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
      names(dots.int_tot[[z]]) <- ploidies
      result.ploidy[[z]] <- ploidies[order(dots.int_tot[[z]], decreasing =T)][1]

      diff.count[[z]] <- dots.int_tot[[z]][order(dots.int_tot[[z]], decreasing =T)][1]
      second[[z]] <- ploidies[order(dots.int_tot[[z]], decreasing =T)][2]
      diff.second[[z]] <- dots.int_tot[[z]][order(dots.int_tot[[z]], decreasing =T)][2]
    } else {
      rownames(dots.int_tot[[z]]) <- ploidies
      result.ploidy[[z]] <- apply(dots.int_tot[[z]], 2, function(x) ploidies[order(x, decreasing =T)][1])
      diff.count[[z]] <- apply(dots.int_tot[[z]], 2, function(x) x[order(x, decreasing =T)][1])
      second[[z]] <- apply(dots.int_tot[[z]], 2, function(x) ploidies[order(x, decreasing =T)][2])
      diff.second[[z]] <- apply(dots.int_tot[[z]], 2, function(x) x[order(x, decreasing =T)][2])
    }
    diff.first.second[[z]] <- diff.count[[z]] - diff.second[[z]]
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
      sd_tot_mt[j, i] <- max_sd_tot[[j]][which(ploidies == result.ploidy[j,i]),i]
      corr_tot_mt[j,i] <- corr_tot[[j]][which(ploidies == result.ploidy[j,i]),i]
      modes_paste_tot_mt[j,i] <- modes_paste_tot[[j]][which(ploidies == result.ploidy[j,i]),i]
      dots.int_tot_mt[j,i] <- dots.int_tot[[j]][which(ploidies == result.ploidy[j,i]),i]
    }
  }

  result.ploidy <- t(result.ploidy)
  diff.first.second <- t(diff.first.second)
  sd_tot_mt <- t(sd_tot_mt)
  corr_tot_mt <- t(corr_tot_mt)
  modes_paste_tot_mt <- t(modes_paste_tot_mt)
  dots.int_tot_mt <- t(dots.int_tot_mt)

  colnames(result.ploidy) <- colnames(diff.first.second) <- names(by_chr)
  colnames(sd_tot_mt) <- colnames(corr_tot_mt) <- colnames(modes_paste_tot_mt) <- colnames(dots.int_tot_mt) <- names(by_chr)

  rownames(diff.first.second) <- rownames(result.ploidy)
  rownames(sd_tot_mt) <- rownames(corr_tot_mt) <- rownames(modes_paste_tot_mt) <- rownames(dots.int_tot_mt) <- rownames(result.ploidy)

  # Find homozygous
  idx <- which(result.ploidy == 1)
  if(length(idx) > 0) result.ploidy[idx] <- NA

  if(level == "sample") n.inbred <- length(idx) else {
    n.inbred <- sum(apply(result.ploidy, 1, function(x) all(is.na(x))))
  }

  est.ploidy.chr_df <- list(ploidy = result.ploidy,                  # Estimated ploidy by area method
                            prop_inside_area = dots.int_tot_mt,      # Proportion of dots inside selected area
                            diff_first_second = diff.first.second,   # Difference between first and second place in area method
                            sd_inside_area = sd_tot_mt,              # standard deviation inside area
                            highest_correlation_modes = corr_tot_mt, # Highest correlation
                            modes_inside_area = modes_paste_tot_mt,
                            tested = ploidies,
                            ploidy.sep = result.ploidy,
                            chr = unique(data_sample$Chr),
                            n.inbred = n.inbred)  # Modes inside areas

  return(structure(est.ploidy.chr_df, class = "qploidy_area_ploidy_estimation"))
}

#' Merges chromosome-arm level analysis results into chromosome level format
#'
#' @param x object of class qploidy_area_ploidy_estimation
#' @param filter_diff filter by difference on area proportion between first and second place
#'
#' @export
merge_arms_format <- function(x, filter_diff = NULL){

  ploidy <- x$ploidy
  if(!is.null(filter_diff)){
    ploidy[which(x$diff_first_second < filter_diff)] <- NA
  }

  chr <- sapply(strsplit(colnames(ploidy), "[.]"), "[[", 1)
  result.ploidy.up <- vector()
  for(i in 1:length(unique(chr))){
    if(!is.null(dim(ploidy[,which(chr == unique(chr)[i])])[1])){
      new.col <- apply(ploidy[,which(chr == unique(chr)[i])],1, function(x) if(length(unique(x)) == 1)unique(x) else paste0(x, collapse = "/"))
    } else new.col <-  ploidy[,which(chr == unique(chr)[i])]
    result.ploidy.up <- cbind(result.ploidy.up, new.col)
  }
  colnames(result.ploidy.up) <- unique(chr)

  x_up <- x
  x_up$ploidy <- result.ploidy.up

  return(x_up)
}



#' print qploidy_area_ploidy_estimation object
#'
#' @param x qploidy_area_ploidy_estimation object
#' @param ... print parameters
#'
#' @method print qploidy_area_ploidy_estimation
#'
#' @export
print.qploidy_area_ploidy_estimation <- function(x, ...){

  count_aneu <- get_aneuploids(x$ploidy)

  df <- data.frame(c1 = c("Number of samples:",
                          "Chromosomes:",
                          "Tested ploidies:",
                          "Number of euploid samples:",
                          "Number of potential aneuploid samples:",
                          "Number of highly inbred samples:"),
                   c2 = c(dim(x$ploidy)[1],
                          {if(!is.null(colnames(x$ploidy))) paste0(colnames(x$ploidy), collapse = ",") else paste0(x$chr, collapse = ",")},
                          paste0(x$tested, collapse = ","),
                          sum(!count_aneu, na.rm = TRUE),
                          sum(count_aneu, na.rm = TRUE),
                          x$n.inbred
                   ))

  colnames(df) <- rownames(df) <- NULL

  cat("Object of class qploidy_area_ploidy_estimation")
  print(format(df, justify = "left", digits = 2))
}


#' indexes for aneuploids
#'
#' @param ploidy_df ploidy table (chromosome in columns and individuals in rows)
#'
#' @export
get_aneuploids  <- function(ploidy_df){

  if(any(grepl("/NA", ploidy_df) | grepl("NA/",ploidy_df))){
    ploidy_df <- gsub("/NA", "", ploidy_df)
    ploidy_df <- gsub("NA/", "", ploidy_df)
  }

  count_aneu <- !apply(ploidy_df, 1, function(y) {
    if(any(is.na(y))) {
      temp <- length(unique(y[-which(is.na(y))]))
      if(temp == 1) TRUE else if(temp == 0) NA else FALSE
    } else length(unique(y)) == 1
  })
  return(count_aneu)
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
