globalVariables(c("color", "xmax", "xmin", "value", "name",
                  "Var1", "ploidy", "chromosome", "aneuploidy_ploidy",
                  "#individuals*#chrom", "Freq", "Freq_all", "#individuals",
                  "LG", "MarkerName", "R", "SampleName", "baf", "sd", "cor",
                  "rnom", "X", "Y", "runif", "prob"))

#' Estimate ploidy using area method
#'
#' This function estimates ploidy using the area method. It evaluates the
#' number of copies by chromosome, sample, or chromosome arm. Note that this
#' function does not have optimal performance, and visual inspection of the
#' plots is required to confirm the estimated ploidy.
#'
#' @param qploidy_standardization Object of class qploidy_standardization.
#' @param samples If "all", all samples contained in the qploidy_standardization
#' object will be evaluated. If a vector with sample names is provided, only
#' those will be evaluated.
#' @param level Character identifying the level of the analysis. Must be one of
#' "chromosome", "sample", or "chromosome-arm". If `chromosome-arm`, the
#' analysis will be performed by chromosome arm (only if `centromeres` argument
#' is defined).
#' @param ploidies Vector of ploidy levels to test. This parameter must be
#' defined.
#' @param area Area around the expected peak to be considered. Default is 0.75.
#' @param centromeres Vector with centromere genomic positions in bp. The vector
#' should be named with the chromosome IDs. This information will only be used
#' if `chromosome-arm` level is defined.
#'
#' @return A list of class `qploidy_area_ploidy_estimation` containing:
#' \itemize{
#'   \item \code{ploidy}: Estimated ploidy by area method.
#'   \item \code{prop_inside_area}: Proportion of dots inside selected area.
#'   \item \code{diff_first_second}: Difference between first and second place in area method.
#'   \item \code{sd_inside_area}: Standard deviation inside area.
#'   \item \code{highest_correlation_modes}: Highest correlation.
#'   \item \code{modes_inside_area}: Modes inside areas.
#'   \item \code{tested}: Tested ploidies.
#'   \item \code{ploidy.sep}: Separated ploidy results.
#'   \item \code{chr}: Unique chromosomes in the dataset.
#'   \item \code{n.inbred}: Number of highly inbred samples.
#' }
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
                                 centromeres = NULL) {

  # Input checks
  if (!level %in% c("sample", "chromosome", "chromosome-arm")) {
    stop("Invalid level. Must be one of 'sample', 'chromosome', or 'chromosome-arm'.")
  }

  if (is.null(ploidies)) {
    stop("The parameter 'ploidies' must be defined.")
  }

  if (!is.null(centromeres)) {
    if (!is.vector(centromeres) || is.null(names(centromeres))) {
      stop("The 'centromeres' parameter must be a named vector with chromosome IDs as names.")
    }

    chr_ids <- unique(qploidy_standardization$data$Chr)
    if (!all(names(centromeres) %in% chr_ids)) {
      stop("Chromosome names in 'centromeres' vector do not match the ones in the dataset.")
    }
  }

  # Filter data based on samples
  if (all(samples == "all")) {
    data_sample <- qploidy_standardization$data[, c(9, 10, 2, 8)] %>%
      filter(!is.na(baf))
  } else {
    data_sample <- qploidy_standardization$data[, c(9, 10, 2, 8)] %>%
      filter(SampleName %in% samples & !is.na(baf))
  }

  if (nrow(data_sample) == 0) stop("No data available for the selected samples.")

  # Check for duplicated marker positions
  one_sample <- data_sample$Position[which(data_sample$SampleName %in%
                                           data_sample$SampleName[1])]
  if (any(duplicated(one_sample))) {
    warning("There are duplicated marker positions. Only the first one will be kept.")
    data_sample <- data_sample %>% group_by(SampleName) %>%
      filter(!duplicated(Position))
  }

  data_sample <- pivot_wider(data_sample, names_from = SampleName,
                             values_from = baf)

  if (is.null(data_sample)) stop("Define data_sample argument.")

  # Split data by level
  if (level == "chromosome" || level == "chromosome-arm") {
    by_chr <- split(data_sample, data_sample$Chr)

    if (!is.null(centromeres) && level == "chromosome-arm") {
      chrs <- match(names(centromeres), names(by_chr))
      if (length(chrs) == 0) stop("Chromosome names in centromeres vector do not match the ones in dataset.")
      for (i in seq_along(chrs)) {
        idx <- which(by_chr[[chrs[i]]]$Position <= centromeres[i])
        if(length(idx) == 0) stop("Verify your centromere position. No data was found for one of the arms.")
        by_chr[[chrs[i]]]$Chr[idx] <- paste0(by_chr[[chrs[i]]]$Chr[idx], ".1")
        idx <- which(by_chr[[chrs[i]]]$Position > centromeres[i])
        if(length(idx) == 0) stop("Verify your centromere position. No data was found for one of the arms.")
        by_chr[[chrs[i]]]$Chr[idx] <- paste0(by_chr[[chrs[i]]]$Chr[idx], ".2")
        by_chr[[chrs[i]]] <- split(by_chr[[chrs[i]]], by_chr[[chrs[i]]]$Chr)
      }
      names(by_chr) <- NULL
      skip <- which(sapply(by_chr, length) != 2)
      if (length(skip) != 0) {
        for (i in seq_along(skip)) {
          by_chr[[i]] <- list(by_chr[[i]])
          names(by_chr[[i]]) <- unique(by_chr[[i]][[1]]$Chr)
        }
      }

      by_chr <- unlist(by_chr, recursive = F)
    }
  } else if (level == "sample") {
    by_chr <- list(data_sample)
  }

  freq <- pascalTriangle(max(ploidies))
  freq <- freq[-1]
  ploidies <- unique(c(1, ploidies))
  dots.int_tot <- corr_tot <- max_sd_tot <- modes_paste_tot <-
    as.list(rep(NA, length(by_chr)))
  for (p in ploidies) {
    # Area method
    ymin <- seq(0, 1, 1 / p) - (area / (p * 2))
    ymax <- seq(0, 1, 1 / p) + (area / (p * 2))

    ymin[which(ymin < 0)] <- 0
    ymax[which(ymax > 1)] <- 1
    rets <- data.frame(ymin, ymax)

    prop_tot <- modes_tot <- sd_areas_tot <- as.list(rep(NA, length(by_chr)))
    for (z in seq_along(by_chr)) {
      for (i in seq_len(nrow(rets))) {
        prop <- apply(by_chr[[z]][, -c(1, 2)], 2, function(x)
          sum(x >= rets$ymin[i] & x <= rets$ymax[i], na.rm = TRUE))
        modes <- apply(by_chr[[z]][, -c(1, 2)], 2, function(x)
          mode(x[x >= rets$ymin[i] & x <= rets$ymax[i]]))
        sd_areas <- apply(by_chr[[z]][, -c(1, 2)], 2, function(x)
          sd(x[x >= rets$ymin[i] & x <= rets$ymax[i]], na.rm = TRUE))

        prop_tot[[z]] <- rbind(prop_tot[[z]], prop)
        modes_tot[[z]] <- rbind(modes_tot[[z]], modes)
        sd_areas_tot[[z]] <- rbind(sd_areas_tot[[z]], sd_areas)
      }
    }

    # remove NAs
    prop_tot <- lapply(prop_tot, function(x) x[-1, ])
    for (z in seq_along(prop_tot)) {
      if (is.null(ncol(prop_tot[[z]]))) dots.int <- sum(prop_tot[[z]],
        na.rm = TRUE) / dim(by_chr[[z]])[1] else
        dots.int <- apply(prop_tot[[z]], 2, function(x) sum(x, na.rm = TRUE) /
          dim(by_chr[[z]])[1])
      dots.int_tot[[z]] <- rbind(dots.int_tot[[z]], dots.int)

      modes_paste <- apply(modes_tot[[z]], 2, function(x)
        paste0(round(x[-1], 3), collapse = "/"))
      modes_paste_tot[[z]] <- rbind(modes_paste_tot[[z]], modes_paste)

      corr <- apply(modes_tot[[z]], 2, function(x)
        cor(x = x[-1], y = seq(0, 1, 1 / p)))
      corr_tot[[z]] <- rbind(corr_tot[[z]], corr)

      max_sd <- apply(sd_areas_tot[[z]], 2, function(x) max(x[-1]))
      max_sd_tot[[z]] <- rbind(max_sd_tot[[z]], max_sd)
    }
  }

  dots.int_tot <- lapply(dots.int_tot, function(x) x[-1, ])
  modes_paste_tot <- lapply(modes_paste_tot, function(x) x[-1, ])
  corr_tot <- lapply(corr_tot, function(x) x[-1, ])
  max_sd_tot <- lapply(max_sd_tot, function(x) x[-1, ])

  # Area method
  result.ploidy <- diff.count <- second <- diff.second <-
    diff.first.second <- list()
  for (z in seq_along(dots.int_tot)) {
    if (is.null(rownames(dots.int_tot[[z]]))) {
      names(dots.int_tot[[z]]) <- ploidies
      result.ploidy[[z]] <- ploidies[order(dots.int_tot[[z]],
        decreasing = T)][1]

      diff.count[[z]] <- dots.int_tot[[z]][order(dots.int_tot[[z]],
        decreasing = T)][1]
      second[[z]] <- ploidies[order(dots.int_tot[[z]], decreasing = T)][2]
      diff.second[[z]] <- dots.int_tot[[z]][order(dots.int_tot[[z]],
        decreasing = T)][2]
    } else {
      rownames(dots.int_tot[[z]]) <- ploidies
      result.ploidy[[z]] <- apply(dots.int_tot[[z]], 2, function(x)
        ploidies[order(x, decreasing = T)][1])
      diff.count[[z]] <- apply(dots.int_tot[[z]], 2, function(x)
        x[order(x, decreasing = T)][1])
      second[[z]] <- apply(dots.int_tot[[z]], 2, function(x)
        ploidies[order(x, decreasing = T)][2])
      diff.second[[z]] <- apply(dots.int_tot[[z]], 2, function(x)
        x[order(x, decreasing = T)][2])
    }
    diff.first.second[[z]] <- diff.count[[z]] - diff.second[[z]]
  }

  result.ploidy <- do.call(rbind, result.ploidy)
  diff.first.second <- do.call(rbind, diff.first.second)

  if (is.vector(max_sd_tot[[1]])) { # only one individual selected
    max_sd_tot <- lapply(max_sd_tot, as.matrix)
    corr_tot <- lapply(corr_tot, as.matrix)
    modes_paste_tot <- lapply(modes_paste_tot, as.matrix)
    dots.int_tot <- lapply(dots.int_tot, as.matrix)
  } else if (dim(max_sd_tot[[1]])[1] == 1) {
    max_sd_tot <- lapply(max_sd_tot, t)
    corr_tot <- lapply(corr_tot, t)
    modes_paste_tot <- lapply(modes_paste_tot, t)
    dots.int_tot <- lapply(dots.int_tot, t)
  }

  sd_tot_mt <- corr_tot_mt <- dots.int_tot_mt <- modes_paste_tot_mt <-
    matrix(NA, nrow = dim(result.ploidy)[1], ncol = dim(result.ploidy)[2])
  for (i in seq_len(dim(result.ploidy)[2])) { # Ind
    for (j in seq_len(dim(result.ploidy)[1])) { # Chr
      sd_tot_mt[j, i] <- max_sd_tot[[j]][which(ploidies == result.ploidy[j, i]), i]
      corr_tot_mt[j, i] <- corr_tot[[j]][which(ploidies == result.ploidy[j, i]), i]
      modes_paste_tot_mt[j, i] <- modes_paste_tot[[j]][which(ploidies == result.ploidy[j, i]), i]
      dots.int_tot_mt[j, i] <- dots.int_tot[[j]][which(ploidies == result.ploidy[j, i]), i]
    }
  }

  result.ploidy <- t(result.ploidy)
  diff.first.second <- t(diff.first.second)
  sd_tot_mt <- t(sd_tot_mt)
  corr_tot_mt <- t(corr_tot_mt)
  modes_paste_tot_mt <- t(modes_paste_tot_mt)
  dots.int_tot_mt <- t(dots.int_tot_mt)

  colnames(result.ploidy) <- colnames(diff.first.second) <- names(by_chr)
  colnames(sd_tot_mt) <- colnames(corr_tot_mt) <- colnames(modes_paste_tot_mt) <-
    colnames(dots.int_tot_mt) <- names(by_chr)

  rownames(diff.first.second) <- rownames(result.ploidy)
  rownames(sd_tot_mt) <- rownames(corr_tot_mt) <- rownames(modes_paste_tot_mt) <-
    rownames(dots.int_tot_mt) <- rownames(result.ploidy)

  # Find homozygous
  idx <- which(result.ploidy == 1)
  if (length(idx) > 0) result.ploidy[idx] <- NA

  if (level == "sample") n.inbred <- length(idx) else {
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
#' @param filter_diff filter by difference on area proportion between first and
#' second place
#'
#' @return An updated object of class `qploidy_area_ploidy_estimation` with the following modifications:
#'
#' - `ploidy`: A matrix where chromosome-arm level results are merged into chromosome-level format.
#'    If `filter_diff` is provided, ploidy values with differences below the threshold are set to `NA`.
#'
#' The structure of the returned object remains consistent with the input, but with updated ploidy information.
#'
#' @export
merge_arms_format <- function(x, filter_diff = NULL){

  ploidy <- x$ploidy
  if(!is.null(filter_diff)){
    ploidy[which(x$diff_first_second < filter_diff)] <- NA
  }

  chr <- sapply(strsplit(colnames(ploidy), "[.]"), "[[", 1)
  result.ploidy.up <- vector()
  for (i in seq_along(unique(chr))) {
    if (!is.null(dim(ploidy[, which(chr == unique(chr)[i])])[1])) {
      new.col <- apply(ploidy[, which(chr == unique(chr)[i])], 1,
        function(x) if (length(unique(x)) == 1) unique(x) else
          paste0(x, collapse = "/"))
    } else new.col <- ploidy[, which(chr == unique(chr)[i])]
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
#' @return No return value, called for side effects.
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
                          {if(!is.null(colnames(x$ploidy)))
                            paste0(colnames(x$ploidy), collapse = ",") else
                            paste0(x$chr, collapse = ",")},
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
#' @return A logical vector where each element corresponds to an individual in the
#'         input ploidy table. The value is `TRUE` if the individual is identified
#'         as potentially aneuploid, and `FALSE` otherwise.
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

