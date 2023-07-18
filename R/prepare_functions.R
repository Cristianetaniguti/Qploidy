#' Edit Axiom Summary file
#'
#' Remove consecutive A allele probe from the summary file
#'
#' @param summary_df data.frame with A and B probes intensities
#'
#' @return list with A and B probes
#'
#' @export
clean_summary <- function(summary_df){
  summaries_filt <- summary_df
  if(length(grep("Contig", summary_df$probeset_id)) > 0) summaries_filt <- summary_df[-grep("Contig", summaries_filt$probeset_id),]
  if(length(grep("contig", summaries_filt$probeset_id)) > 0) summaries_filt <- summaries_filt[-grep("contig", summaries_filt$probeset_id),]
  if(length(grep("comp", summaries_filt$probeset_id)) > 0) summaries_filt <- summaries_filt[-grep("comp", summaries_filt$probeset_id),]

  A_probes <- summaries_filt[seq(1,dim(summaries_filt)[1], 2),]
  B_probes <- summaries_filt[seq(2,dim(summaries_filt)[1], 2),]
  return(list(A_probes=A_probes, B_probes = B_probes))
}

#' Get R and theta values from summary file with option to do standard normalization by plate and markers
#'
#' @param cleaned_summary summary object from clean_summary function. List with A and B intensities.
#' @param ind.names data.frame with individual names correspondence
#' @param sd.normalization logical. If TRUE performs standard normalization
#' @param atan logical. If TRUE calculates the theta using atan2
#'
#' @importFrom tidyr pivot_longer
#' @export
get_R_theta <- function(cleaned_summary, ind.names, sd.normalization = TRUE, atan = FALSE){
  R_all <- cleaned_summary$A_probes[,-1] + cleaned_summary$B_probes[,-1]
  if(atan){
    theta_all <- as.data.frame((atan2(as.matrix(cleaned_summary$B_probes[,-1]), as.matrix(cleaned_summary$A_probes[,-1])))/(pi/2))
  }  else {
    theta_all <- cleaned_summary$B_probes[,-1]/(cleaned_summary$B_probes[,-1] + cleaned_summary$A_probes[,-1])
  }

  probes_names <-  cleaned_summary$A_probes[,1]
  probes_names <- gsub("-A","", probes_names)

  R_all <- cbind(MarkerName = probes_names, R_all)
  theta_all <- cbind(MarkerName = probes_names, theta_all)

  array.id <- colnames(R_all)
  ind.names <- as.data.frame(ind.names)
  ids <- ind.names[,2][match(colnames(R_all), ind.names[,1])]
  colnames(R_all) <- ids

  ids <- ind.names[,2][match(colnames(theta_all), ind.names[,1])]
  colnames(theta_all) <- ids

  colnames(theta_all)[1] <- colnames(R_all)[1] <- "MarkerName"

  # R Quantile normalization
  if(sd.normalization){
    R_melt <- pivot_longer(R_all, cols = -1, names_to = "ind", values_to = "R")

    exp <- sapply(strsplit(array.id[-1], "-"), "[[", 4)
    exp <- sapply(strsplit(exp, "_"), "[[", 1)

    array_ids <- data.frame(array=exp, inds = ids[-1])
    R_melt$array <- array_ids$array[match(R_melt$ind, array_ids$inds)]

    R_z <- R_melt %>% group_by(MarkerName, array) %>% mutate(zscore = (R-mean(R))/sd(R))
    R_all <- pivot_wider(R_z[,c(1,2,5)], names_from = ind, values_from = zscore)
  }

  return(list(R_all=R_all, theta_all = theta_all))
}

#' Convert Summary file to fitpoly format
#'
#' @param R_all data.frame with R values
#' @param theta_all data.frame with theta values
#' @param ind.names data.frame first column with the plate name and second with sample name
#'
#' @importFrom tidyr pivot_longer
#' @export
summary_to_fitpoly <- function(R_all, theta_all){
  Y_all <- R_all[,-1]*theta_all[,-1]
  X_all <- R_all[,-1] - Y_all

  X_all <- cbind(MarkerName= R_all[,1], X_all)
  Y_all <- cbind(MarkerName = R_all[,1], Y_all)

  R <- pivot_longer(R_all, cols = 2:ncol(R_all), values_to = "R", names_to = "SampleName")
  theta <- pivot_longer(theta_all, cols = 2:ncol(theta_all), values_to = "ratio", names_to = "SampleName")
  X <- pivot_longer(X_all, cols = 2:ncol(X_all), values_to = "X", names_to = "SampleName")
  Y <- pivot_longer(Y_all, cols = 2:ncol(Y_all), values_to = "Y", names_to = "SampleName")

  fitpoly_input <- cbind(X, Y = Y[,3], R = R[,3], ratio = theta[,3])

  return(fitpoly_input)
}
