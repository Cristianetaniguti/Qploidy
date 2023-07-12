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


#' Convert Summary file to fitpoly format
#'
#' @param cleaned_summary list with A and B intensities
#' @param ind.names data.frame first column with the plate name and second with sample name
#'
#' @importFrom tidyr pivot_longer
#' @export
summary_to_fitpoly <- function(cleaned_summary, ind.names, geno.pos){
  R_all <- cleaned_summary$A_probes[,-1] + cleaned_summary$B_probes[,-1]
  theta_all <- cleaned_summary$B_probes[,-1]/(cleaned_summary$B_probes[,-1] + cleaned_summary$A_probes[,-1])
  #theta_all <- as.data.frame((atan2(as.matrix(cleaned_summary$B_probes[,-1]), as.matrix(cleaned_summary$A_probes[,-1])))/(pi/2))
  probes_names <-  cleaned_summary$A_probes[,1]
  probes_names <- gsub("-A","", as.vector(probes_names$probeset_id))

  X_all <- cbind(MarkerName = probes_names, cleaned_summary$A_probes[,-1])
  Y_all <- cbind(MarkerName = probes_names, cleaned_summary$B_probes[,-1])

  R_all <- cbind(MarkerName = probes_names, R_all)
  theta_all <- cbind(MarkerName = probes_names, theta_all)

  R <- pivot_longer(R_all, cols = 2:ncol(R_all), values_to = "R", names_to = "SampleName")
  theta <- pivot_longer(theta_all, cols = 2:ncol(theta_all), values_to = "ratio", names_to = "SampleName")
  X <- pivot_longer(X_all, cols = 2:ncol(X_all), values_to = "X", names_to = "SampleName")
  Y <- pivot_longer(Y_all, cols = 2:ncol(Y_all), values_to = "Y", names_to = "SampleName")

  fitpoly_input <- cbind(X, Y = Y[,3], R = R_all[,3], ratio = theta_all[,3])

  ind.names <- as.data.frame(ind.names)
  ids <- ind.names[,2][match(fitpoly_input$SampleName, ind.names[,1])]
  fitpoly_input$SampleName <- ids

  ids <- ind.names[,2][match(colnames(R_all), ind.names[,1])]
  colnames(R_all) <- ids

  ids <- ind.names[,2][match(colnames(theta_all), ind.names[,1])]
  colnames(theta_all) <- ids

  colnames(theta_all)[1] <- colnames(R_all)[1] <- "MarkerName"

  return(list(fitpoly_input = fitpoly_input, R= R_all, theta = theta_all))
}
