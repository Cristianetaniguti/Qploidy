
#' Estimate Centers for Standardization Using Updog Bias
#'
#' This function calculates the centers for standardization based on the estimated
#' bias from the `updog` package. It identifies genotype dosage clusters and determines
#' whether markers should be retained or removed based on the number of clusters.
#'
#' @param multidog_obj An object of class `multidog` (from the `updog` package),
#' containing information about SNPs, ploidy, sequencing error rates, and bias.
#' @param threshold.n.clusters An integer specifying the minimum number of dosage
#' clusters (heterozygous classes) required for a marker to be retained for
#' standardization. Default is `2`.
#' @param rm.mks A logical vector indicating which markers should be removed.
#' The names of the vector correspond to the marker names.
#'
#' @return A named list where each element corresponds to a marker and contains:
#'   - `rm`: An integer flag indicating whether the marker is retained (`0`) or removed (`1`).
#'   - `centers_theta`: A numeric vector of cluster centers (sorted in descending order).
#'   - `MarkerName`: The name of the marker.
#'   - `n.clusters`: The number of clusters identified for the marker.
#'
#' @details The function uses the `xi_fun` to calculate the cluster centers for each marker
#' based on the ploidy, sequencing error rate, and bias. Markers with fewer clusters than
#' the specified threshold are flagged for removal.
#'
#' @import tidyr
#' @import dplyr
#' @importFrom magrittr "%>%"
#'
#' @export
updog_centers <- function(multidog_obj, threshold.n.clusters=2, rm.mks){

  n.clusters.df <- multidog_obj$inddf %>%
    filter(!(snp %in% rm.mks)) %>%
    group_by(snp) %>%
    summarize(n.clusters = length(table(geno)))

  if(length(rm.mks) >0) snpdf <- multidog_obj$snpdf[-which(multidog_obj$snpdf$snp %in% rm.mks),] else snpdf <- multidog_obj$snpdf
  ploidy <- snpdf$ploidy
  seq <- snpdf$seq
  bias <- snpdf$bias

  result <- list()
  for(i in seq_along(bias)){
    centers_theta <- xi_fun(p = (0:ploidy[i])/ploidy[i], eps = seq[i], h = bias[i])

    result[[i]] <- list(rm = if(n.clusters.df[i,]$n.clusters >= threshold.n.clusters) 0 else 1,
                        centers_theta = sort(1 - centers_theta),
                        MarkerName = n.clusters.df[i,]$snp,
                        n.clusters = n.clusters.df[i,]$n.clusters)
  }

  names(result) <- snpdf$snp
  return(result)
}

#' Calculate B-Allele Frequency (BAF) from Theta Values
#'
#' This function calculates the B-allele frequency (BAF) from normalized theta values,
#' using cluster centers that represent genotype classes. BAF is computed by linearly
#' interpolating the theta values between adjacent genotype cluster centroids.
#'
#' The approach is based on the methodology described by Wang et al. (2007), and is
#' commonly used in SNP genotyping to infer allele-specific signal intensities.
#'
#' @param theta_subject A numeric vector of theta values to be standardized. These typically
#' represent allelic ratios or normalized intensity values for a set of samples.
#'
#' @param centers_theta A numeric vector of length `ploidy + 1`, representing the estimated
#' cluster centers (centroids) for each genotype class. These values should be sorted in
#' increasing order from homozygous reference to homozygous alternative.
#'
#' @param ploidy An integer indicating the ploidy level of the organism (e.g., `2` for diploid).
#'
#' @return A numeric vector of BAF values ranging from 0 to 1
#'
#' @note The `centers_theta` vector must contain exactly `ploidy + 1` values, and must be
#' sorted in ascending order. If `theta_subject` values fall outside the range, BAFs are
#' capped at 0 or 1 accordingly.
#'
#' @references Wang, K., Li, M., Hadley, D., Liu, R., Glessner, J., Grant, S. F. A., Hakonarson, H., & Bucan, M. (2007). PennCNV: An integrated hidden Markov model designed for high-resolution copy number variation detection in whole-genome SNP genotyping data. \emph{Genome Research, 17}(11), 1665–1674. \doi{10.1101/gr.6861907}
#'
#' @examples
#' theta <- c(0.1, 0.35, 0.6, 0.95)
#' centers <- c(0.1, 0.5, 0.9)
#' get_baf(theta, centers, ploidy = 2)
#'
#' @export
get_baf <- function(theta_subject, centers_theta, ploidy){

  if (length(centers_theta) != (ploidy + 1)) {
    stop("centers_theta must contain exactly ploidy + 1 values.")
  }

  baf <- rep(NA, length(theta_subject))
  ploidy_freq <- seq(0,1,1/(ploidy))
  ploidy_freq_multi <- 1/ploidy
  centers_theta <- sort(centers_theta, decreasing = F)
  idx <- which(theta_subject <= centers_theta[1])
  if(length(idx) > 0) baf[idx] <- 0
  for(i in 2:(ploidy+1)){
    idx <- which(theta_subject > centers_theta[i-1] & theta_subject <= centers_theta[i])
    D2 <- centers_theta[i] - centers_theta[i-1]
    D1 <- as.numeric(theta_subject[idx]) - centers_theta[i-1]
    baf[idx] <- ploidy_freq[i-1] + (D1/D2)*ploidy_freq_multi
  }
  idx <- which(theta_subject >= centers_theta[ploidy+1])
  if(length(idx) > 0) baf[idx] <- 1

  return(baf)
}

#' To create baf in parallel
#'
#' @param par_all_item list containing R and theta matrices, and clusters models
#' @param ploidy integer defining ploidy
#'
#' @return A list of numeric vectors, where each vector contains the BAF values
#'         for a corresponding row in the input `par_all_item` matrices. Each BAF
#'         vector has values ranging from 0 to 1, representing the standardized
#'         allelic ratios for the respective samples or markers.
#'
#' @export
get_baf_par <- function(par_all_item, ploidy=2){
  baf <- list()
  i <- 1
  for(i in seq_len(nrow(par_all_item[[1]]))){
    baf[[i]] <- get_baf(theta_subject = as.numeric(par_all_item[[1]][i,-1]),
                        centers_theta = as.numeric(par_all_item[[2]][i,-1]),
                        ploidy = ploidy)
  }

  return(baf)
}

#' Estimate Cluster Centers for Genotype Dosage Classes
#'
#' This function estimates the cluster centers for each genotype dosage class
#' based on the `theta` values (e.g., allelic ratios or normalized signal intensities).
#' It supports imputing missing clusters and optionally removing outliers.
#'
#' @param ratio_geno A data.frame containing the following columns:
#'   - `MarkerName`: Identifier for each marker.
#'   - `SampleName`: Identifier for each sample.
#'   - `theta`: Numeric variable representing allelic ratio or signal intensity.
#'   - `geno`: Integer dosage (e.g., 0, 1, 2 for diploids).
#'
#' @param ploidy Integer specifying the organism ploidy (e.g., 2 for diploid).
#'
#' @param n.clusters.thr Integer specifying the minimum number of genotype clusters
#'   required for a marker to be retained. If fewer clusters are found, missing ones
#'   can be imputed depending on the `type`.
#'   Defaults to `ploidy + 1` if `NULL`.
#'
#' @param type Character string indicating the data source type:
#'   - `"intensities"`: For array-based allele intensities.
#'   - `"counts"`: For sequencing read counts.
#'   Default is `"intensities"`.
#'
#' @param rm_outlier Logical; if `TRUE`, outlier samples within genotype clusters
#'   will be identified and removed prior to center calculation (default: `TRUE`).
#'
#' @param cluster_median Logical; if `TRUE`, cluster centers are calculated using
#'   the median of `theta` values. If `FALSE`, the mean is used (default: `TRUE`).
#'
#' @return A named list with the following elements:
#'   - `rm`: Integer flag: `0` (retained), `1` (no clusters found), or `2` (too few clusters).
#'   - `centers_theta`: A numeric vector of cluster center positions on the theta scale.
#'   - `MarkerName`: Marker identifier.
#'   - `n.clusters`: Number of clusters (including imputed ones if applicable).
#'
#' @import dplyr
#' @import tidyr
#' @importFrom magrittr "%>%"
#'
#' @export
get_centers <- function(ratio_geno,
                        ploidy,
                        n.clusters.thr = NULL,
                        type = c("intensities", "counts"),
                        rm_outlier = TRUE,
                        cluster_median = TRUE){

  if(all(is.na(ratio_geno$theta))) centers <- which(!is.na(ratio_geno$theta)) else {
    if(is.null(n.clusters.thr)) n.clusters.thr <- ploidy + 1

    # Adjust codification
    ad <- ratio_geno %>% filter(!is.na(theta) & !is.na(geno)) %>% group_by(geno) %>% summarise(mean = mean(theta))
    if(ad$mean[1] > ad$mean[nrow(ad)]) {
      genos = data.frame(geno = 0:ploidy, geno.new = ploidy:0)
      ratio_geno$geno <- genos$geno.new[match(ratio_geno$geno,genos$geno)]
    }

    plot_data_split <- split(ratio_geno, ratio_geno$geno)

    if(rm_outlier) plot_data_split <- lapply(plot_data_split, function(x) rm_outlier(x))

    if(cluster_median){
      centers <- lapply(plot_data_split, function(x) apply(x[,3:4], 2, function(y) median(y, na.rm = TRUE)))
    } else centers <- lapply(plot_data_split, function(x) apply(x[,3:4], 2, function(y) mean(y, na.rm = TRUE)))
  }

  if(length(centers) == 0 || length(centers) < n.clusters.thr) {
    return(list(rm= {if(length(centers) == 0) 1 else if(length(centers) < n.clusters.thr) 2},
                centers_theta = centers,
                MarkerName = unique(ratio_geno$MarkerName),
                n.clusters = length(centers)))
  } else {
    centers_df <- data.frame(dose = names(centers), do.call(rbind, centers))

    # If one or more cluster is missing, this code will input using the mean distance between the other clusters
    if(length(centers) < ploidy + 1){
      doses <- 0:ploidy
      new_centers_df <- data.frame(dose = doses)
      new_centers_df$theta <- NA
      new_centers_df$theta[match(names(centers), doses)] <- centers_df$theta
      mis <- which(is.na(new_centers_df$theta))
      wi <- which(!is.na(new_centers_df$theta))

      if(length(type == 2)) select_type <- match.arg(type)
      if(select_type == "counts"){
        loop <- sort(mis, decreasing = T)
        if(any(loop == 1)) loop <- c(1,loop[-which(loop==1)]) # 1 and ploidy + 1 first in the loop
        for(miss_i in loop){
          if(miss_i == 1) { # Avoid numbers < 0 and > 1 as theta
            new_centers_df[1,2] <- 0
          } else if(miss_i == nrow(new_centers_df)) {
            new_centers_df[nrow(new_centers_df),2] <- 1
          } else {
            # find the interval
            before <- after <- miss_i
            before <- before - 1
            after <- after + 1
            while(is.na(new_centers_df$theta[before])) {
              before <- before - 1
            }
            while(is.na(new_centers_df$theta[after])) {
              after <- after + 1
            }
            new_centers_df$theta[miss_i] <- new_centers_df$theta[before] + (miss_i-before)*((new_centers_df$theta[after] - new_centers_df$theta[before])/(after-before))
          }
        }
      } else if(select_type == "intensities"){
        input <- mean(diff(centers_df$theta)/diff(as.numeric(names(centers))))
        for (i in seq_along(mis)) {
          if (mis[i] < wi[1]) {
            new_centers_df$theta[mis[i]] <- new_centers_df$theta[wi[1]] -
              (wi[1] - mis[i]) * input
          }
          if (mis[i] > wi[1]) {
            new_centers_df$theta[mis[i]] <- new_centers_df$theta[wi[length(wi)]] +
              diff(c(wi[length(wi)], mis[i])) * input
          }
        }
      }
      centers_df <- new_centers_df
    }

    centers_df <- cbind(centers_df, cluster = seq_len(nrow(centers_df)))

    return(list(rm = 0,
                centers_theta = centers_df$theta,
                MarkerName = unique(ratio_geno$MarkerName),
                n.clusters = length(centers)))

  }
}

#' Calculate Z-Scores for Allele Intensities or Counts
#'
#' This function computes per-marker Z-scores based on the total signal intensity (R),
#' which typically represents the sum of reference (X) and alternative (Y) allele signals.
#' The Z-score measures how much each sample deviates from the mean intensity of that marker.
#'
#' The function also merges positional metadata from the `geno.pos` input, adding chromosome
#' and physical position for each marker.
#'
#' @param data A data.frame containing signal intensity and ratio values with the following columns:
#' \describe{
#'   \item{MarkerName}{Marker identifiers.}
#'   \item{SampleName}{Sample identifiers.}
#'   \item{X}{Reference allele intensity or count.}
#'   \item{Y}{Alternative allele intensity or count.}
#'   \item{R}{Total signal or depth (i.e., \code{X + Y}).}
#'   \item{ratio}{Allelic ratio, typically \code{Y / (X + Y)}.}
#' }
#'
#' @param geno.pos A data.frame with marker genomic positions, containing the following columns:
#' \describe{
#'   \item{MarkerName}{Marker identifiers.}
#'   \item{Chromosome}{Chromosome identifier where the marker is located.}
#'   \item{Position}{Genomic position (base-pair coordinate) of the marker.}
#' }
#'
#' @return A data.frame containing the following columns:
#' \describe{
#'   \item{MarkerName}{Marker ID.}
#'   \item{Chr}{Chromosome corresponding to the marker.}
#'   \item{Position}{Genomic position (bp).}
#'   \item{SampleName}{Sample ID.}
#'   \item{z}{Z-score computed per marker across all samples.}
#' }
#' Markers with missing chromosome or position information are excluded from the final output.
#'
#' @import dplyr
#' @import tidyr
#' @importFrom magrittr "%>%"
#'
#' @examples
#' data <- data.frame(
#'   MarkerName = rep("m1", 5),
#'   SampleName = paste0("S", 1:5),
#'   X = c(100, 110, 90, 95, 85),
#'   Y = c(200, 190, 210, 205, 215),
#'   R = c(300, 300, 300, 300, 300),
#'   ratio = c(0.67, 0.63, 0.70, 0.68, 0.72)
#' )
#' geno.pos <- data.frame(MarkerName = "m1", Chromosome = "1", Position = 123456)
#' get_zscore(data, geno.pos)
#'
#' @export
get_zscore <- function(data = NULL,
                       geno.pos = NULL){

  stopifnot(all(c("MarkerName", "SampleName", "R") %in% colnames(data)))

  data$R[data$R == 0] <- NA

  zscore <- data %>% group_by(MarkerName) %>%
    mutate(z = (R - mean(R, na.rm = TRUE))/sd(R, na.rm = TRUE)) %>% select(MarkerName, SampleName, z)

  colnames(geno.pos)[1] <- c("MarkerName")

  chr <- geno.pos$Chromosome[match(zscore$MarkerName,geno.pos$MarkerName)]
  pos <- geno.pos$Position[match(zscore$MarkerName,geno.pos$MarkerName)]

  zscore <- cbind(MarkerName=zscore$MarkerName, Chr = chr, Position = pos, zscore[,-1])

  if(length(which(is.na(zscore$Chr)))> 0)
    zscore <- zscore[-which(is.na(zscore$Chr)),]

  return(zscore)
}


#' Identify and Remove Outliers Based on Bonferroni-Holm Adjusted P-values
#'
#' This function detects and removes outlier observations from a vector of `theta` values
#' using externally studentized residuals and the Bonferroni-Holm adjustment for multiple testing.
#' It is typically used during genotype cluster center estimation to clean noisy values.
#'
#' The method fits a constant model (`theta ~ 1`) and computes standardized residuals.
#' Observations with significant deviation are flagged using the Bonferroni-Holm procedure
#' and removed if their adjusted p-value is below the defined `alpha` threshold.
#'
#' This function was originally developed by **Kaio Olympio** and incorporated into the Qploidy workflow.
#'
#' @param data A data.frame containing a `theta` column.
#'   This is usually a subset of the full dataset, representing samples within a single genotype class.
#'
#' @param alpha Significance level for identifying outliers (default is `0.05`).
#'   Observations with adjusted p-values below this threshold will be removed.
#'
#' @param z Logical; if `TRUE`, evaluates outliers in the `z` column and sets outlier `z` values to `NA`,
#'   returning the full data.frame. Otherwise, original behavior is preserved.
#'
#' @return A data.frame containing only the non-outlier observations from the input.
#'   If fewer than two non-NA `theta` values are present or if all values are identical,
#'   the input is returned unmodified.
#'
#' @author Kaio Olympio
#'
#' @importFrom stats lm pnorm sd
#' @import multtest
#'
#' @export
rm_outlier <- function(data, alpha=0.05, z=FALSE){
  # If z=TRUE, operate on z column; else, on theta
  col <- if (z) "z" else "theta"
  vec <- data[[col]]
  rm.na <- which(is.na(vec))
  if(length(rm.na) > 0) vec <- vec[-rm.na]
  if(length(vec) < 2 || length(unique(vec)) == 1) {
    if(z) {
      data$outlier <- FALSE
    }
    return(data)
  } else {
    lm.object <- lm(vec ~ 1)
    resid <- lm.object$residuals
    studresid <- resid/sd(resid, na.rm=TRUE)
    rawp.BHStud = 2 * (1 - pnorm(abs(studresid)))
    test.BHStud <- multtest::mt.rawp2adjp(rawp.BHStud,proc=c("Holm"),alpha = alpha)
    adjp = cbind(test.BHStud[[1]][,1])
    bholm = cbind(test.BHStud[[1]][,2])
    index = test.BHStud$index
    out_flag = ifelse(bholm<alpha, "OUTLIER ", ".")
    BHStud_test = cbind(adjp,bholm,index,out_flag)
    BHStud_test2 = BHStud_test[order(index),]
    names = c("rawp","bholm","index","out_flag")
    colnames(BHStud_test2) <- names
    outliers_BH <- as.numeric(BHStud_test2[which(BHStud_test2[,"out_flag"]!="."),"index"])
    if(z) {
      # Add outlier column: TRUE for outliers, FALSE otherwise
      data$outlier <- FALSE
      if(length(outliers_BH) > 0) {
        data[["z"]][outliers_BH] <- NA
        data$outlier[outliers_BH] <- TRUE
      }
      return(data)
    } else {
      if(length(outliers_BH) >0) new.vec <- data[-outliers_BH,] else new.vec <- data
      return(new.vec)
    }
  }
}