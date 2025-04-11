globalVariables(c("theta", "R", "geno", "Var1", "array.id",
                  "ids", "sd", "median", "z", "lm", "sd",
                  "snp", "centers_theta", "pnorm"))

#' Get centers for standardization using updog estimated bias
#'
#' @param multidog_obj object of class multidog (from updog)
#' @param threshold.n.clusters minimum number of dosage clusters (heterozygous classes) to account with the marker for standardization
#' @param rm.mks vector for logical indicating which markers should be removed, names of the vector are names of the markers
#'
#'@export
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
  for(i in 1:length(bias)){
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
#' @export
get_baf_par <- function(par_all_item, ploidy=2){
  baf <- list()
  i <- 1
  for(i in 1:nrow(par_all_item[[1]])){
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

  if(length(centers) == 0 | length(centers) < n.clusters.thr) {
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
        for(i in 1:length(mis)){
          if(mis[i] < wi[1]) new_centers_df$theta[mis[i]] <- new_centers_df$theta[wi[1]] - (wi[1]-mis[i])*input
          if(mis[i] > wi[1]) new_centers_df$theta[mis[i]] <- new_centers_df$theta[wi[length(wi)]] + diff(c(wi[length(wi)],mis[i]))*input
        }
      }
      centers_df <- new_centers_df
    }

    centers_df <- cbind(centers_df, cluster = 1:nrow(centers_df))

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
#' @return A data.frame containing only the non-outlier observations from the input.
#'   If fewer than two non-NA `theta` values are present or if all values are identical,
#'   the input is returned unmodified.
#'
#' @importFrom stats lm pnorm sd
#' @import multtest
#'
#' @export
rm_outlier <- function(data, alpha=0.05){
  # Produce externally standardized residuals
  theta <- data$theta
  rm.na <- which(is.na(theta))
  if(length(rm.na) > 0) theta <- theta[-rm.na]
  if(length(theta) < 2 | length(unique(theta)) == 1) return(data) else {
    lm.object <- lm(theta ~ 1)
    resid <- lm.object$residuals
    studresid <- resid/sd(resid, na.rm=TRUE)
    # Calculate adjusted p-values
    rawp.BHStud = 2 * (1 - pnorm(abs(studresid)))
    #Produce a Bonferroni-Holm tests for the adjusted p-values
    #The output is a list
    test.BHStud<-mt.rawp2adjp(rawp.BHStud,proc=c("Holm"),alpha = alpha)
    #Create vectors/matrices out of the list of the BH tests
    adjp = cbind(test.BHStud[[1]][,1])
    bholm = cbind(test.BHStud[[1]][,2])
    index = test.BHStud$index
    # Condition to flag outliers according to the BH test
    out_flag = ifelse(bholm<alpha, "OUTLIER ", ".")
    #Create a matrix with all the output of the BH test
    BHStud_test = cbind(adjp,bholm,index,out_flag)
    #Order the file by index
    BHStud_test2 = BHStud_test[order(index),]
    #Label colums
    names = c("rawp","bholm","index","out_flag")
    colnames(BHStud_test2) <- names
    #Create a final file, with the data and the test and the labels for the outliers

    # Take a look at the outliers
    outliers_BH <- as.numeric(BHStud_test2[which(BHStud_test2[,"out_flag"]!="."),"index"])
    if(length(outliers_BH) >0) new.theta <- data[-outliers_BH,] else new.theta <- data

    return(new.theta)
  }
}

#' Standardize Allelic Ratio Data and Compute BAF and Z-Scores
#'
#' This function performs signal standardization of genotype data by aligning `theta` values
#' (allelic ratios or normalized intensities) to expected genotype clusters. It outputs
#' standardized BAF (B-allele frequency) and Z-scores per sample and marker.
#'
#' Reference genotypes are used to estimate cluster centers either from dosage data
#' (e.g., via `fitpoly` or `updog`) or using an `updog` `multidog` object directly.
#' This function supports both array-based (intensity) and sequencing-based (count) data.
#'
#' It applies marker and genotype-level quality filters, uses parallel computing
#' to estimate BAF, and generates a final annotated output suitable for CNV or
#' dosage variation analyses.
#'
#' @param data A `data.frame` containing the full dataset with the following columns:
#' \describe{
#'   \item{MarkerName}{Marker identifiers.}
#'   \item{SampleName}{Sample identifiers.}
#'   \item{X}{Reference allele intensity or count.}
#'   \item{Y}{Alternative allele intensity or count.}
#'   \item{R}{Total signal intensity or read depth (X + Y).}
#'   \item{ratio}{Allelic ratio, typically Y / (X + Y).}
#' }
#'
#' @param genos A `data.frame` containing genotype dosage information for the reference panel.
#' This should include samples of known ploidy and ideally euploid individuals.
#' Required columns:
#' \describe{
#'   \item{MarkerName}{Marker identifiers.}
#'   \item{SampleName}{Sample identifiers.}
#'   \item{geno}{Estimated dosage (0, 1, 2, ...).}
#'   \item{prob}{Genotype call probability (used for filtering low-confidence genotypes).}
#' }
#'
#' @param geno.pos A `data.frame` with marker position metadata. Required columns:
#' \describe{
#'   \item{MarkerName}{Marker identifiers.}
#'   \item{Chromosome}{Chromosome names.}
#'   \item{Position}{Base-pair positions on the genome.}
#' }
#'
#' @param threshold.missing.geno Numeric (0–1). Maximum fraction of missing genotype data allowed per marker.
#' Markers with a higher fraction will be removed.
#'
#' @param threshold.geno.prob Numeric (0–1). Minimum genotype call probability threshold.
#' Genotypes with lower probability will be treated as missing.
#'
#' @param ploidy.standardization Integer. The ploidy level of the reference panel used for standardization.
#'
#' @param threshold.n.clusters Integer. Minimum number of expected dosage clusters per marker.
#' For diploid data, this is typically 3 (corresponding to genotypes 0, 1, and 2).
#'
#' @param n.cores Integer. Number of cores to use in parallel computations (e.g., for cluster center estimation and BAF generation).
#'
#' @param type Character. Type of data used for clustering:
#' \describe{
#'   \item{"intensities"}{For array-based allele intensity data.}
#'   \item{"counts"}{For sequencing data.}
#'   \item{"updog"}{Automatically set when `multidog_obj` is provided.}
#' }
#'
#' @param multidog_obj Optional. An object of class `multidog` from the `updog` package, containing model fits and estimated biases.
#' If provided, this will override the `type` parameter and use `updog`'s expected cluster positions.
#'
#' @param out_filename Optional. Path to save the final standardized dataset to disk as a CSV file (suitable for Qploidy).
#'
#' @param parallel.type Character. Parallel backend to use (`"FORK"` or `"PSOCK"`). `"FORK"` is faster but only works on Unix-like systems.
#'
#' @param rm_outlier Logical. If `TRUE`, uses Bonferroni-Holm corrected residuals to remove outliers before estimating cluster centers.
#'
#' @param cluster_median Logical. If `TRUE`, uses the median of theta values to estimate cluster centers. If `FALSE`, uses the mean.
#'
#' @param verbose Logical. If `TRUE`, prints progress and filtering information to the console.
#'
#' @return An object of class `"qploidy_standardization"` (list) with the following components:
#' \describe{
#'   \item{info}{Named vector of standardization parameters.}
#'   \item{filters}{Named vector summarizing how many markers were removed at each filtering step.}
#'   \item{data}{A data.frame containing merged BAF, Z-score, and genotype information by marker and sample.}
#' }
#'
#'
#' @import dplyr
#' @import tidyr
#' @import vroom
#' @import parallel
#' @import lme4
#' @import emmeans
#'
#' @export
standardize <- function(data = NULL,
                        genos = NULL,
                        geno.pos = NULL,
                        threshold.missing.geno=0.90,
                        threshold.geno.prob=0.8,
                        ploidy.standardization = NULL,
                        threshold.n.clusters = NULL,
                        n.cores =1,
                        out_filename = NULL,
                        type = "intensities",
                        multidog_obj = NULL,
                        parallel.type = "PSOCK",
                        verbose = TRUE,
                        rm_outlier = TRUE,
                        cluster_median = TRUE){

  if(is.null(data) | is.null(genos) | is.null(geno.pos) | is.null(ploidy.standardization) | is.null(threshold.n.clusters)) stop("Not all required inputs were defined.")

  if(!all(colnames(data) %in% c("MarkerName", "SampleName", "X", "Y", "R", "ratio"))) stop("Column names of the provided data object does not match the required.")
  if(!all(colnames(genos) %in% c("MarkerName", "SampleName", "geno", "prob"))) stop("Column names of the provided genos object does not match the required.")
  if(!all(colnames(geno.pos) %in% c("MarkerName", "SampleName", "Chromosome", "Position"))) stop("Column names of the provided geno.pos object does not match the required.")

  dose <- max(genos$geno, na.rm = T)
  if(dose != ploidy.standardization) stop("Ploidy of the provided reference samples do not match with the one defined in the ploidy.standardization parameter.")

  if(verbose) cat("Generating standardize BAFs...\n")
  if(is.null(threshold.n.clusters)) threshold.n.clusters <- ploidy.standardization + 1

  ## Filter by prob
  idx <- genos$prob < threshold.geno.prob
  prob.rm <- table(idx)
  if(!is.na(prob.rm["TRUE"])) prob.rm <- round(as.numeric(prob.rm["TRUE"]/sum(prob.rm)*100),2) else prob.rm <- 0
  idx <- which(idx)
  if(length(idx) > 0) genos$geno[idx] <- NA

  if(verbose) print(paste0("Percentage of genotypes turned into missing data because of low genotype probability:", prob.rm))

  ## Filter by missing data
  n.na <- genos %>% group_by(MarkerName) %>% summarize(n.na = (sum(is.na(geno))/length(geno)))
  rm.mks <- n.na$MarkerName[which(n.na$n.na > threshold.missing.geno)]
  mis.rm <- length(rm.mks)

  if(verbose) print(paste0("Markers remove because of excess of missing data:", mis.rm))

  if(length(rm.mks) > 0){
    genos_filt <- genos[which(!(genos$MarkerName %in% rm.mks)),]
    data_filt <- data[which(!(data$MarkerName %in% rm.mks)),]
  } else {
    genos_filt <- genos
    data_filt <- data
  }

  # Organize data
  data_standardization <- inner_join(data_filt[,-c(3,4,5)], genos_filt[,-4], by = c("MarkerName", "SampleName"))

  if(dim(data_standardization)[1] == 0) stop("Individuals in `data` and `genos` don't have same ID.")

  colnames(data_standardization)[3] <- "theta"

  rm.na <- which(is.na(data_standardization$geno))
  if(length(rm.na) > 0){
    data_standardization$theta[rm.na] <- NA
  }

  if(is.null(multidog_obj)){
    lst_standardization <- split(data_standardization, data_standardization$MarkerName)

    if(verbose) cat("Going to parallel mode...\n")
    clust <- makeCluster(n.cores, type = parallel.type)
    clusterExport(clust, c("get_centers"))
    clusters <- parLapply(clust, lst_standardization, get_centers,
                          ploidy= ploidy.standardization,
                          n.clusters.thr = threshold.n.clusters,
                          type = type,
                          rm_outlier = rm_outlier,
                          cluster_median = cluster_median)

    stopCluster(clust)

    if(verbose) cat("Back to single core usage\n")

    gc(verbose = FALSE)

  } else { # centers defined using updog bias
    clusters <- updog_centers(multidog_obj, threshold.n.clusters = threshold.n.clusters, rm.mks = rm.mks)
  }

  # Filter by number of clusters
  rm.mks <- sapply(clusters, function(x) x$rm != 0)
  clusters.rm <- sum(rm.mks)

  if(verbose) print(paste0("Markers remove because of smaller number of clusters than set threshold:",clusters.rm))

  if(length(which(rm.mks)) > 0)  clusters_filt <- clusters[-which(rm.mks)] else clusters_filt <- clusters

  if(length(clusters_filt) == 0) stop("All markers were filtered, adapt filtering parameters.")

  keep.mks <- sapply(clusters_filt, function(x) x$MarkerName)

  # Getting BAF for complete dataset
  theta_filt <- pivot_wider(data[which(data$MarkerName %in% keep.mks),-c(3:5)], names_from = "SampleName", values_from = "ratio")
  centers_filt <- t(sapply(clusters_filt, function(x) x$centers_theta))
  centers_filt <- data.frame(mks=keep.mks, centers_filt)
  centers_filt <- centers_filt[match(theta_filt$MarkerName, centers_filt$mks),]

  par <- rep(1:n.cores, each=round((nrow(theta_filt)/n.cores)+1,0))[1:nrow(theta_filt)]

  par_theta <- split.data.frame(theta_filt, par)
  par_clusters_filt <- split(centers_filt, par)

  if(length(par_theta) < n.cores) n.cores <- length(par_theta)
  par_all <- list()
  for(i in 1:n.cores){
    par_all[[i]] <- list()
    par_all[[i]][[1]] <- as.data.frame(par_theta[[i]])
    par_all[[i]][[2]] <- as.data.frame(par_clusters_filt[[i]])
  }

  # Get BAF
  if(verbose) cat("Going to parallel mode...\n")

  clust <- makeCluster(n.cores, type = parallel.type)
  clusterExport(clust, c("get_baf_par"))
  bafs <- parLapply(clust, par_all, get_baf_par, ploidy = ploidy.standardization)

  stopCluster(clust)

  if(verbose) cat("Back to single core usage\n")

  gc()

  bafs_lt <- unlist(bafs, recursive = F)
  bafs_m <- do.call(rbind, bafs_lt)
  rownames(bafs_m) <- theta_filt$MarkerName
  colnames(bafs_m) <- colnames(theta_filt)[-1]
  bafs_df <- as.data.frame(bafs_m)
  bafs_df <- cbind(mks = theta_filt$MarkerName, bafs_df)

  # Add chr and pos info
  chr <- geno.pos$Chromosome[match(bafs_df$mks,geno.pos$MarkerName)]
  pos <- geno.pos$Position[match(bafs_df$mks,geno.pos$MarkerName)]

  bafs_join <- cbind(MarkerName=bafs_df$mks, Chr = chr, Position = pos, bafs_df[,-1])

  no.geno.info <- length(which(is.na(bafs_join$Chr)))
  if(no.geno.info > 0)
    bafs_join <- bafs_join[-which(is.na(bafs_join$Chr)),]

  baf_melt <- pivot_longer(bafs_join, cols = 4:ncol(bafs_join), names_to = "SampleName", values_to = "baf")
  if(verbose) cat("BAFs ready!\n")

  # Z score
  if(verbose) cat("Generating z scores...\n")
  zscore <- get_zscore(data, geno.pos)
  if(verbose) cat("Z scores ready!\n")

  if(verbose) cat("Merging results into qploidy_standardization object...\n")
  qploidy_data <- full_join(data, data_standardization[,-3], c("MarkerName", "SampleName"))
  qploidy_data <- full_join(qploidy_data,baf_melt, c("MarkerName", "SampleName"))
  qploidy_data <- full_join(qploidy_data[,-c(8,9)], zscore, c("MarkerName", "SampleName"))

  result <- structure(list(info = c(threshold.missing.geno = threshold.missing.geno,
                                    threshold.geno.prob = threshold.geno.prob,
                                    ploidy.standardization = ploidy.standardization,
                                    threshold.n.clusters = threshold.n.clusters,
                                    out_filename = out_filename,
                                    type = if(!is.null(multidog_obj)) "updog" else type),
                           filters = c(n.markers.start = length(unique(data$MarkerName)),
                                       geno.prob.rm = prob.rm,
                                       miss.rm = mis.rm,
                                       clusters.rm= clusters.rm,
                                       no.geno.info.rm = no.geno.info,
                                       n.markers.end = length(unique(baf_melt$MarkerName))),
                           data = qploidy_data), class = "qploidy_standardization")

  if(!is.null(out_filename)) {
    if(verbose) cat(paste0("Writting Qploidy app input file:", out_filename))
    info = data.frame(t(result$info))
    filters = data.frame(t(result$filters))
    vroom_write(info, file = out_filename, col_names = T)
    vroom_write(filters, file = out_filename, append = TRUE, col_names = T)
    vroom_write(result$data, file = out_filename, append = TRUE, col_names = T)
  }

  result <- structure(list(info = c(threshold.missing.geno = threshold.missing.geno,
                                    threshold.geno.prob = threshold.geno.prob,
                                    ploidy.standardization = ploidy.standardization,
                                    threshold.n.clusters = threshold.n.clusters,
                                    out_filename = out_filename,
                                    type = if(!is.null(multidog_obj)) "updog" else type),
                           filters = c(n.markers.start = length(unique(data$MarkerName)),
                                       geno.prob.rm = prob.rm,
                                       miss.rm = mis.rm,
                                       clusters.rm= clusters.rm,
                                       no.geno.info.rm = no.geno.info,
                                       n.markers.end = length(unique(baf_melt$MarkerName))),
                           data = qploidy_data), class = "qploidy_standardization")
  if(verbose) cat("Done!\n")
  return(result)
}


#' Print method for object of class 'qploidy_standardization'
#'
#' @param x object of class 'qploidy_standardization'
#' @param ... print parameters
#'
#' @return printed information about Qploidy standardization process
#'
#' @method print qploidy_standardization
#'
#' @export
print.qploidy_standardization <- function(x, ...){

  info <- data.frame(c1 = c("standardization type:", "Ploidy:",
                            "Minimum number of heterozygous classes (clusters) present:",
                            "Maximum number of missing genotype by marker:",
                            "Minimum genotype probability:"),
                     c2 = c(x$info["type"], x$info["ploidy.standardization"],
                            x$info["threshold.n.clusters"],
                            1 - as.numeric(x$info["threshold.missing.geno"]),
                            x$info["threshold.geno.prob"]))

  format.df <- data.frame(c1 = c("Number of markers at raw data:",
                                 "Percentage of filtered genotypes by probability threshold:",
                                 "Number of markers filtered by missing data:",
                                 "Number of markers filtered for not having the minimum number of clusters:",
                                 "Number of markers filtered for not having genomic information:",
                                 "Number of markers with estimated BAF:"),
                          c2 = c(x$filters["n.markers.start"],
                                 "-",
                                 x$filters["miss.rm"],
                                 x$filters["clusters.rm"],
                                 x$filters["no.geno.info.rm"],
                                 x$filters["n.markers.end"]),
                          c3 = c("(100%)",
                                 paste0("(",x$filters["geno.prob.rm"], " %)"),
                                 paste0("(",round(x$filters["miss.rm"]/x$filters["n.markers.start"]*100,2)," %)"),
                                 paste0("(", round(x$filters["clusters.rm"]/x$filters["n.markers.start"]*100,2)," %)"),
                                 paste0("(", round(x$filters["no.geno.info.rm"]/x$filters["n.markers.start"]*100,2)," %)"),
                                 paste0("(", round(x$filters["n.markers.end"]/x$filters["n.markers.start"]*100,2)," %)")))

  colnames(info) <- rownames(info) <- colnames(format.df) <- rownames(format.df) <- NULL

  cat("This is on object of class 'ploidy_standardization'\n")
  cat("--------------------------------------------------------------------\n")
  cat("Parameters\n")
  print(format(info, justify="left", digit = 2))
  cat("--------------------------------------------------------------------\n")
  cat("Filters\n")
  print(format(format.df, justify="left", digit = 2))
}

#' Read file generated by standardize function
#'
#' @param qploidy_standardization_file path to file generated by standardize function
#'
#' @import vroom
#'
#' @export
read_qploidy_standardization <- function(qploidy_standardization_file){
  info <- vroom(qploidy_standardization_file,n_max = 1)
  info_v <- as.character(info)
  names(info_v) <- colnames(info)
  filters <- vroom(qploidy_standardization_file, skip = 2, n_max=1)
  filters_v <- as.numeric(filters)
  names(filters_v) <- colnames(filters)
  data <- vroom(qploidy_standardization_file, skip = 4)

  return(structure(list(info = info_v,
                        filters = filters_v,
                        data = data), class = "qploidy_standardization"))
}
