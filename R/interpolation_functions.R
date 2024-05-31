globalVariables(c("theta", "R", "geno", "Var1"))

#' Get centers for standardization using updog estimated bias
#'
#' @param multidog_obj object of class multidog (from updog)
#' @param threshold.n.clusters minimum number of dosage clusters (heterozygous classes) to account with the marker for standardization
#' @param rm.mks vector for logical indicating which markers should be removed, names of the vector are names of the markers
#'
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
    centers_theta <- updog:::xi_fun(p = (0:ploidy[i])/ploidy[i], eps = seq[i], h = bias[i])

    result[[i]] <- list(rm = if(n.clusters.df[i,]$n.clusters >= threshold.n.clusters) 0 else 1,
                        centers_theta = sort(1 - centers_theta),
                        MarkerName = n.clusters.df[i,]$snp,
                        n.clusters = n.clusters.df[i,]$n.clusters)
  }

  names(result) <- snpdf$snp
  return(result)
}

#' Create BAF According to Wang 2007
#'
#' @param theta_subject numeric theta to be standardized
#' @param centers_theta theta centroids defined by clusterization
#' @param ploidy integer defining ploidy
#'
#' @export
get_baf <- function(theta_subject, centers_theta, ploidy){
  baf <- vector()
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
  for(i in 1:nrow(par_all_item[[1]])){
    baf[[i]] <- get_baf(theta_subject = as.numeric(par_all_item[[1]][i,-1]),
                        centers_theta = as.numeric(par_all_item[[2]][i,-1]),
                        ploidy = ploidy)
  }

  return(baf)
}

#' Define cluster centers
#'
#' @param ratio_geno ratio and genotype information for each marker and sample. data.frame with columns: 1) mks; 2) ind; 3) theta; 4) geno.
#' @param ploidy integer defining ploidy
#' @param n.clusters.thr minimum number of clusters required. If clusters < ploidy + 1, the missing clusters will be imputed
#' @param type if array alleles intensities data "intensities", if read counts data "counts"
#' @param rm_outlier if TRUE, it remove the outlier markers before defining the cluster centers
#'
#' @import tidyr
#' @import dplyr
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
      if(type == "counts"){
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
      } else if(type == "intensities"){
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

#' Calculates z scores
#'
#' @param data data.frame with columns: 1) MarkeName: markers IDs; 2) SampleName: Samples IDs; 3) X: reference allele intensities or counts; 4) Y: alternative allele intensities or counts, 5) R: sum of the intensities; and 6) ratio: Y/(X+Y)
#'
#' @param geno.pos data.frame with columns: 1) MarkerName: markers IDs; 2) Chromosome: chromosome where the marker is located; 3) Position: position on the chromosome the marker is located (bp).
#'
#' @param zscore_file_name zscore output filename
#'
#' @import tidyr
#' @import dplyr
#'
#' @export
get_zscore <- function(data = NULL,
                       geno.pos = NULL){

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


##' Identify outliers using Bonferroni-Holm tests for the adjusted p-values and remove them from the input vector
##' @param resid residuals (e.g. lm.object$residuals)
##'
##' @import multtest
##'
##' @export
rm_outlier <- function(plot_data_split_one, alpha=0.05){
  # Produce externally standardized residuals
  theta <- plot_data_split_one$theta
  rm.na <- which(is.na(theta))
  if(length(rm.na) > 0) theta <- theta[-rm.na]
  if(length(theta) < 2 | length(unique(theta)) == 1) return(plot_data_split_one) else {
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
    if(length(outliers_BH) >0) new.theta <- plot_data_split_one[-outliers_BH,] else new.theta <- plot_data_split_one

    return(new.theta)
  }
}

#' Performs standardization and returns BAF and z score
#'
#' @param data data.frame with columns: 1) MarkeName: markers IDs; 2) SampleName: Samples IDs; 3) X: reference allele intensities or counts; 4) Y: alternative allele intensities or counts; 5) R: sum of the intensities; and 6) ratio: Y/(X+Y)
#'
#' @param genos data.frame with genotype information for individuals to be used as reference for standardization. We suggest to select for these individuals that are euploid and all have same ploidy. For array technologies, we suggest using fitpoly for obtaining the dosages, and for sequencing technologies, updog. This file has as columns: 1) MarkerName: markers IDs; 2) SampleName: Samples IDs; 3) geno: dosage.
#'
#' @param geno.pos data.frame with columns: 1) MarkerName: markers IDs; 2) Chromosome: chromosome where the marker is located; 3) Position: position on the chromosome the marker is located (bp).
#'
#' @param threshold.missing.geno fraction of missing genotype information allowed by marker. Markers with higher fraction are discarted
#'
#' @param threshold.geno.prob minimum genotype probability allowed. Genotypes with lower probability will be replaced by NA
#'
#' @param ploidy.standardization ploidy of the reference samples defined in `genos`
#'
#' @param threshold.n.clusters minimum number of dosage clusters (heterozygous classes) to account with the marker for standardization
#'
#' @param n.cores number of cores to be used in parallelized processes
#'
#' @param type method to determine the clusters centers for each marker. It can be "intensities" if array data, "counts" if sequencing data or "updog" if multidog object is provided in multidog_obj argument
#' @param multidog_obj object of class multidog from updog package analysis
#' @param out_filename output file name
#'
#' @param verbose If TRUE display informative messages
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
  prob.rm <- round(as.numeric(prob.rm["TRUE"]/sum(prob.rm)*100),2)
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

  # bugfix
  legths <- sapply(bafs_lt, length)
  keep.mks.n <- which(legths == names(which.max(table(legths))))

  bafs_lt <- bafs_lt[keep.mks.n]
  bafs_m <- do.call(rbind, bafs_lt)
  rownames(bafs_m) <- theta_filt$MarkerName[keep.mks.n]
  colnames(bafs_m) <- colnames(theta_filt)[-1]
  bafs_df <- as.data.frame(bafs_m)
  bafs_df <- cbind(mks = theta_filt$MarkerName[keep.mks.n], bafs_df)

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
#' @param ...
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
#'
#' @importFrom ggpubr ggarrange annotate_figure text_grob
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#'
#' @return printed information about Qploidy standardization process
#'
#' @method plot qploidy_standardization
#'
#' @export
plot.qploidy_standardization <- function(x,
                                         sample = NULL,
                                         chr = NULL,
                                         type = c("all", "het", "BAF","zscore","BAF_hist", "BAF_hist_overall", "Ratio_hist_overall"),
                                         area_single = 0.75,
                                         ploidy = 4,
                                         dot.size = 1,
                                         add_estimated_peaks = FALSE,
                                         add_expected_peaks = FALSE,
                                         centromeres = NULL,
                                         add_centromeres = FALSE,
                                         colors = FALSE,
                                         window_size = 2000000,
                                         het_interval = 0.1,
                                         rm_homozygous = FALSE){

  if(!inherits(x, "qploidy_standardization")) stop("Object is not of class qploidy_standardization")

  if(is.null(sample)) stop("Define sample ID")

  if(is.numeric(chr)) chr <- sort(unique(x$data$Chr))[chr] else if(is.null(chr)) chr <- sort(unique(x$data$Chr))

  data_sample <- x$data %>% select(MarkerName, SampleName, Chr, Position, baf, z, ratio) %>% filter(SampleName == sample & Chr %in% chr)

  if(nrow(data_sample) == 0) stop("Sample or chromosome not found.")

  baf_sample <- data_sample %>% pivot_wider(names_from = SampleName, values_from = baf)
  zscore_sample <- data_sample %>% pivot_wider(names_from = SampleName, values_from = z)

  colnames(baf_sample)[ncol(baf_sample)] <-  "sample"

  baf_point <- baf_hist <- p_z <- raw_ratio <- het_rate <- baf_hist_overall <- NULL

  if(any(type == "all" | type == "BAF")){
    baf_point <- plot_baf(baf_sample,
                          area_single,
                          ploidy,
                          dot.size = 1,
                          add_estimated_peaks,
                          add_expected_peaks,
                          centromeres,
                          add_centromeres,
                          colors)
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
                              rm_homozygous = rm_homozygous)
  }

  if(any(type == "all" | type == "BAF_hist_overall")){
    baf_hist_overall <- plot_baf_hist(data_sample = baf_sample,
                              area_single,
                              ploidy,
                              colors,
                              add_estimated_peaks,
                              add_expected_peaks,
                              BAF_hist_overall = TRUE,
                              rm_homozygous = rm_homozygous)
  }

  if(any(type == "all" | type == "Ratio_hist_overall")){
    baf_hist_overall <- plot_baf_hist(data_sample = baf_sample,
                                      area_single,
                                      ploidy,
                                      colors,
                                      add_estimated_peaks,
                                      add_expected_peaks,
                                      BAF_hist_overall = TRUE,
                                      ratio = TRUE,
                                      rm_homozygous = rm_homozygous)
  }

  if(any(type == "zscore")){
    colnames(zscore_sample)[ncol(zscore_sample)] <- "z"
    p_z <- zscore_sample  %>%
      ggplot(aes(x = Position , y = z)) +
      facet_grid(.~Chr, scales = "free") +
      geom_smooth(method = "gam") +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      geom_hline(yintercept=median(zscore_sample$z), linetype="dashed")

  }

  if(any(type == "ratio")){
    raw_ratio <- data_sample %>%  ggplot(aes(x = Position, y = ratio)) + geom_point(alpha =0.7, size=dot.size) +
      facet_grid(.~Chr, scales = "free") + theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
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

    het_rate <- data_sample %>%  group_by(Chr, window) %>% summarise(prop_het = sum(het_rate, na.rm = TRUE)/length(which(het_rate == 1 | het_rate == 0))) %>%
      ggplot(aes(x = window, y = prop_het, color = prop_het)) + geom_line() +
      scale_color_gradient(low = "black", high = "red") +
      facet_grid(.~Chr, scales = "free") + theme_bw() + ylab("proportion of heterozygous loci") + xlab("Position")+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
  }

  p_all <- list(het_rate, raw_ratio, baf_point, baf_hist,baf_hist_overall, p_z)

  rm <- which(sapply(p_all, is.null))
  if(length(rm) != 0) p_all <- p_all[-rm]

  p_result <- ggarrange(plotlist = p_all, ncol = 1)

  p_result <- annotate_figure(p_result, top = text_grob(sample, face = "bold", size = 14))

  return(p_result)
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
