globalVariables(c("theta", "R", "geno", "Var1"))

#' Get centers for standardization using updog estimated bias
#'
#' @param multidog.obj object of class multidog (from updog)
#' @param threshold.n.clusters minimum number of dosage clusters (heterozygous classes) to account with the marker for standardization
#' @param rm.mks vector for logical indicating which markers should be removed, names of the vector are names of the markers
#'
updog_centers <- function(multidog.obj, threshold.n.clusters=2, rm.mks){

  n.clusters.df <- multidog.obj$inddf %>%
    filter(!(snp %in% rm.mks)) %>%
    group_by(snp) %>%
    summarize(n.clusters = length(table(geno)))

  if(length(rm.mks) >0) snpdf <- multidog.obj$snpdf[-which(multidog.obj$snpdf$snp %in% rm.mks),] else snpdf <- multidog.obj$snpdf
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
                        rm_outlier = TRUE){

  if(is.null(n.clusters.thr)) n.clusters.thr <- ploidy + 1

  # Adjust codification
  ad <- ratio_geno %>% filter(!is.na(theta) & !is.na(geno)) %>% group_by(geno) %>% summarise(mean = mean(theta))
  if(ad$mean[1] > ad$mean[nrow(ad)]) {
    genos = data.frame(geno = 0:ploidy, geno.new = ploidy:0)
    ratio_geno$geno <- genos$geno.new[match(ratio_geno$geno,genos$geno)]
  }

  plot_data_split <- split(ratio_geno, ratio_geno$geno)

  if(rm_outlier) plot_data_split <- lapply(plot_data_split, function(x) rm_outlier(x))

  centers <- lapply(plot_data_split, function(x) apply(x[,3:4], 2, function(y) median(y, na.rm = TRUE)))

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
                        multidog.obj = NULL,
                        parallel.type = "PSOCK",
                        verbose = TRUE){

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

  if(is.null(multidog.obj)){
    lst_standardization <- split(data_standardization, data_standardization$MarkerName)

    if(verbose) cat("Going to parallel mode...\n")
    clust <- makeCluster(n.cores, type = parallel.type)
    clusterExport(clust, c("get_centers"))
    clusters <- parLapply(clust, lst_standardization, get_centers,
                          ploidy= ploidy.standardization,
                          n.clusters.thr = threshold.n.clusters,
                          type = type)

    stopCluster(clust)

    if(verbose) cat("Back to single core usage\n")

    gc(verbose = FALSE)

  } else { # centers defined using updog bias
    clusters <- updog_centers(multidog.obj, threshold.n.clusters = threshold.n.clusters, rm.mks = rm.mks)
  }

  # Filter by number of clusters
  rm.mks <- sapply(clusters, function(x) x$rm != 0)
  clusters.rm <- sum(rm.mks)

  if(verbose) print(paste0("Markers remove because of smaller number of clusters than set threshold:",clusters.rm))

  if(length(which(rm.mks)) > 0)  clusters_filt <- clusters[-which(rm.mks)] else clusters_filt <- clusters

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

  if(!is.null(out_filename)) {
    if(verbose) cat(paste0("Writting Qploidy app input file:", out_filename))
    vroom_write(qploidy_data, file = out_filename)
  }

  result <- structure(list(info = c(threshold.missing.geno = threshold.missing.geno,
                                    threshold.geno.prob = threshold.geno.prob,
                                    ploidy.standardization = ploidy.standardization,
                                    threshold.n.clusters = threshold.n.clusters,
                                    out_filename = out_filename,
                                    type = if(!is.null(multidog.obj)) "updog" else type),
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
  cat("This is on object of class 'ploidy_standardization'\n")
  cat("--------------------------------------------------------------------\n")
  cat("The following parameters were used to generate it:\n")
  cat("standardization type:", x$info["type"], "\n")
  cat("Ploidy:", x$info["ploidy.standardization"], "\n")
  cat("Minimum number of heterozygous classes (clusters) present:", x$info["threshold.n.clusters"], "\n")
  cat("Maximum number of missing genotype by marker:",  1 - as.numeric(x$info["threshold.missing.geno"]), "\n")
  cat("Minimum genotype probability:", x$info["threshold.geno.prob"], "\n")
  cat("--------------------------------------------------------------------\n")
  cat("Number of markers after filters applied:\n")
  cat("Number of markers at raw data:", x$filters["n.markers.start"], " (100%)\n")
  cat("Percentage of filtered genotypes by probability threshold:", x$filters["geno.prob.rm"], " % \n")
  cat("Number of markers filtered by missing data:", x$filters["miss.rm"], " (",
      round(x$filters["miss.rm"]/x$filters["n.markers.start"]*100,2),"%)", "\n")
  cat("Number of markers filtered for not having the minimum number of clusters:", x$filters["clusters.rm"], " (",
      round(x$filters["clusters.rm"]/x$filters["n.markers.start"]*100,2),"%)", "\n")
  cat("Number of markers filtered for not having genomic information:", x$filters["no.geno.info.rm"], " (",
      round(x$filters["no.geno.info.rm"]/x$filters["n.markers.start"]*100,2),"%)", "\n")
  cat("Number of markers with estimated BAF:", x$filters["n.markers.end"], " (",
      round(x$filters["n.markers.end"]/x$filters["n.markers.start"]*100,2),"%)", "\n")
}


#' PLot method for object of class 'qploidy_standardization'
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
                                         type = c("all", "het","ratio","BAF","zscore","BAF_hist", "BAF_hist_overall"),
                                         area_single = 0.75,
                                         ploidy = 4,
                                         dot.size = 1,
                                         add_estimated_peaks = FALSE,
                                         add_expected_peaks = FALSE,
                                         centromeres = NULL,
                                         add_centromeres = FALSE,
                                         colors = FALSE,
                                         window_size = 2000000){

  if(!inherits(x, "qploidy_standardization")) stop("Object is not of class qploidy_standardization")

  if(is.null(sample)) stop("Define sample ID")

  if(is.numeric(chr)) chr <- sort(unique(x$data$Chr))[chr] else if(is.null(chr)) chr <- sort(unique(x$data$Chr))

  data_sample <- x$data %>% select(MarkerName, SampleName, Chr, Position, baf, z, ratio) %>% filter(SampleName == sample & Chr %in% chr)
  baf_sample <- data_sample %>% pivot_wider(names_from = SampleName, values_from = baf)
  zscore_sample <- data_sample %>% pivot_wider(names_from = SampleName, values_from = z)

  colnames(baf_sample)[ncol(baf_sample)] <-  "sample"

  baf_point <- baf_hist <- p_z <- raw_ratio <- het_rate <- NULL

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

  if(any(type == "BAF_hist")){
    baf_hist <- plot_baf_hist(baf_sample,
                              area_single,
                              ploidy,
                              colors,
                              add_estimated_peaks,
                              add_expected_peaks,
                              BAF_hist_overall = FALSE)

  } else if(any(type == "all" | type == "BAF_hist_overall")){
    baf_hist <- plot_baf_hist(baf_sample,
                              area_single,
                              ploidy,
                              colors,
                              add_estimated_peaks,
                              add_expected_peaks,
                              BAF_hist_overall = TRUE)
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
    data_sample$het_rate[data_sample$ratio <= 0.01 | data_sample$ratio >= 0.99] <- 0
    data_sample$het_rate[data_sample$ratio > 0.01 & data_sample$ratio < 0.99] <- 1

    data_sample_lst <- split.data.frame(data_sample, data_sample$Chr)

    for(j in 1:length(data_sample_lst)){
      intervals_start <- c(seq(1,max(data_sample_lst[[j]]$Position), window_size))
      intervals_end <- c(seq(1,max(data_sample_lst[[j]]$Position), window_size) - 1, max(data_sample_lst[[j]]$Position))
      intervals_end <- intervals_end[-1]
      data_sample_lst[[j]]$window <- NA
      for(i in 1:length(intervals_start)){
        data_sample_lst[[j]]$window[which(data_sample_lst[[j]]$Position >= intervals_start[i] & data_sample_lst[[j]]$Position <= intervals_end[i])] <- intervals_start[i]
      }
    }

    data_sample <- do.call(rbind, data_sample_lst)

    het_rate <- data_sample %>%  group_by(Chr, window) %>% summarise(n_het = sum(het_rate, na.rm = TRUE)) %>%
      ggplot(aes(x = window, y = n_het, color = n_het)) + geom_line() +
      scale_color_gradient(low = "black", high = "red") +
      facet_grid(.~Chr, scales = "free") + theme_bw() + ylab("# heterozygous locus") + xlab("Position")+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
  }

  p_all <- list(het_rate, raw_ratio, baf_point, baf_hist, p_z)

  rm <- which(sapply(p_all, is.null))
  if(length(rm) != 0) p_all <- p_all[-rm]

  p_result <- ggarrange(plotlist = p_all, ncol = 1)

  annotate_figure(p_result, top = text_grob(sample, face = "bold", size = 14))

  return(p_result)
}
