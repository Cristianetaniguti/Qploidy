globalVariables(c("theta", "R", "geno", "Var1"))

# R_subject <- R_tetra[1,]
# theta_subject <- theta[1,]
# centers_theta <- clusters_filt[[1]]$centers_theta
# clusters_filt[[1]]$plot
# ploidy <- 4

#' Create BAF According to Wang 2007
#'
#' @param theta_subject numeric theta to be interpolated
#' @param centers_theta theta centroids defined by clusterization
#' @param ploidy integer defining ploidy
#'
#' @export
get_baf <- function(theta_subject, centers_theta, ploidy){
  pos <- vector()
  centers_theta <- sort(centers_theta, decreasing = F)
  idx <- which(theta_subject <= centers_theta[1])
  if(length(idx) > 0) pos[idx] <- 1
  for(i in 2:(ploidy+1)){
    idx <- which(theta_subject > centers_theta[i-1] & theta_subject <= centers_theta[i])
    if(length(idx) > 0) pos[idx] <- i
  }
  idx <- which(theta_subject >= centers_theta[ploidy+1])
  if(length(idx) > 0) pos[idx] <- i+1

  ploidy_freq <- seq(0,1,1/(ploidy))
  #ploidy_freq_multi <- c(ploidy_freq[2],ploidy_freq[2:length(ploidy_freq)])
  #ploidy_freq <- c(ploidy_freq[1],rep(ploidy_freq[-c(1,length(ploidy_freq))], each = 2), ploidy_freq[length(ploidy_freq)])
  ploidy_freq_multi <- 1/ploidy
  baf <- rep(NA, length(pos))
  for(i in 1:(ploidy + 2)){
    idx <- which(pos == i)
    if(i == 1 & length(idx) > 0) {
      baf[idx] <- 0
    } else if(i != 1 & i != ploidy + 2 & length(idx) > 0){
      D2 <- centers_theta[i] - centers_theta[i-1]
      D1 <- as.numeric(theta_subject[idx]) - centers_theta[i-1]
      baf[idx] <- ploidy_freq[i-1] + (D1/D2)*ploidy_freq_multi
    } else if(i == ploidy + 2 & length(idx) > 0){
      baf[idx] <- 1
    }
  }
  baf <- unlist(baf)
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
    baf[[i]] <- get_baf(theta_subject = par_all_item[[1]][i,],
                        centers_theta = par_all_item[[2]][[i]]$centers_theta,
                        ploidy = ploidy)
  }
  return(baf)
}


#' Clusterization using fitpoly
#'
#' @param scores_temp fitpoly scores output file
#' @param ploidy integer defining ploidy
#' @param plot logical to create or not the interpolation plot
#' @param n.clusters.thr minimum number of clusters required. If clusters < ploidy + 1, the missing clusters will be imputed.
#'
#' @export
par_fitpoly_interpolation <- function(scores_temp, ploidy, n.clusters.thr = 0){

  if(is.null(n.clusters.thr)) n.clusters.thr <- ploidy + 1

  plot_data_split <- split(scores_temp, scores_temp$geno)

  scores_temp$geno <- as.factor(scores_temp$geno)

  centers <- lapply(plot_data_split, function(x) apply(x[,3:4], 2, mean))

  if(length(centers) == 0 | length(centers) < n.clusters.thr) {
    return(list(rm= {if(length(centers) == 0) 1 else
      if(length(centers) < n.clusters.thr) 2 else 0},
      probe_name = unique(scores_temp$mks),
      n.clusters = length(centers)))
  } else {
    centers_df <- data.frame(dose = names(centers), do.call(rbind, centers))

    # If one or more cluster is missing, this code will input using the mean distance between the other clusters
    if(length(centers) != ploidy + 1){
      doses <- 0:ploidy
      input <- mean(diff(centers_df$theta)/diff(as.numeric(names(centers))))
      new_centers_df <- data.frame(dose = doses)
      new_centers_df$theta <- NA
      new_centers_df$theta[match(names(centers), doses)] <- centers_df$theta
      mis <- which(is.na(new_centers_df$theta))
      wi <- which(!is.na(new_centers_df$theta))
      for(i in 1:length(mis)){
        if(mis[i] < wi[1]) new_centers_df$theta[mis[i]] <- new_centers_df$theta[wi[1]] - (wi[1]-mis[i])*input
        if(mis[i] > wi[1]) new_centers_df$theta[mis[i]] <- new_centers_df$theta[wi[length(wi)]] + diff(c(wi[length(wi)],mis[i]))*input
      }
      centers_df <- new_centers_df
    }

    centers_df <- cbind(centers_df, cluster = 1:nrow(centers_df))
    return(list(rm = 0,
                centers_theta = centers_df$theta,
                probe_name = unique(scores_temp$mks),
                n.clusters = length(centers)))
  }
}

#' Interpolation plot
#'
#' @param scores_temp data.frame of fitpoly scores output
#' @param ploidy integer defining ploidy
#'
#' @export
plot_one_marker <- function(scores_temp, ploidy){
  plot_data_split <- split(scores_temp, scores_temp$geno)
  scores_temp$geno <- as.factor(scores_temp$geno)

  centers <- lapply(plot_data_split, function(x) apply(x[,3:4], 2, mean))
  if(length(centers) != ploidy + 1 | any(sapply(plot_data_split, nrow) < 2)) {

    p <- ggplot(scores_temp, aes(x=theta, y=R, color = geno))  +
      geom_point() + ggtitle(paste0("Marker:", unique(scores_temp$mks))) + theme_bw()

  } else {
    centers_df <- data.frame(do.call(rbind, centers))
    centers_df <- cbind(centers_df, cluster = 1:(ploidy + 1))

    p <- ggplot(scores_temp, aes(x=theta, y=R, color = geno))  +
      geom_point() +
      geom_point(data = centers_df, aes(x= theta, y = R), color = "black") +
      geom_line(data = centers_df, aes(x= theta, y = R), color = "black") + theme_bw() +
      ggtitle(paste0("Marker:", unique(scores_temp$mks)))
  }
  return(p)
}


#' Performs interpolation and returns normalized BAFs
#'
#' @param data data.frame with columns:1) MarkeName: markers IDs;2) SampleName: Samples IDs;3) X: reference allele intensities or counts;4) Y: alternative allele intensities or counts,5) R: sum of the intensities; and 6)ratio: Y/(X+Y)
#'
#' @param genos data.frame with genotype information for individuals to be used as reference for interpolation. We suggest to select for these individuals that are euploid and all have same ploidy. For array technologies, we suggest using fitpoly for obtaining the dosages, and for sequencing technologies, updog. This file has as columns: 1) MarkerName: markers IDs; 2) SampleName: Samples IDs; 3) geno: dosage.
#'
#' @param geno.pos data.frame with columns: 1) MarkerName: markers IDs; Chromosome: chromosome where the marker is located; Position: position on the chromosome the marker is located (bp).
#'
#' @param threshold.missing.geno fraction of missing genotype information allowed by marker. Markers with higher fraction are discarted
#'
#' @import dplyr
#' @import tidyr
#' @import vroom
#' @import parallel
#'
#' @export
interpolate_BAFs <- function(data = NULL,
                             genos = NULL,
                             geno.pos = NULL,
                             threshold.missing.geno=0.5,
                             ploidy.interpolation = NULL,
                             threshold.n.clusters = NULL,
                             n.cores =1,
                             baf_file_name,
                             verbose = TRUE){

  if(is.null(threshold.n.clusters)) threshold.n.clusters <- ploidy.interpolation + 1

  n.na <- genos %>% group_by(MarkerName) %>% summarize(n.na = (sum(is.na(geno))/length(geno)))

  rm.mks <- n.na$MarkerName[which(n.na$n.na > threshold.missing.geno)]

  if(verbose) print(paste0("Markers remove because of excess of missing data:", length(rm.mks)))

  if(length(rm.mks) > 0){
    genos_filt <- genos[-which(genos$MarkerName %in% rm.mks),]
  } else genos_filt <- genos

  data_id <- paste0(data$SampleName, data$MarkerName)
  genos_id <- paste0(genos$SampleName, genos$MarkerName)

  data_filt <- data[match(genos_id, data_id),]

  theta <- data_filt$ratio

  rm.na <- which(is.na(genos_filt$geno))
  if(length(rm.na) > 0){
    theta[rm.na] <- NA
  }

  data_interpolation <- data.frame(mks = data_filt$MarkerName,
                                   ind = data_filt$SampleName,
                                   theta = theta,
                                   geno = genos_filt$geno)

  lst_interpolation <- split(data_interpolation, data_interpolation$mks)

  # Generate clusters
  clust <- makeCluster(n.cores)
  #clusterExport(clust, c("par_fitpoly_interpolation", "ploidy.interpolation", "threshold.n.clusters"))
  clusterExport(clust, c("par_fitpoly_interpolation"))
  clusters <- parLapply(clust, lst_interpolation, function(x) {
    par_fitpoly_interpolation(x, ploidy= ploidy.interpolation, n.clusters.thr = threshold.n.clusters)
  })
  stopCluster(clust)

  gc()

  # Filter by number of clusters
  rm.mks <- sapply(clusters, function(x) x$rm != 0)

  if(verbose) print(paste0("Markers remove because of smaller number of clusters than set threshold:", length(rm.mks)))

  if(length(which(rm.mks)) > 0)  clusters_filt <- clusters[-which(rm.mks)] else clusters_filt <- clusters

  keep.mks <- names(clusters_filt)
  # Getting BAF for complete dataset
  theta_filt <- pivot_wider(data[which(data$MarkerName %in% keep.mks),-c(3:5)], names_from = "SampleName", values_from = "ratio")

  par <- rep(1:n.cores, each=round((nrow(theta_filt)/n.cores)+1,0))[1:nrow(theta_filt)]

  par_theta <- split.data.frame(theta_filt[,-1], par)
  par_clusters_filt <- split(clusters_filt, par)

  if(length(par_theta) < n.cores) n.cores <- length(par_theta)
  par_all <- list()
  for(i in 1:n.cores){
    par_all[[i]] <- list()
    par_all[[i]][[1]] <- par_theta[[i]]
    par_all[[i]][[2]] <- par_clusters_filt[[i]]
  }

  # Get BAF
  clust <- makeCluster(n.cores)
  #clusterExport(clust, c("get_baf", "get_baf_par", "ploidy.interpolation"))
  clusterExport(clust, c("get_baf", "get_baf_par"))
  bafs <- parLapply(clust, par_all, function(x) {
    get_baf_par(x, ploidy = ploidy.interpolation)
  })
  stopCluster(clust)

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
  bafs_df <- cbind(mks=rownames(bafs_df), bafs_df)

  # Add chr and pos info
  chr <- geno.pos$Chromosome[match(bafs_df$mks,geno.pos$MarkerName)]
  pos <- geno.pos$Position[match(bafs_df$mks,geno.pos$MarkerName)]

  baf <- cbind(Name=bafs_df$mks, Chr = chr, Position = pos, bafs_df[,-1])

  if(length(which(is.na(baf$Chr))) > 0)
    baf <- baf[-which(is.na(baf$Chr)),]

  vroom_write(baf, file = baf_file_name)
}

#' @param data data.frame with columns:1) MarkeName: markers IDs;2) SampleName: Samples IDs;3) X: reference allele intensities or counts;4) Y: alternative allele intensities or counts,5) R: sum of the intensities; and 6)ratio: Y/(X+Y)
#'
#' @param geno.pos data.frame with columns: 1) MarkerName: markers IDs; Chromosome: chromosome where the marker is located; Position: position on the chromosome the marker is located (bp).
#'
#' @import tidyr
#' @import dplyr
#'
#' @export
get_zscore <- function(data = NULL,
                       geno.pos = NULL,
                       zscore_file_name = NULL){

  zscore <- data %>% group_by(MarkerName) %>%
    mutate(z = (R - mean(R, na.rm = TRUE))/sd(R, na.rm = TRUE)) %>% select(MarkerName, SampleName, z)

  colnames(geno.pos)[1] <- c("MarkerName")

  chr <- geno.pos$Chromosome[match(zscore$MarkerName,geno.pos$MarkerName)]
  pos <- geno.pos$Position[match(zscore$MarkerName,geno.pos$MarkerName)]

  zscore <- cbind(Name=zscore$MarkerName, Chr = chr, Position = pos, zscore[,-1])

  if(length(which(is.na(zscore$Chr)))> 0)
    zscore <- zscore[-which(is.na(zscore$Chr)),]

  vroom_write(zscore, file = zscore_file_name)
}
