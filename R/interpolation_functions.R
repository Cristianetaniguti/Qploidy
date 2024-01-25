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
par_fitpoly_interpolation <- function(scores_temp, ploidy, plot=TRUE, n.clusters.thr = 0){

  if(is.null(n.clusters.thr)) n.clusters.thr <- ploidy + 1

  plot_data_split <- split(scores_temp, scores_temp$geno)

  scores_temp$geno <- as.factor(scores_temp$geno)

  centers <- lapply(plot_data_split, function(x) apply(x[,3:4], 2, mean))
  #if(length(centers) != ploidy + 1 | any(sapply(plot_data_split, nrow) < 2)| any(is.na(sapply(centers, "[[", 2)))) {
  # scores_temp$theta <- (scores_temp$theta - min(scores_temp$theta, na.rm = TRUE))/(max(theta, na.rm = TRUE) -min(scores_temp$theta, na.rm = TRUE))

  if(length(centers) == 0 | length(centers) < n.clusters.thr) {
    if(plot){
      p <- ggplot(scores_temp, aes(x=theta, y=R, color = geno))  +
        geom_point()
    } else p <- NULL

    return(list(rm= {if(length(centers) == 0) 1 else
      if(length(centers) < n.clusters.thr) 2 else 0},
      probe_name = unique(scores_temp$mks),
      plot=p,
      n.clusters = length(centers)))
  } else {
    centers_df <- data.frame(dose = names(centers), do.call(rbind, centers))

    # If one or more cluster is missing, this code will input using the mean distance between the other clusters
    if(length(centers) != ploidy + 1){
      doses <- 0:ploidy
      input <- mean(diff(centers_df$theta)/diff(as.numeric(names(centers))))
      new_centers_df <- data.frame(dose = doses)
      new_centers_df$R <- NA
      new_centers_df$R[match(names(centers), doses)] <- centers_df$R
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

    if(plot){
      p <- ggplot(scores_temp, aes(x=.data$theta, y=.data$R, color = .data$geno))  +
        geom_point() +
        geom_point(data = centers_df, aes(x= .data$theta, y = .data$R), color = "black") +
        geom_line(data = centers_df, aes(x= .data$theta, y = .data$R), color = "black") + theme_bw()
    } else p <- NULL

    return(list(rm = 0,
                centers_theta = centers_df$theta,
                certers_R = centers_df$R,
                plot = p,
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
