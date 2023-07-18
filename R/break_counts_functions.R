#' Count number of recombination breakpoints in parents gametes
#'
#' @param object of class mappoly.homoprob
#' @param n.graphics define number of graphics
#' @param ncol define number of columns
#'
#' @export
count_breaks_mappoly <- function(homoprob, aneuploids = NULL, n.graphics=NULL, ncol=NULL, by_LG = TRUE){

  homoprob$LG <- paste0("Chr",homoprob$LG)
  homoprob$most_likely <- homoprob$probability
  homoprob$most_likely[which(homoprob$probability > 0.5)] <- 1
  homoprob$most_likely[which(homoprob$probability <= 0.5)] <- 0
  homoprob$parent <- NA
  ploidy <- length(unique(homoprob$homolog))/2
  homoprob$parent[which(homoprob$homolog %in% letters[1:ploidy])] <- "P1"
  homoprob$parent[which(homoprob$homolog %in% letters[(ploidy+1):(2*ploidy)])] <- "P2"

  counts <- homoprob %>% group_by(individual, LG, parent, homolog) %>%
    summarize(counts = sum(sequence(rle(as.character(most_likely))$length) == 1) - 1) %>%
    group_by(individual, LG, parent) %>%
    summarize(total_counts = sum(counts))

  if(!is.null(aneuploids)){
    aneuploids <- melt(aneuploids)
    colnames(aneuploids) <- c("individual", "LG", "ploidy")
    counts <- merge(counts, aneuploids, by = c("individual", "LG"))
  }

  p <- list()
  n.ind <- length(unique(counts$individual))

  if(is.null(n.graphics) & is.null(ncol)){
    if(n.ind/25 <= 1) {
      n.graphics = 1
      ncol=1
    }else {
      n.graphics = round(n.ind/25,0)
      ncol=round(n.ind/25,0)
    }
  }

  size <-n.ind
  if(size%%n.graphics == 0){
    div.n.graphics <- rep(1:n.graphics, each= size/n.graphics)
  } else {
    div.n.graphics <-   c(rep(1:n.graphics, each = round(size/n.graphics,0)),rep(n.graphics, size%%n.graphics))
  }

  div.n.graphics <- div.n.graphics[1:n.ind]
  div.n.graphics <- rep(div.n.graphics, each = 2*length(unique(counts$LG)))

  counts$individual <- factor(as.character(counts$individual), levels = sort(as.character(unique(counts$individual))))
  counts$LG <- as.factor(counts$LG)
  counts$ploidy <- as.factor(counts$ploidy)

  ploidycolors <- brewer.pal(9, "Set1")[1:length(unique(counts$ploidy))]
  mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(counts$LG)))
  set.seed(20)
  mycolors <- sample(mycolors)
  ploidycolors <- sample(ploidycolors)
  names(ploidycolors) <- levels(counts$ploidy)
  names(mycolors) <- levels(counts$LG)

  p_list <- counts %>% ungroup() %>%  mutate(div.n.graphics = div.n.graphics) %>%
    split(., .$div.n.graphics) %>%
    lapply(., function(x) ggplot(x, aes(x=paste0(individual, "_", parent), y=total_counts, fill = LG)) +
             {if(by_LG) geom_bar(stat = "identity", aes(fill = LG)) else geom_bar(stat = "identity", aes(fill = ploidy))}+
             coord_flip() +
             {if(by_LG) scale_fill_manual(values=mycolors) else scale_fill_manual(values=ploidycolors)} +
             theme(axis.title.y = element_blank(),
                   axis.title.x = element_blank(),
                   axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                   legend.key.size = unit(1, 'cm')) +
             {if(!by_LG) facet_grid(. ~ LG)} +
             {if(by_LG) labs(fill="groups") else labs(fill="ploidy")})

  p <- ggarrange(plotlist = p_list, common.legend = T, ncol = ncol, nrow = round(n.graphics/ncol,0))

  return(p)
}


#' Plots mappoly.homoprob from MAPpoly
#'
#' @param x an object of class \code{mappoly.homoprob}
#'
#' @param stack logical. If \code{TRUE}, probability profiles of all homologues
#'              are stacked in the plot (default = FALSE)
#'
#' @param lg indicates which linkage group should be plotted. If \code{NULL}
#'           (default), it plots the first linkage group. If
#'           \code{"all"}, it plots all linkage groups
#'
#' @param ind indicates which individuals should be plotted. It can be the
#'            position of the individuals in the dataset or it's name.
#'            If \code{NULL} (default), the function plots the first
#'            individual
#'
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
#'
#' @param ... unused arguments
#'
#' @import ggplot2
#'
#' @export
plot.mappoly.homoprob <- function(x, stack = FALSE, lg = NULL,
                                  ind = NULL,
                                  verbose = TRUE, ...){

  all.ind <- as.character(unique(x$homoprob$individual))
  #### Individual handling ####
  if(length(ind) > 1){
    if (verbose) message("More than one individual provided: using the first one")
    ind <- ind[1]
  }
  if(is.null(ind)){
    ind <- as.character(all.ind[1])
    df.pr1 <- subset(x$homoprob, individual  ==  ind)
  } else if(is.numeric(ind)) {
    if(ind > length(all.ind))
      stop("Please chose an individual number between 1 and ", length(all.ind))
    ind <- as.character(all.ind[ind])
    df.pr1 <- subset(x$homoprob, individual  ==  ind)
  } else if (is.character(ind)){
    if(!ind%in%all.ind)
      stop(safeError("Invalid individual name"))
  } else stop(safeError("Invalid individual name"))

  #### LG handling ####
  if(is.null(lg))
    lg <- 1
  if(all(lg == "all"))
    lg <- unique(x$homoprob$LG)
  LG <- individual <- map.position <- probability <- homolog <- NULL
  if(length(lg) > 1 & !stack)
  {
    if (verbose) message("Using 'stack = TRUE' to plot multiple linkage groups")
    stack <- TRUE
  }
  if(stack){
    ##subset linkage group
    if(!is.null(lg)){
      df.pr1 <- subset(x$homoprob, LG%in%lg)
      df.pr1 <- subset(df.pr1, individual  ==  ind)
    } else
      df.pr1 <- subset(x$homoprob, individual  ==  ind)
    p <- ggplot(df.pr1, aes(x = map.position, y = probability, fill = homolog, color  = homolog)) +
      geom_density(stat = "identity", alpha = 0.7, position = "stack") +
      ggtitle(ind) +
      facet_grid(rows = vars(LG)) +
      ylab(label = "Homologs probabilty") +
      xlab(label = "Map position") +  theme_minimal()
  } else {
    ##subset linkage group
    if(is.null(lg)){
      lg <- 1
      df.pr1 <- subset(x$homoprob, LG %in% lg)
    } else df.pr1 <- subset(x$homoprob, LG %in% lg)
    df.pr1 <- subset(df.pr1, individual  ==  ind)
    p <- ggplot(df.pr1, aes(x = map.position, y = probability, fill = homolog, color  = homolog)) +
      geom_density(stat = "identity", alpha = 0.7) +
      ggtitle(paste(ind, "   LG", lg)) +
      facet_grid(rows = vars(homolog)) +
      theme_minimal() +
      ylab(label = "Homologs probabilty") +
      xlab(label = "Map position")
  }
  return(p)
}
