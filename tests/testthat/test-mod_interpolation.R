testServer(
  mod_interpolation_server,
  # Add here your module params
  args = list()
  , {
    ns <- session$ns
    expect_true(
      inherits(ns, "function")
    )
    expect_true(
      grepl(id, ns(""))
    )
    expect_true(
      grepl("test", ns("test"))
    )

    session$setInputs(load_summary = list(datapath = system.file("summary_example.txt", package = "Qploidy")),
                      load_ind_names = list(datapath = system.file("ind.names_example.txt", package = "Qploidy")),
                      load_geno_pos = list(datapath = system.file("geno.pos_example.tsv.gz", package = "Qploidy")),
                      fitpoly_scores = list(datapath = system.file("tetraploids_refs_scores.tsv.gz", package = "Qploidy")),
                      zscores = list(datapath = system.file("zscore_example.tsv.gz", package = "Qploidy")),
                      refs = paste0("Tetra_", 1:50),
                      ploidy = 4,
                      n.cores = 1)

    # Prepare Axiom file
    # library(vroom)
    # library(dplyr)
    # library(tidyr)
    # library(testthat)
    # library(Qploidy)
    #
    # input <- list()
    # input$load_summary$datapath <- system.file("summary_example.txt", package = "Qploidy")
    # input$load_ind_names$datapath <- system.file("ind.names_example.txt", package = "Qploidy")
    # input$load_geno_pos$datapath <- system.file("geno.pos_example.tsv.gz", package = "Qploidy")
    # input$fitpoly_scores$datapath <- system.file("tetraploids_refs_scores.tsv.gz", package = "Qploidy")
    # input$zscore$datapath <- system.file("zscore_example.tsv.gz", package = "Qploidy")
    # input$refs <- paste0("Tetra_", 1:50)
    # input$ploidy <- 4
    # input$n.cores <- 1

    summary <- vroom(input$load_summary$datapath, show_col_types = FALSE)
    cleaned_summary <- clean_summary(summary_df = summary)

    expect_true(length(cleaned_summary) == 2)
    expect_true(round(sum(cleaned_summary$A_probes[,4]),0) == 3925567)

    ind.names <- vroom(input$load_ind_names$datapath, show_col_types = FALSE)

    R_theta <- get_R_theta(cleaned_summary, ind.names)

    R_all <- R_theta[[1]]
    theta_all <- R_theta[[2]]

    fitpoly_input <- summary_to_fitpoly(R_all, theta_all)

    expect_equal(round(sum(fitpoly_input$X),0), 239122013)
    expect_equal(round(sum(fitpoly_input$ratio),0), 82226)

    # Export file for fitpoly
    fitpoly_input_sele <- fitpoly_input %>% filter(SampleName %in% input$refs)

    # library(fitPoly)
    # saveMarkerModels(ploidy=4,
    #                  data=fitpoly_input_sele,
    #                  p.threshold=0.9,
    #                  filePrefix="tetraploids_refs",
    #                  ncores=20)

    # After fitpoly
    scores <- vroom(input$fitpoly_scores$datapath, show_col_types = FALSE)
    expect_equal(sum(scores$geno, na.rm = TRUE), 203403)

    # Filters
    n.na <- scores %>% group_by(MarkerName) %>% summarize(n.na = (sum(is.na(geno))/length(geno))*100)
    rm.mks <- n.na$MarkerName[which(n.na$n.na > 25)]
    missing.data <- length(rm.mks)
    scores_filt <- scores[-which(scores$MarkerName %in% rm.mks),]

    keep.mks <- match(paste0(scores_filt$SampleName, "_", scores_filt$MarkerName),
                      paste0(fitpoly_input$SampleName, "_", fitpoly_input$MarkerName))

    fitpoly_input_filt <- fitpoly_input[keep.mks,]
    theta <- fitpoly_input_filt$ratio
    R <- fitpoly_input_filt$R

    rm.na <- which(is.na(scores_filt$geno))
    if(length(rm.na) > 0){
      R[rm.na] <- NA
      theta[rm.na] <- NA
    }

    data_interpolation <- data.frame(mks = fitpoly_input_filt$MarkerName,
                                     ind = fitpoly_input_filt$SampleName,
                                     R = R,
                                     theta = theta,
                                     geno = scores_filt$geno)

    lst_interpolation <- split(data_interpolation, data_interpolation$mks)

    # Generate clusters
    library(parallel)
    clust <- makeCluster(input$n.cores)
    ploidy <- input$ploidy
    clusterExport(clust, c("par_fitpoly_interpolation"))
    #clusterExport(clust, c("par_fitpoly_interpolation", "ploidy"))
    clusters <- parLapply(clust, lst_interpolation, function(x) {
      library(ggplot2)
      par_fitpoly_interpolation(x, ploidy= ploidy)
    })
    stopCluster(clust)

    # Filter by number of clusters
    rm.mks <- sapply(clusters, function(x) x$rm != 0)
    wrong_n_clusters <- sum(rm.mks)
    clusters_filt <- clusters[-which(rm.mks)]

    # Filtered markers table
    filtered.markers <- data.frame(n.mk.start = nrow(cleaned_summary[[1]]),
                                   missing.data,
                                   wrong_n_clusters,
                                   n.mk.selected = length(clusters_filt))

    expect_equal(c(nrow(cleaned_summary[[1]]),
                   missing.data,
                   wrong_n_clusters,
                   length(clusters_filt)), c(2379, 294, 741, 1332))

    mks <- names(clusters)
    mks[which(rm.mks)] <- paste(mks[which(rm.mks)], "(discarded)")

    # Add discarted
    data_interpolation_disc.mks <- do.call(rbind, lst_interpolation)
    data_interpolation_disc.mks$mks.disc <- mks[match(data_interpolation_disc.mks$mks, names(clusters))]

    # Plot markers interpolation
    p <- plot_one_marker(lst_interpolation[[1]], ploidy = 4)
    # Discarded marker:
    p <- plot_one_marker(lst_interpolation[[2]], ploidy = 4)


    keep.mks <- names(clusters_filt)
    # Getting logR and BAF for entire dataset
    theta_filt <- theta_all[match(keep.mks, theta_all$MarkerName),]

    par <- rep(1:input$n.cores, each=round((nrow(theta_filt)/input$n.cores)+1,0))[1:nrow(theta_filt)]

    par_theta <- split.data.frame(theta_filt[,-1], par)
    par_clusters_filt <- split(clusters_filt, par)

    par_all <- list()
    for(i in 1:input$n.cores){
      par_all[[i]] <- list()
      par_all[[i]][[1]] <- par_theta[[i]]
      par_all[[i]][[2]] <- par_clusters_filt[[i]]
    }

    # Get BAF
    clust <- makeCluster(input$n.cores)
    clusterExport(clust, c("get_baf", "get_baf_par"))
    #clusterExport(clust, c("get_baf", "get_baf_par", "ploidy"))
    bafs_diplo <- parLapply(clust, par_all, function(x) {
      get_baf_par(x, ploidy = ploidy)
    })
    stopCluster(clust)

    bafs_diplo_lt <- unlist(bafs_diplo, recursive = F)
    bafs_diplo_m <- do.call(rbind, bafs_diplo_lt)
    rownames(bafs_diplo_m) <- keep.mks
    colnames(bafs_diplo_m) <- colnames(theta_filt)[-1]
    bafs_diplo_df <- as.data.frame(bafs_diplo_m)
    bafs_diplo_df <- cbind(mks=rownames(bafs_diplo_df), bafs_diplo_df)

    expect_equal(sum(bafs_diplo_df[,-1]), 25343, tolerance = 1)

    geno.pos <- vroom(input$load_geno_pos$datapath, show_col_types = FALSE)

    chr <- geno.pos$Chr[match(bafs_diplo_df$mks,geno.pos$Name)]
    pos <- geno.pos$Position[match(bafs_diplo_df$mks,geno.pos$Name)]

    baf <- cbind(Name=bafs_diplo_df$mks, Chr = chr, Position = pos, bafs_diplo_df[,-1])

    colnames(geno.pos)[2] <- "Chromosome"
    zscore <- get_zscore(data = fitpoly_input, geno.pos = geno.pos)
  })

test_that("module ui works", {
  ui <- mod_interpolation_ui(id = "test")
  golem::expect_shinytaglist(ui)
  # Check that formals have not been removed
  fmls <- formals(mod_interpolation_ui)
  for (i in c("id")){
    expect_true(i %in% names(fmls))
  }
})

