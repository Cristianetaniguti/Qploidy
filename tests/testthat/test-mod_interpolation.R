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

    library(vroom)
    library(dplyr)
    library(tidyr)
    session$setInputs(load_summary = list(datapath = system.file("summary_example.txt", package = "Qploidy")),
                      load_ind_names = list(datapath = system.file("ind.names_example.txt", package = "Qploidy")),
                      load_geno_pos = list(datapath = system.file("geno.pos_example.txt", package = "Qploidy")),
                      fitpoly_scores = list(datapath = system.file("fitpoly_out_306_scores.txt", package = "Qploidy")),
                      load_baf = list(datapath = system.file("baf_sub_roses_texas.txt", package = "Qploidy")),
                      load_logR = list(datapath = system.file("logR_sub_roses_texas.txt", package = "Qploidy")),
                      refs = paste0("Tetra_", 1:50),
                      ploidy = 4,
                      n.cores = 1)

    #Prepare Axiom file
    input <- list()
    input$load_summary$datapath <- system.file("summary_example.txt", package = "Qploidy")
    input$load_ind_names$datapath <- system.file("ind.names_example.txt", package = "Qploidy")
    input$load_geno_pos$datapath <- system.file("geno.pos_example.txt", package = "Qploidy")
    input$fitpoly_scores$datapath <- system.file("fitpoly_out_2780_scores.dat", package = "Qploidy")
    input$refs <- paste0("Tetra_", 1:50)
    input$ploidy <- 4
    input$n.cores <- 1

    summary <- vroom(input$load_summary$datapath)
    cleaned_summary <- clean_summary(summary_df = summary)

    expect_true(length(cleaned_summary) == 2)
    expect_true(round(sum(cleaned_summary$A_probes[,4]),0) == 3925567)

    ind.names <- vroom(input$load_ind_names$datapath)

    fitpoly_input <- summary_to_fitpoly(cleaned_summary = cleaned_summary, ind.names, geno.pos)
    R_all <- fitpoly_input[[2]]
    theta_all <- fitpoly_input[[3]]
    fitpoly_input <- fitpoly_input[[1]]

    expect_equal(round(sum(fitpoly_input$X),0), 239122013)
    expect_equal(round(sum(fitpoly_input$ratio),0), 82728)

    # Export file for fitpoly
    fitpoly_input_sele <- fitpoly_input %>% filter(SampleName %in% input$refs)

    # After fitpoly
    scores <- vroom(input$fitpoly_scores$datapath)
    expect_equal(sum(scores$geno, na.rm = TRUE), 169607)

    # Filters
    n.na <- scores %>% group_by(MarkerName) %>% summarize(n.na = (sum(is.na(geno))/length(geno))*100)
    rm.mks <- n.na$MarkerName[which(n.na$n.na > 25)]
    missing.data <- length(rm.mks)
    scores_filt <- scores[-which(scores$MarkerName %in% rm.mks),]

    keep.mks <- match(paste0(scores_filt$MarkerName, "_",scores_filt$SampleName),
                      paste0(fitpoly_input_sele$MarkerName, "_", fitpoly_input_sele$SampleName))

    fitpoly_input_filt <- fitpoly_input_sele[keep.mks,]
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
    clusterExport(clust, c("par_fitpoly_interpolation", "input"))
    clusters <- parLapply(clust, lst_interpolation, function(x) {
      library(ggplot2)
      par_fitpoly_interpolation(x, ploidy= input$ploidy, plot = FALSE)
    })
    stopCluster(clust)

    # Plot markers interpolation
    p <- plot_one_marker(lst_interpolation[[50]], ploidy = 4)
    # Discarded marker:
    p <- plot_one_marker(lst_interpolation[[1]], ploidy = 4)

    # Filter by number of clusters
    rm.mks <- sapply(clusters, function(x) is.null(x$mod))
    wrong_n_clusters <- sum(rm.mks)
    clusters_filt <- clusters[-which(rm.mks)]

    # Filtered markers table
    filtered.markers <- data.frame(nrow(cleaned_summary[[1]]),
                                   missing.data,
                                   wrong_n_clusters,
                                   length(clusters_filt))

    colnames(filtered.markers) <- c("n.mk.start", "missing.data", "wrong_n_clusters", "n.mk.selected")
    expect_equal(c(nrow(cleaned_summary[[1]]),
                   missing.data,
                   wrong_n_clusters,
                   length(clusters_filt)), c(2379, 484, 1824, 71))

    keep.mks <- names(clusters_filt)
    # Getting logR and BAF for entire dataset
    R_filt <- R_all[match(keep.mks, R_all$MarkerName),]
    theta_filt <- theta_all[match(keep.mks, theta_all$MarkerName),]

    input$n.cores <- 2
    par <- rep(1:input$n.cores, each=round((nrow(R_filt)/input$n.cores)+1,0))[1:nrow(R_filt)]

    par_R <- split.data.frame(R_filt[,-1], par)
    par_theta <- split.data.frame(theta_filt[,-1], par)
    par_clusters_filt <- split(clusters_filt, par)

    par_all <- list()
    for(i in 1:input$n.cores){
      par_all[[i]] <- list()
      par_all[[i]][[1]] <- par_R[[i]]
      par_all[[i]][[2]] <- par_theta[[i]]
      par_all[[i]][[3]] <- par_clusters_filt[[i]]
    }

    clust <- makeCluster(input$n.cores)
    clusterExport(clust, c("get_logR", "get_logR_par", "input"))
    logRs_diplo <- parLapply(clust, par_all, function(x) {
      get_logR_par(x, ploidy = input$ploidy)
    })
    stopCluster(clust)

    logRs_diplod_lt <- unlist(logRs_diplo, recursive = F)
    logRs_diplod_m <- do.call(rbind, logRs_diplod_lt)
    rownames(logRs_diplod_m) <- keep.mks
    logRs_diplod_m <- cbind(mks=rownames(logRs_diplod_m), logRs_diplod_m)

    # Get BAF
    clust <- makeCluster(input$n.cores)
    clusterExport(clust, c("get_baf", "get_baf_par", "input"))
    bafs_diplo <- parLapply(clust, par_all, function(x) {
      get_baf_par(x, ploidy = input$ploidy)
    })
    stopCluster(clust)

    bafs_diplo_lt <- unlist(bafs_diplo, recursive = F)
    bafs_diplo_m <- do.call(rbind, bafs_diplo_lt)
    rownames(bafs_diplo_m) <- keep.mks
    colnames(bafs_diplo_m) <- colnames(logRs_diplod_m)[-1]
    bafs_diplo_df <- as.data.frame(bafs_diplo_m)
    bafs_diplo_df <- cbind(mks=rownames(bafs_diplo_df), bafs_diplo_df)

    expect_equal(sum(bafs_diplo_df[,-1]), 2628.072)
    expect_equal(round(sum(logRs_diplod_m[,-1]),0), -921)

    # geno.pos <- vroom(input$load_geno_pos$datapath)
    # geno.pos <- as.data.frame(geno.pos)
    #
    # new.id <- paste0(geno.pos[,2], "_", geno.pos[,3])
    # fitpoly_input$MarkerName <- geno.pos[match(fitpoly_input$MarkerName, geno.pos[,1]),]
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

