testServer(
  mod_all_server,
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
    # packages
    # library(vroom)
    # library(tidyr)

    session$setInputs(samples = c("Diplo_1", "Diplo_2","Tetra_1","Tetra_2","Unknow_1"),
                      ploidys = c(2,5),
                      area = 0.75,
                      filter_corr = 0,
                      filter_diff = 0,
                      graphics = c("Unknow_1"),
                      area_single = 0.75,
                      dot.size = 1,
                      colors = TRUE,
                      add_lines = TRUE,
                      ploidy = 4)

    #inputs
    # input <- list()
    # input$samples <- c("Diplo_1", "Diplo_2","Tetra_1","Tetra_2","Unknow_1")
    # input$ploidys <- c(2,5)
    # input$area <- 0.75
    # input$filter_corr <- 0
    # input$filter_diff <- 0
    # input$graphics <- c("Unknow_1")
    # input$area_single <- 0.75
    # input$dot.size <- 1
    # input$colors <- TRUE
    # input$add_lines <- TRUE
    # input$ploidy <- 4

    # upload files
    baf <- vroom(system.file("baf.example.txt", package = "Qploidy"), show_col_types = FALSE)
    logR <- vroom(system.file("logR.example.txt", package = "Qploidy"), show_col_types = FALSE)

    logR_baf <- list(logR, baf)

    # Get overall ploidy estimation tables
    data_sample <- logR_baf[[2]][,c(2,3,which(colnames(logR_baf[[2]]) %in% c(input$samples)))]

    data_sample <- data_sample[order(data_sample$Chr, data_sample$Position),]

    est.ploidy.chr_df <- area_estimate_ploidy_by_chr(data_sample, ploidys = input$ploidys, area = input$area)

    # Get overall ploidy estimation graphics
    ## Aneuploid individuals
    ps <- plots_overall(est.ploidy.chr_df,
                        filter_diff = input$filter_diff,
                        filter_corr = input$filter_corr)

    # Overal break counts
    ## Load
    aneuploids <- vroom(system.file("aneuploids.ex.txt", package = "Qploidy"), show_col_types = FALSE)
    aneuploids$X <- c("Unknow_1", "Unknow_3")

    temp <- load(system.file("mappoly.homoprob.ex.RData", package = "Qploidy"))
    haplo_mappoly <- get(temp)
    polyorigin  <- vroom(system.file("genofile_sub.csv", package = "Qploidy"), show_col_types = FALSE)
    f1.codes <- vroom(system.file("F1codes.polyorigin.txt", package = "Qploidy"), show_col_types = FALSE)
    haplo_polyorigin <- get_probs_polyorigin(polyorigin,
                                             f1.codes = f1.codes,
                                             ploidy = 4, n.cores = 2)

    #     haplo_polyorigin <- Qploidy:::get_probs_polyorigin_sd(polyorigin,
    #                                                           f1.codes = f1.codes,
    #                                                           ploidy = 4, n.cores = 2)

    ## MAPpoly
    p_m_df <- count_breaks_df(homoprob = haplo_mappoly$homoprob,
                              inds = input$samples,
                              aneuploids = aneuploids)

    ## polyOrigin
    p_p_df <- count_breaks_df(homoprob = haplo_polyorigin$homoprob, aneuploids = aneuploids)

    # Single individual analysis
    data_sample <- logR_baf[[2]][,c(2,3,which(colnames(logR_baf[[2]]) %in% c(input$graphics)))]
    colnames(data_sample)[3] <- "sample"

    p_baf <- plot_baf(data_sample, input$area_single, input$ploidy, input$dot.size, input$add_lines, input$colors)

    p_hist <- plot_baf_hist(data_sample, input$area_single, input$ploidy, input$colors, input$add_lines)

    if(any(unique(haplo_mappoly$homoprob$individual) %in% input$graphics)){
      haplo_lst <- list()
      for(i in 1:length(unique(haplo_mappoly$homoprob$LG))){
        haplo_lst[[i]] <- plot(haplo_mappoly, lg = unique(haplo_mappoly$homoprob$LG)[i],
                               ind = input$graphics,
                               use.plotly = FALSE)
      }

      all_haplo_mappoly <- ggarrange(plotlist = haplo_lst,ncol = 1, common.legend = TRUE)
    } else all_haplo_mappoly <- NULL

    if(any(unique(haplo_polyorigin$homoprob$individual) %in% input$graphics)){
      haplo_lst <- list()
      for(i in 1:length(unique(haplo_polyorigin$homoprob$LG))){
        haplo_lst[[i]] <- plot(haplo_polyorigin, lg = unique(haplo_polyorigin$homoprob$LG)[i],
                               ind = input$graphics,
                               use.plotly = FALSE)
      }

      all_haplo_polyorigin <- ggarrange(plotlist = haplo_lst,ncol = 1, common.legend = TRUE)
    } else all_haplo_polyorigin <- NULL

    ## Check with diaQTL
    # library(diaQTL)
    # data <- read_data(genofile = "genofile_sub.csv",
    #                   ploidy = 4,
    #                   pedfile = "pedfile_sub.csv",
    #                   n.core = 2)
    #
    # p_d <- haplo_plot(data = data,
    #                 id = "Unknow_1",
    #                 chrom = 1,
    #                 position = "cM")


  })

test_that("module ui works", {
  ui <- mod_all_ui(id = "test")
  golem::expect_shinytaglist(ui)
  # Check that formals have not been removed
  fmls <- formals(mod_all_ui)
  for (i in c("id")){
    expect_true(i %in% names(fmls))
  }
})

