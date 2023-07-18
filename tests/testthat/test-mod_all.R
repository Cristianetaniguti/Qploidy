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

    # Parei aqui! Testar os graficos
    baf <- system.file("baf.example.txt", package = "Qploidy")
    logR <- system.file("logR.example.txt", package = "Qploidy")


    # Test break counts
    ## from mappoly
    aneuploids <- vroom(system.file("aneuploids.ex.txt", package = "Qploidy"))
    load(system.file("mappoly.homoprob.ex.RData", package = "Qploidy"))
    p_m1 <- plot(mappoly.homoprob, ind = 1)
    p_m2 <- count_breaks_mappoly(mappoly.homoprob$homoprob, aneuploids, by_LG = FALSE)

    ## from polyOrigin
    f1.codes <- vroom(system.file("F1codes.polyorigin.txt", package = "Qploidy"))
    df  <- vroom(system.file("genofile_sub.csv", package = "Qploidy"))
    homoprob <- get_probs_polyorigin(df,
                                     f1.codes = f1.codes,
                                     ploidy = 4, n.cores = 2)

    p_p1 <- plot(x = homoprob, lg = 1, ind = 2)
    p_p2 <- count_breaks_mappoly(homoprob = homoprob$homoprob, aneuploids = aneuploids, by_LG = FALSE)

    ## Check with diaQTL
    # library(diaQTL)
    # data <- read_data(genofile = "genofile_sub.csv",
    #                   ploidy = 4,
    #                   pedfile = "pedfile_sub.csv",
    #                   n.core = 2)
    #
    # p_d <- haplo_plot(data = data,
    #                 id = "16400_N080",
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

