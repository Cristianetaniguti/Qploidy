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

    session$setInputs(summary = list(datapath = system.file("fitpoly_input.txt", package = "Qploidy")),
                      refs = paste0("Tetra_", 1:50),
                      ploidy = 4,
                      n.cores = 1)

    input <- list()
    input$summary$datapath <- system.file("fitpoly_input.txt", package = "Qploidy")
    input$refs <- paste0("Tetra_", 1:50)
    input$ploidy <- 4
    input$n.cores <- 1

    fitpoly_input <- vroom(input$summary$datapath)
    refs_fitpoly_inputs <- as.data.frame(fitpoly_input[which(fitpoly_input$SampleName %in% input$refs),])

    library(fitPoly)
    out <- sample(1:1000, 1)
    saveMarkerModels(ploidy= input$ploidy,
                     data=refs_fitpoly_inputs,
                     p.threshold=0.1,
                     filePrefix= paste0("fitpoly_out_",out),
                     ncores=input$n.cores)

    scores <- vroom(paste0("fitpoly_out_",out, "_scores.dat"))
    expect_true(input$x == 1)

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

