testServer(
  mod_upload_server,
  # Add here your module params
  args = list()
  , {
    library(vroom)
    # Here are some examples of tests you can
    # run on your module
    # - Testing the setting of inputs
    session$setInputs(load_summary = list(datapath = system.file("summary_example.txt", package = "Qploidy")),
                      load_ind_names = list(datapath = system.file("ind.names_example.txt", package = "Qploidy")),
                      load_geno_pos = list(datapath = system.file("geno.pos_example.txt", package = "Qploidy")),
                      load_baf = list(datapath = system.file("baf_sub_roses_texas.txt", package = "Qploidy")),
                      load_logR = list(datapath = system.file("logR_sub_roses_texas.txt", package = "Qploidy")),
                      example = "roses_texas")

    #Prepare Axiom file
    input <- list()
    input$load_summary$datapath <- system.file("summary_example.txt", package = "Qploidy")
    input$load_ind_names$datapath <- system.file("ind.names_example.txt", package = "Qploidy")
    input$load_geno_pos$datapath <- system.file("geno.pos_example.txt", package = "Qploidy")

    summary <- vroom(input$load_summary$datapath)
    cleaned_summary <- clean_summary(summary_df = summary)

    expect_true(length(cleaned_summary) == 2)
    expect_true(round(sum(cleaned_summary$A_probes[,4]),0) == 3998947)

    ind.names <- vroom(input$load_ind_names$datapath)
    geno.pos <- vroom(input$load_geno_pos$datapath)

    fitpoly_input <- summary_to_fitpoly(cleaned_summary = cleaned_summary, ind.names, geno.pos)

    expect_true(round(sum(fitpoly_input$X),0) == 108290827)
    expect_true(round(sum(fitpoly_input$ratio),0) == 35339)
})

test_that("module ui works", {
  ui <- mod_upload_ui(id = "test")
  golem::expect_shinytaglist(ui)
  # Check that formals have not been removed
  fmls <- formals(mod_upload_ui)
  for (i in c("id")){
    expect_true(i %in% names(fmls))
  }
})

