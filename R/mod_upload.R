#' upload UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_upload_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
             radioButtons(ns("load"), label = h3("Choose data set"),
                          choices = list("PennCNV first try" = "PennCNV_first",
                                         "Axiom CNV Sumary tool" = "Axiom_CNV",
                                         "PennCNV second try" = "PennCNV_second",
                                         "In house logR and BAF - diploids" = "in_house_diploids",
                                         "In house logR and BAF - tetraploids" = "in_house_tetraploids",
                                         "In house logR and BAF - others" = "in_house_others",
                                         "In house logR and BAF - tetra model - diploids" = "in_house_tetra_diploids",
                                         "In house logR and BAF - tetra model - tetraploids" = "in_house_tetra_tetraploids",
                                         "In house logR and BAF - tetra model - others" = "in_house_tetra_others",
                                         "In house logR and BAF - tetra model - potato" = "potato"),
                          selected = "PennCNV_first"),
      )
    )

  )
}

#' upload Server Functions
#'
#' @noRd
mod_upload_server <- function(input, output, session, data, parent_session){
  ns <- session$ns
  data <- eventReactive(input$load, {
    # Read file from PennCNV
    if(input$load == "PennCNV_first"){
      baf <- vroom("data/joint_logR_BAF/baf.txt")
      logR <- vroom("data/joint_logR_BAF/logR.txt")
      # Read file from Axiom CNV Summary Tool
    } else if(input$load == "Axiom_CNV"){
      logR <- vroom("data/Axiom_graph_results/logR2.txt")
      baf <- vroom("data/Axiom_graph_results/baf3.txt")
      # Read file from PennCNV second try
    } else if(input$load == "PennCNV_second") {
      logR <- vroom("data/script_second_round/logR.txt")
      baf <- vroom("data/script_second_round/baf.txt")
    } else if(input$load == "in_house_diploids"){
      logR <- vroom("data/interpolation_results/logR_diploids_filt.txt")
      baf <- vroom("data/interpolation_results/baf_diploids_filt.txt")
    } else if(input$load == "in_house_tetraploids"){
      logR <- vroom("data/interpolation_results/logR_tetra_filt.txt")
      baf <- vroom("data/interpolation_results/baf_tetra_filt.txt")
    } else if(input$load == "in_house_others"){
      logR <- vroom("data/interpolation_results/logR_others_filt.txt")
      baf <- vroom("data/interpolation_results/baf_others_filt.txt")
    } else if(input$load == "in_house_tetra_diploids"){
      logR <- vroom("data/interpolation_results/logR_tetra_diploids_filt2.txt")
      baf <- vroom("data/interpolation_results/baf_tetra_diploids_filt2.txt")
    } else if(input$load == "in_house_tetra_tetraploids"){
      logR <- vroom("data/interpolation_results/logR_tetra_tetra_filt.txt")
      baf <- vroom("data/interpolation_results/baf_tetra_tetra_filt.txt")
    } else if(input$load == "in_house_tetra_others"){
      logR <- vroom("data/interpolation_results/logR_tetra_others_filt.txt")
      baf <- vroom("data/interpolation_results/baf_tetra_others_filt.txt")
    } else if(input$load == "potato"){
      logR <- vroom("data/interpolation_results/logR_tetra_potato_filt_genos2.txt")
      baf <- vroom("data/interpolation_results/baf_tetra_potato_filt_genos2.txt")
    }
    data <- list(logR, baf)
    data
  })

  return(list(data = reactive(data())))
}

## To be copied in the UI
# mod_upload_ui("upload_1")

## To be copied in the server
# mod_upload_server("upload_1")
