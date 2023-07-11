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
             radioButtons(ns("example_data"), label = "Choose example data set",
                          choices = c("Roses Texas" = "roses_texas",
                                      "Roses France" = "roses_france",
                                      "Potatoes Texas" = "potatoes"),
                          selected = "roses_texas")
      )

    )
  )
}

#' upload Server Functions
#'
#' @import vroom
#'
#' @noRd
mod_upload_server <- function(input, output, session, parent_session){
  ns <- session$ns

  # Reset buttons
  values <- reactiveValues(
    upload_state_summary = 0,
    upload_state_logr_baf = 0
  )

  observeEvent(input$reset_all, {
    values$upload_state_summary = 0
    values$upload_state_logr_baf = 0
  })

  observeEvent(input$submit_summary, {
    values$upload_state_summary <- 'uploaded'
  })

  observeEvent(input$submit_logr_baf, {
    values$upload_state_logr_baf <- 'uploaded'
  })



  input_logr_baf <- reactive({
    if(values$upload_state_logr_baf == 0 | values$upload_state_summary != 0){
      return(NULL)
    } else if (values$upload_state_logr_baf == "uploaded"){
      validate(
        need(!is.null(input$load_logr) | !is.null(input$load_baf), "Upload both logR and BAF files before submit")
      )
      logR <- vroom(input$load_logR$datapath)
      baf <- vroom(input$load_baf$datapath)

      return(list(logR, baf))
    }
  })

  loadExample <- reactive({
    if(is.null(input_logr_baf()) &
       is.null(input_summary())){
      print(input$example_data)
      if(input$example_data == "roses_texas"){
        baf <- vroom(system.file("baf_sub_roses_texas.txt", package = "Qploidy"))
        logR <- vroom(system.file("logr_sub_roses_texas.txt", package = "Qploidy"))
        summary <- vroom(system.file("fitpoly_input.txt", package = "Qploidy"))
      } else if(input$example_data == "roses_france"){
        logR <- vroom(system.file("baf_roses_texas.txt", package = "Qploidy"))
        baf <- vroom(system.file("logr_roses_texas.txt", package = "Qploidy"))
      } else if(input$example_data == "potatoes") {
        logR <- vroom(system.file("baf_roses_texas.txt", package = "Qploidy"))
        baf <- vroom(system.file("logr_roses_texas.txt", package = "Qploidy"))
      }
      data <- list(logR, baf, summary)
      data
    } else NULL
  })

  logR_BAF = reactive({
    if(is.null(input_logr_baf())){
      return(loadExample()[1:2])
    } else return(input_logr_baf())
  })

  summary = reactive({
    if(is.null(input_summary())){
      return(loadExample()[[3]])
    } else return(input_summary())
  })

  return(list(summary = reactive(summary()),
              logR_BAF = reactive(logR_BAF())))
}

## To be copied in the UI
# mod_upload_ui("upload_1")

## To be copied in the server
# mod_upload_server("upload_1")
