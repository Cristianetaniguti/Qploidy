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
             fileInput("load_summary", label = "Upload Axiom summary file"), hr(),
             p("or"),
             fileInput("load_baf", label = "Upload BAF file"),
             fileInput("load_logR", label = "Upload logR file"), hr(),
             p("or"),
             radioButtons(ns("example"), label = "Choose example data set",
                          choices = list("Roses Texas" = "roses_texas",
                                         "Roses France" = "roses_france",
                                         "Potatoes Texas" = "potatoes"),
                          selected = "roses_texas"),
      )
    )

  )
}

#' upload Server Functions
#'
#' @noRd
mod_upload_server <- function(input, output, session, data, parent_session){
  ns <- session$ns

  # Reset buttons
  values <- reactiveValues(
    upload_state_summary = 0,
    upload_state_baf = 0,
    upload_state_logr = 0,
  )

  observeEvent(input$reset_all, {
    values$upload_state_summary = 0
    values$upload_state_baf = 0
    values$upload_state_logr = 0
  })

  observeEvent(input$submit_summary, {
    values$upload_state_summary <- 'uploaded'
    values$upload_state_summary <- 0
  })

  observeEvent(input$submit_logr_baf, {
    values$upload_state_baf <- 'uploaded'
    values$upload_state_baf <- 0
  })

  input_summary <- reactive({
    if(values$upload_state_summary == 0){
      return(NULL)
    } else if (values$upload_state_summary == "uploaded"){
      validate(
        need(!is.null(input$load_summary), "Upload summary file before submit")
      )
      summary <- vroom(input$load_summary$datapath)

      return(summary)
    }
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
      if(input$example == "roses_texas"){
        baf <- vroom(system.file("baf_roses_texas.txt", package = "Qploidy"))
        logR <- vroom(system.file("logr_roses_texas.txt", package = "Qploidy"))
      } else if(input$example == "roses_france"){
        logR <- vroom(system.file("baf_roses_texas.txt", package = "Qploidy"))
        baf <- vroom(system.file("logr_roses_texas.txt", package = "Qploidy"))
      } else if(input$example == "potatoes") {
        logR <- vroom(system.file("baf_roses_texas.txt", package = "Qploidy"))
        baf <- vroom(system.file("logr_roses_texas.txt", package = "Qploidy"))
      }
      data <- list(logR, baf)
      data
    } else NULL
  })

  logR_BAF = reactive({
    if(is.null(input_logr_baf()))
      return(loadExample())
    else return(input_logr_baf())
  })

  return(list(summary = reactive(input_summary()),
              logR_BAF = reactive(logR_BAF())))
}

## To be copied in the UI
# mod_upload_ui("upload_1")

## To be copied in the server
# mod_upload_server("upload_1")
