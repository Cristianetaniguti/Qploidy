#' interpolation UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_interpolation_ui <- function(id){
  ns <- NS(id)
  tagList(
 
  )
}
    
#' interpolation Server Functions
#'
#' @noRd 
mod_interpolation_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
 
  })
}
    
## To be copied in the UI
# mod_interpolation_ui("interpolation_1")
    
## To be copied in the server
# mod_interpolation_server("interpolation_1")
