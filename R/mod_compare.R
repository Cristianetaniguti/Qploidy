#' compare UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_compare_ui <- function(id){
  ns <- NS(id)
  tagList(
 
  )
}
    
#' compare Server Functions
#'
#' @noRd 
mod_compare_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
 
  })
}
    
## To be copied in the UI
# mod_compare_ui("compare_1")
    
## To be copied in the server
# mod_compare_server("compare_1")
