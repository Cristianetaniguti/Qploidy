#' help UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList includeMarkdown
mod_help_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      column(width=12),
      column(width=12,
             box(title="Ploidy Estimation", id = "Qploidy_box",width = 12, collapsible = TRUE, collapsed = TRUE, status = "info", solidHeader = TRUE,
                 "This tab uses Qploidy package to standardize marker allele counts to estimate ploidy or aneuploidy.",
                 br(), br(),
                 bs4Dash::tabsetPanel(id = "Qploidy_tabset",
                                      tabPanel("Parameters description", value = "Qploidy_par", br(),
                                               includeMarkdown(system.file("help_files/Qploidy_par.Rmd", package = "Qploidy"))
                                      ),
                                      tabPanel("Results description", value = "Qploidy_results", br(),
                                               includeMarkdown(system.file("help_files/Qploidy_res.Rmd", package = "Qploidy"))
                                      ),
                                      tabPanel("How to cite", value = "Qploidy_cite", br(),
                                               includeMarkdown(system.file("help_files/Qploidy_cite.Rmd", package = "Qploidy"))
                                      ))
             )
      ),
      column(width=2)
      # Add Help content here
    )
  )
}

#' help Server Functions
#'
#' @noRd
mod_help_server <- function(input, output, session, parent_session){

  ns <- session$ns

}

## To be copied in the UI
# mod_help_ui("help_1")

## To be copied in the server
# mod_help_server("help_1")
