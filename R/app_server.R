#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  # Your application server logic

  mod_interpolation_server("interpolation_1")

  # callModule(mod_all_server,
  #            "all_ui_1",
  #            parent_session = session)
  #
  # callModule(mod_single_server,
  #            "single_ui_1",
  #            parent_session = session)
}
