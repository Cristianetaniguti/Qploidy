#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  # Your application server logic

  datas <- callModule(mod_upload_server,
                      "upload_ui_1",
                      parent_session=session)

  callModule(mod_all_server,
             "all_ui_1",
             data = datas$data,
             parent_session = session)

  callModule(mod_single_server,
             "single_ui_1",
             data = datas$data,
             parent_session = session)
}
