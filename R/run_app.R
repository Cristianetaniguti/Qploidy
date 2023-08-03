#' Run the Shiny Application
#'
#' @param ... arguments to pass to golem_opts. See `?get_golem_options` for more details.
#' @param onStart to be documented
#' @param options to be documented
#' @param enableBookmarking to be documented
#' @param uiPattern to be documented
#'
#' @inheritParams shinyApp
#'
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
#'
#' @export
run_app <- function(..., onStart = NULL,  options = list(),  enableBookmarking = NULL,  uiPattern = "/") {
  with_golem_options(
    app = shinyApp(
      ui = app_ui,
      server = app_server,
      onStart = onStart,
      options = options,
      enableBookmarking = enableBookmarking,
      uiPattern = uiPattern
    ),
    golem_opts = list(...)
  )
}
