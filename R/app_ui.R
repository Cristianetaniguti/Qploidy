#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinythemes
#' @import shinyWidgets
#' @import shinydashboard
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @import reshape2
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic
    navbarPage("Qploidy",
               theme = shinythemes::shinytheme("flatly"),
               tabPanel("About", value = "about",
                        includeMarkdown(system.file("about.Rmd", package = "Qploidy")),
               ),
               tabPanel("Upload data",
                        mod_upload_ui("upload_1")),
               navbarMenu("Data interpolation",
                          tabPanel("Reference", "Summary tab contents..."),
                          tabPanel("Graphics", "Summary tab contents...")
               ),
               tabPanel("All individuals",
                        mod_all_ui("all_1")
               ),
               tabPanel("Single individual",
                        mod_single_ui("single_1")
               )

    ),

  )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "Qploidy"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
