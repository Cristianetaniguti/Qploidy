#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#'
#' @import shiny
#' @import shinythemes
#' @import shinyWidgets
#' @import shinydashboard
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#'
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic
    dashboardPage(
      dashboardHeader(disable = TRUE),
      dashboardSidebar(disable = TRUE),
      dashboardBody(
        # Lab colors
        tags$head(
          tags$style(
            HTML('
         a.action-button {
                color: #6c81c0;
         }
         a.action-button.red {
                color: #6c81c0;
         }
        .navbar-static-top {background-color: green;}

        .navbar-default .navbar-nav>.active>a {background-color: green;}

        .body {
            background-color: #22284c;
        }

        .box {margin-bottom: 40px;}

        .box.box-solid > .box-header > .box-tools .btn {
           position: relative;
           bottom: 5px;
           box-shadow: none;
        }

        .box.box-solid.box-primary > .box-header a,
        .box.box-solid.box-primary > .box-header .btn {
            color: #ffffff;
        }

        .box.box-solid.box-primary>.box-header {
          height: 50px;
          color:#fff;
          background:#6c81c0
        }

        .box.box-solid.box-primary{
          border-bottom-color:#6c81c0;
          border-left-color:#6c81c0;
          border-right-color:#6c81c0;
          border-top-color:#6c81c0;
        }

        .box.box-solid.box-info>.box-header {
          height: 50px;
          color:#fff;
          background:#22284c
        }

        .box.box-solid.box-info{
          border-bottom-color:#22284c;
          border-left-color:#22284c;
          border-right-color:#22284c;
          border-top-color:#22284c;
        }

        .box.box-solid.box-warning>.box-header {
          color:#fff;
          background:#a91021ff
        }

        .box.box-solid.box-warning{
          border-bottom-color:#a91021ff;
          border-left-color:#a91021ff;
          border-right-color:#a91021ff;
          border-top-color:#a91021ff;
        }
                              '))),
        tags$head(tags$style(HTML('.navbar-static-top {background-color: #22284c;}',
                                  '.navbar-default .navbar-nav>.active>a {background-color: #22284c;}'))),
        # Your application UI logic
        navbarPage("Qploidy",
                   theme = shinythemes::shinytheme("flatly"),
                   tabPanel("About", value = "about",
                            includeMarkdown(system.file("about.Rmd", package = "Qploidy")),
                   ),
                   tabPanel("Step 1 - data interpolation",
                            mod_interpolation_ui("interpolation_1")
                   ),
                   tabPanel("Step 2 - ploidy estimation",
                            mod_all_ui("all_1")
                            )
        )
      )
    )
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
