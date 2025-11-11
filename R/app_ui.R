#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom bs4Dash bs4Badge bs4DashSidebar bs4DashNavbar bs4DashPage sidebarMenu menuItem menuSubItem dashboardBody tabItems tabItem box dashboardFooter
#' @import shinyWidgets
#'
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic
    bs4DashPage(
      skin = "black",
      bs4DashNavbar(
        title = tagList(
          tags$img(src = 'www/Qploidy_logo.png', height = '40', width = '35')
        ),
        rightUi = tags$li(
          class = "dropdown",
          tags$a(
            href = "#",
            class = "nav-link",
            `data-toggle` = "dropdown",
            icon("info-circle")
          ),
          tags$div(
            class = "dropdown-menu dropdown-menu-right",
            tags$a(
              class = "dropdown-item",
              href = "#",
              "Session Info",
              onclick = "Shiny.setInputValue('session_info_button', Math.random())"
            )
          )
        )
      ),
      help = NULL, #This is the default bs4Dash button to control the presence of tooltips and popovers, which can be added as a user help/info feature.
      bs4DashSidebar(
        skin="light",
        status = "info",
        fixed=TRUE,
        #minified = F,
        expandOnHover = TRUE,
        sidebarMenu(id = "MainMenu",
                    flat = FALSE,
                    tags$li(class = "header", style = "color: grey; margin-top: 10px; margin-bottom: 10px; padding-left: 15px;", "Menu"),
                    menuItem("Home", tabName = "welcome", icon = icon("house"),startExpanded = FALSE),
                    menuItem(
                      span("Ploidy Estimation"),
                      tabName = "qploidy",
                      icon = icon("dna")),
                    tags$li(class = "header", style = "color: grey; margin-top: 18px; margin-bottom: 10px; padding-left: 15px;", "Information"),
                    menuItem("Source Code", icon = icon("circle-info"), href = "https://www.github.com/Cristianetaniguti/Qploidy"),
                    menuItem("Help", tabName = "help", icon = icon("circle-question"))
        )
      ),
      footer = dashboardFooter(
        right = div(
          style = "display: flex; align-items: center;",  # Align text and images horizontally
          div(
            style = "display: flex; flex-direction: column; margin-right: 15px; text-align: right;",
            div("2025 Breeding Insight"),
            div("Funded by USDA through Cornell University")
          ),
          div(
            a(
              img(src = "www/usda-logo-color.png", height = "45px"),
              style = "margin-right: 15px;"
            ),
            a(
              img(src = "www/cornell_seal_simple_web_b31b1b.png", height = "45px")
            )
          )
        ),
        left = div(
          style = "display: flex; align-items: center; height: 100%;",
          sprintf("v%s", as.character(utils::packageVersion("Qploidy")))
        )
      ),
      dashboardBody(
        disconnectMessage(), #Adds generic error message for any error if not already accounted for
        tags$style(
          HTML(
            ".main-footer {
            background-color: white;
            color: grey;
            height: 65px;
            padding-top: 5px;
            padding-bottom: 5px;
          }
          .main-footer a {
            color: grey;
          }"
          )
        ),
        tabItems(
          tabItem(
            tabName = "welcome", mod_Home_ui("Home_1")
          ),
          tabItem(
            tabName = "qploidy", mod_qploidy_ui("qploidy_1")
          ),
          tabItem(
            tabName = "help", mod_help_ui("help_1")
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
