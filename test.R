library(shiny)
library(shinydashboard)

ui <- dashboardPage(
  dashboardHeader(), dashboardSidebar(),
  dashboardBody(valueBoxOutput("qc_box"))
)

server <- function(input, output, session){
  total  <- 1234
  passed <- 987
  failed <- total - passed

  output$qc_box <- renderValueBox({
    valueBox(
      value = tagList(
        div(style="font-size:1.6em; line-height:1.1;", formatC(total, big.mark=",")),
        div(style="font-size:0.95em; opacity:0.9;",
            span("Passed: "), strong(passed), " • ",
            span("Failed: "), strong(failed))
      ),
      subtitle = "QC summary",
      icon = icon("flask"),
      color = if (failed == 0) "green" else "orange",
      width = 4
    )
  })
}

shinyApp(ui, server)
