#' all UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_all_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      box(width= 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,  status="primary", title = "Upload files", label = tags$b("Upload files"),
          fileInput(ns("load_logR"), label = "Upload logR file"),
          fileInput(ns("load_baf"), label = "Upload BAF file")
      ),
      box(width= 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,  status="primary", title = "Choose examples", label = tags$b("Choose examples"),
          radioButtons(ns("example_data"), label = "Choose example data set",
                       choices = c("Roses Texas" = "roses_texas",
                                   "Roses France" = "roses_france",
                                   "Potatoes Texas" = "potatoes"),
                       selected = "roses_texas")),
      sliderInput(ns("ploidys"), label = "select ploidy", min = 2, max = 8, value = c(2,5), step = 1),
      numericInput(ns("area"), label = "Total area", value = 0.75, step = 0.1),
      numericInput(ns("filter_diff"), label = "Minimum filter difference", value = 0, step = 0.01),
      numericInput(ns("filter_corr"), label = "Minimum correlation between estimated and expected peaks", value = 0, step = 0.01),
      pickerInput(ns("samples"),
                  label = "Select samples for overall analysis",
                  choices = "This will be updated with files in data/joint_logR_BAF",
                  selected = "This will be updated with files in data/joint_logR_BAF",
                  options = pickerOptions(
                    size = 8,
                    `selected-text-format` = "count > 3",
                    `live-search`=TRUE,
                    actionsBox = TRUE,
                    dropupAuto = FALSE
                  ),
                  multiple = TRUE),
      actionButton(ns("run_overal"), "Run")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Tables",
                 p("The table include the ploidy estimated by chromosome and overall using the area peak method.
            The measures of quality are the difference between first and second place in the area method and the correlation coefficient (Pearson) between estimated and expected peak position."),
                 column(12,
                        box(width= 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,  status="primary", title = "Ploidy estimations by chromosome table", label = tags$b("Ploidy estimations by chromosome table"),
                            downloadButton(ns('result.ploidy_download'), "Download"), br(), hr(),
                            DT::dataTableOutput(ns("result.ploidy")), br()
                        ),
                        box(width= 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,  status="primary", title = "Proportion of dots inside selected area", label = tags$b("Proportion of dots inside selected area"),
                            downloadButton(ns('dots.int_tot_mt_download'), "Download"), br(), hr(),
                            DT::dataTableOutput(ns("dots.int_tot_mt")), br()
                        ),
                        box(width= 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,  status="primary", title = "Difference between first and second place in area method table", label = tags$b("Difference between first and second place in area method table"),
                            downloadButton(ns('diff.first.second_download'), "Download"), br(), hr(),
                            DT::dataTableOutput(ns("diff.first.second")), br()
                        ),
                        box(width= 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,  status="primary", title = "Standard deviation inside area table", label = tags$b("Standard deviation inside area table"),
                            downloadButton(ns('sd_tot_mt_download'), "Download"), br(), hr(),
                            DT::dataTableOutput(ns("sd_tot_mt")), br()
                        ),
                        box(width= 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,  status="primary", title = "Highest correlation table", label = tags$b("Highest correlation table"),
                            downloadButton(ns('corr_tot_mt_download'), "Download"), br(), hr(),
                            DT::dataTableOutput(ns("corr_tot_mt")), br()
                        ),
                        box(width= 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,  status="primary", title = "Modes inside areas table", label = tags$b("Modes inside areas table"),
                            downloadButton(ns('modes_paste_tot_mt_download'), "Download"), br(), hr(),
                            DT::dataTableOutput(ns("modes_paste_tot_mt")), br()
                        )
                 )
        ),
        tabPanel("Graphics",
                 p("The graphics include the ploidy estimated by chromosome and overall using the area peak method.
            The measures of quality are the difference between first and second place in the area method and the correlation coefficient (Pearson) between estimated and expected peak position."),
                 column(12,
                        p("Euploid individuals ploidy:"),
                        htmlOutput(ns("text_ploidy")), br(),
                        p("Number of euploid individuals:"),
                        plotOutput(ns("overall_p1")), br(),
                        p("Aneuploid individuals:"),
                        p("Ploidy of each chromosome in aneuploidy individuals:"),
                        plotOutput(ns("overall_p2")), br(),
                        p("Number chromosomes with each ploidy in aneuploidy individuals:"),
                        plotOutput(ns("overall_p3")),
                        p("Aneuploid individuals:"),
                        htmlOutput(ns("aneuploid_text")), br(), hr(),
                        downloadButton(ns('aneuploid_df_download'), "Download"), br(),
                        DT::dataTableOutput(ns("aneuploid_df")), br(), hr(),
                        p("Euploid weird individuals:"),
                        htmlOutput(ns("ploidy_all_weird")), br(),
                        p("Aneuploid weird individuals:"),
                        htmlOutput(ns("aneuploidy_all_weird")), br(),
                 )
        )
      )
    )
  )
}

#' all Server Functions
#'
#' @noRd
mod_all_server <- function(id){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    input_logR_baf <- reactive({
      if(!is.null(input$load_logR) & !is.null(input$load_baf)){
        logR <- vroom(input$load_logR$datapath)
        baf <- vroom(input$load_baf$datapath)
        return(list(logR, baf))
      } else {
        return(NULL)
      }
    })

    loadExample <- reactive({
      if(is.null(input_logR_baf())){
        if(input$example_data == "roses_texas"){
          logR <- vroom(system.file("logR.example.txt", package = "Qploidy"))
          baf <- vroom(system.file("baf.example.txt", package = "Qploidy"))
          return(list(logR, baf))
        } else if(input$example_data == "roses_france"){
          cat("Developing")
        } else if(input$example_data == "potatoes") {
          cat("Developing")
        }
      } else {
        return(NULL)
      }
    })

    logR_baf <- reactive({
      if(is.null(input_logR_baf())){
        logR_baf <- loadExample()
      } else {
        logR_baf <- input_logR_baf()
      }
      return(logR_baf)
    })

    observe({
      choices_names <- as.list(unique(colnames(logR_baf()[[1]])[-c(1:3)]))
      names(choices_names) <- unique(colnames(logR_baf()[[1]])[-c(1:3)])

      updatePickerInput(session, "samples",
                        label = "Select samples for overall analysis",
                        choices = choices_names,
                        selected=unlist(choices_names)[1])
    })

    est.ploidy.chr_df <- eventReactive(input$run_overal,{
      data_sample <- logR_baf()[[2]][,c(2,3,which(colnames(logR_baf()[[2]]) %in% c(input$samples)))]

      data_sample <- data_sample[order(data_sample$Chr, data_sample$Position),]

      est.ploidy.chr_df <- area_estimate_ploidy_by_chr(data_sample, ploidy = input$ploidys, area = input$area)
      est.ploidy.chr_df
    })

    output$result.ploidy <- DT::renderDataTable({
      DT::datatable(est.ploidy.chr_df()[[1]], extensions = 'Buttons',
                    options = list(
                      scrollX = TRUE,
                      dom = 'Bfrtlp',
                      buttons = c('copy', 'csv', 'excel', 'pdf')
                    ),
                    class = "display")
    })

    output$result.ploidy_download <- downloadHandler(
      filename = "table.csv",
      content = function(file) {
        write.csv(est.ploidy.chr_df()[[1]], file = file)
      }
    )

    output$dots.int_tot_mt <- DT::renderDataTable({
      DT::datatable(est.ploidy.chr_df()[[2]], extensions = 'Buttons',
                    options = list(
                      scrollX = TRUE,
                      dom = 'Bfrtlp',
                      buttons = c('copy', 'csv', 'excel', 'pdf')
                    ),
                    class = "display")
    })

    output$dots.int_tot_mt_download <- downloadHandler(
      filename = "table.csv",
      content = function(file) {
        write.csv(est.ploidy.chr_df()[[2]], file = file)
      }
    )

    output$diff.first.second <- DT::renderDataTable({
      DT::datatable(est.ploidy.chr_df()[[3]], extensions = 'Buttons',
                    options = list(
                      scrollX = TRUE,
                      dom = 'Bfrtlp',
                      buttons = c('copy', 'csv', 'excel', 'pdf')
                    ),
                    class = "display")
    })

    output$diff.first.second_download <- downloadHandler(
      filename = "table.csv",
      content = function(file) {
        write.csv(est.ploidy.chr_df()[[3]], file = file)
      }
    )

    output$sd_tot_mt <- DT::renderDataTable({
      DT::datatable(est.ploidy.chr_df()[[4]], extensions = 'Buttons',
                    options = list(
                      scrollX = TRUE,
                      dom = 'Bfrtlp',
                      buttons = c('copy', 'csv', 'excel', 'pdf')
                    ),
                    class = "display")
    })

    output$sd_tot_mt_download <- downloadHandler(
      filename = "table.csv",
      content = function(file) {
        write.csv(est.ploidy.chr_df()[[4]], file = file)
      }
    )

    output$corr_tot_mt <- DT::renderDataTable({
      DT::datatable(est.ploidy.chr_df()[[5]], extensions = 'Buttons',
                    options = list(
                      scrollX = TRUE,
                      dom = 'Bfrtlp',
                      buttons = c('copy', 'csv', 'excel', 'pdf')
                    ),
                    class = "display")
    })

    output$corr_tot_mt_download <- downloadHandler(
      filename = "table.csv",
      content = function(file) {
        write.csv(est.ploidy.chr_df()[[5]], file = file)
      }
    )

    output$modes_paste_tot_mt <- DT::renderDataTable({
      DT::datatable(est.ploidy.chr_df()[[6]], extensions = 'Buttons',
                    options = list(
                      scrollX = TRUE,
                      dom = 'Bfrtlp',
                      buttons = c('copy', 'csv', 'excel', 'pdf')
                    ),
                    class = "display")
    })

    output$modes_paste_tot_mt_download <- downloadHandler(
      filename = "table.csv",
      content = function(file) {
        write.csv(est.ploidy.chr_df()[[6]], file = file)
      }
    )

    est.ploidy.chr_plots <- reactive({
      # Build graphics to overview the estimations
      ps <- plots_overall(est.ploidy.chr_df,  input$filter_diff, input$filter_corr)
      ps
    })

    output$overall_p1 <- renderPlot({
      est.ploidy.chr_plots()[[1]]
    })

    output$text_ploidy <- renderUI({
      HTML(est.ploidy.chr_plots()[[4]])
    })

    output$overall_p2 <- renderPlot({
      est.ploidy.chr_plots()[[2]]
    })

    output$overall_p3 <- renderPlot({
      est.ploidy.chr_plots()[[3]]
    })

    output$aneuploid_text <- renderUI({
      HTML(est.ploidy.chr_plots()[[5]])
    })

    output$ploidy_all_weird <- renderUI({
      HTML(est.ploidy.chr_plots()[[7]])
    })

    output$aneuploidy_all_weird <- renderUI({
      HTML(est.ploidy.chr_plots()[[8]])
    })

    output$aneuploid_df <- DT::renderDataTable({
      DT::datatable(est.ploidy.chr_plots()[[6]], extensions = 'Buttons',
                    options = list(
                      scrollX = TRUE,
                      dom = 'Bfrtlp',
                      buttons = c('copy', 'csv', 'excel', 'pdf')
                    ),
                    class = "display")
    })

    output$aneuploid_df_download <- downloadHandler(
      filename = "table.csv",
      content = function(file) {
        write.csv(est.ploidy.chr_plots()[[6]], file = file)
      }
    )
  })
}

## To be copied in the UI
# mod_all_ui("all_1")

## To be copied in the server
# mod_all_server("all_1")
