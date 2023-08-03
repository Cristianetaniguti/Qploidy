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
      box(width= 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,  status="info", title = "Upload files", label = tags$b("Upload files"),
          #fileInput(ns("load_logR"), label = "Upload logR file"),
          fileInput(ns("load_baf"), label = "Upload BAF file")
      ),
      box(width= 12, solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,  status="info", title = "Choose examples", label = tags$b("Choose examples"),
          radioButtons(ns("example_data"), label = "Choose example data set",
                       choices = c("Example data" = "example_data",
                                   "Roses Texas" = "roses_texas",
                                   "Roses France" = "roses_france",
                                   "Potatoes Texas" = "potatoes"),
                       selected = "example_data")),
      box(width= 12, solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,  status="info", title = "Overall estimations options", label = tags$b("Overall estimations options"),
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
          actionButton(ns("run_overal"), "Run")), hr(),
      box(width= 12, solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,  status="info", title = "Linkage map", label = tags$b("Linkage map"),
          fileInput(ns("load_mappoly"), label = "Upload mappoly homoprob object"),
          fileInput(ns("load_polyorigin"), label = "Upload polyOrigin output"),
          numericInput(ns("ploidy_polyorigin"), label = "Maps ploidy", value = 4, step = 1),
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Tables - Overall analysis",
                 p("The table include the ploidy estimated by chromosome and overall using the area peak method.
            The measures of quality are the difference between first and second place in the area method and the correlation coefficient (Pearson) between estimated and expected peak position."),
                 box(width= 12, solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,  status="primary", title = "Ploidy estimations by chromosome table", label = tags$b("Ploidy estimations by chromosome table"),
                     downloadButton(ns('result.ploidy_download'), "Download"), br(), hr(),
                     DT::dataTableOutput(ns("result.ploidy")), br()
                 ),
                 box(width= 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,  status="primary", title = "Proportion of dots inside selected area", label = tags$b("Proportion of dots inside selected area"),
                     downloadButton(ns('dots.int_tot_mt_download'), "Download"), br(), hr(),
                     DT::dataTableOutput(ns("dots.int_tot_mt")), br()
                 ),
                 box(width= 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,  status="primary", title = "Difference between first and second place in area method table", label = tags$b("Difference between first and second place in area method table"),
                     downloadButton(ns('diff.first.second_download'), "Download"), br(), hr(),
                     DT::dataTableOutput(ns("diff.first.second")), br()
                 ),
                 box(width= 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,  status="primary", title = "Standard deviation inside area table", label = tags$b("Standard deviation inside area table"),
                     downloadButton(ns('sd_tot_mt_download'), "Download"), br(), hr(),
                     DT::dataTableOutput(ns("sd_tot_mt")), br()
                 ),
                 box(width= 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,  status="primary", title = "Highest correlation table", label = tags$b("Highest correlation table"),
                     downloadButton(ns('corr_tot_mt_download'), "Download"), br(), hr(),
                     DT::dataTableOutput(ns("corr_tot_mt")), br()
                 ),
                 box(width= 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,  status="primary", title = "Modes inside areas table", label = tags$b("Modes inside areas table"),
                     downloadButton(ns('modes_paste_tot_mt_download'), "Download"), br(), hr(),
                     DT::dataTableOutput(ns("modes_paste_tot_mt"))
                 )

        ),
        tabPanel("Graphics - Overall analysis",
                 p("The graphics include the ploidy estimated by chromosome and overall using the area peak method.
            The measures of quality are the difference between first and second place in the area method and the correlation coefficient (Pearson) between estimated and expected peak position."),
                 box(width= 12, solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,  status="primary", title = "Euploid individuals", label = tags$b("Euploid individuals ploidy"),
                     htmlOutput(ns("text_ploidy")), hr(),
                     plotOutput(ns("overall_p1"))),
                 box(width= 12, solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,  status="primary", title = "Aneuploidy individuals", label = tags$b("Aneuploidy individuals"),
                     htmlOutput(ns("aneuploid_text")), hr(),
                     plotOutput(ns("overall_p2")), br(),
                     box(width= 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,  status="primary", title = "Number chromosomes with each ploidy in aneuploidy individuals", label = tags$b("Number chromosomes with each ploidy in aneuploidy individuals"),
                         plotOutput(ns("overall_p3"))),
                     box(width= 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,  status="primary", title = "Aneuploid individuals table", label = tags$b("Aneuploid individuals table"),
                         downloadButton(ns('aneuploid_df_download'), "Download"), br(), br(),
                         DT::dataTableOutput(ns("aneuploid_df")))), br(), hr(),
                 box(width= 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,  status="primary", title = "Euploid weird individuals", label = tags$b("Euploid weird individuals"),
                     htmlOutput(ns("ploidy_all_weird"))),
                 box(width= 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,  status="primary", title = "Aneuploid weird individuals", label = tags$b("Aneuploid weird individuals"),
                     htmlOutput(ns("aneuploidy_all_weird")))

        ),
        tabPanel("Impact in mapping population",
                 box(width= 12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,  status="primary", title = "Table", label = tags$b("Table"),
                     p("MAPpoly"),
                     downloadButton(ns('breaks_mappoly_download'), "Download"), br(), br(),
                     DT::dataTableOutput(ns("breaks_mappoly_df")),
                     p("polyOrigin"),
                     downloadButton(ns('breaks_polyorigin_df_download'), "Download"), br(), br(),
                     DT::dataTableOutput(ns("breaks_polyorigin_df"))
                 ),
                 box(width= 12, solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,  status="primary", title = "Graphic", label = tags$b("Graphic"),
                     p("MAPpoly"),
                     plotOutput(ns("breaks_mappoly_plot")),
                     p("polyOrigin"),
                     plotOutput(ns("breaks_polyorigin_plot"))
                 )
        ),
        tabPanel("Graphics - Single individual analysis",
                 br(),
                 box(width= 12, solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,  status="info", title = "Single individual visualization options", label = tags$b("Single individual visualization options"),
                     numericInput(ns("area_single"), label = "Total area", value = 0.75, step = 0.1),
                     numericInput(ns("ploidy"), label = "Input ploidy", value = 2),
                     pickerInput(ns("graphics"),
                                 label = "Select sample for graphics",
                                 choices = "This will be updated after samples are selected above",
                                 selected = "This will be updated after samples are selected above",
                                 options = pickerOptions(
                                   size = 8,
                                   `selected-text-format` = "count > 3",
                                   `live-search`=TRUE,
                                   actionsBox = TRUE,
                                   dropupAuto = FALSE
                                 ),
                                 multiple = FALSE),
                     checkboxInput(ns("colors"), label = "Color area", value = TRUE),
                     checkboxInput(ns("add_lines"), label = "Add expected and estimated modes", value = TRUE),
                     numericInput(ns("dot.size"), label = "Dot size", value = 1),br(),
                     actionButton(ns("run_individual"), "Run")
                 ),
                 # box(width = 12, solidHeader = TRUE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = "Segmented logR plot",
                 #     column(12,
                 #            plotOutput(ns("plot_logR")), br(),
                 #     )
                 # ),
                 box(width = 12, solidHeader = TRUE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = "BAF plot",
                     column(12,
                            br(),
                            plotOutput(ns("plot_lines")), br(),
                            plotOutput(ns("plot_hist"))
                     ),
                 ),br(),
                 box(width = 12, solidHeader = TRUE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = "Haplotypes",
                     box(width = 12, solidHeader = FALSE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = "MAPpoly probabilities",
                         plotOutput(ns("plot_haplo_mappoly"), height = "1100px"),
                     ),
                     box(width = 12, solidHeader = FALSE, collapsible = TRUE,  collapsed = TRUE, status="primary", title = "PolyOrigin probabilities",
                         plotOutput(ns("plot_haplo_polyorigin"), height = "1100px"))
                 )
        )
      )
    )
  )
}

#' all Server Functions
#'
#' @importFrom stats coef cor lm median sd
#' @importFrom utils data write.csv write.table
#' @import dpseg
#'
#' @noRd
mod_all_server <- function(id){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    input_logR_baf <- reactive({
      withProgress(message = 'Working:', value = 0, {
        incProgress(0.5, detail = paste("Loading BAF file..."))

        #if(!is.null(input$load_logR) & !is.null(input$load_baf)){
        if(!is.null(input$load_baf)){
          #logR <- vroom(input$load_logR$datapath, show_col_types = FALSE)
          logR <- NULL
          baf <- vroom(input$load_baf$datapath, show_col_types = FALSE)
          if(!is.null(input$load_mappoly)){
            haplo_mappoly <- load(input$load_mappoly$datapath)
          } else haplo_mappoly <- NULL

          if(!is.null(input$load_polyorigin)){ #fixme
            genofile <- read_polyOrigin(input$load_polyorigin$datapath)
            f1.codes <- vroom(system.file("F1codes.polyorigin.txt", package = "Qploidy"), show_col_types = FALSE)
            haplo_polyorigin <- get_probs_polyorigin(genofile,
                                                     f1.codes = f1.codes,
                                                     ploidy = input$ploidy_polyorigin, n.cores = 2)
          } else haplo_polyorigin <- NULL
          return(list(logR, baf, haplo_mappoly, haplo_polyorigin))
        } else {
          return(NULL)
        }
      })
    })

    loadExample <- reactive({
      withProgress(message = 'Working:', value = 0, {
        incProgress(0.3, detail = paste("Loading BAF file..."))

        if(is.null(input_logR_baf())){
          if(input$example_data == "example_data"){
            #logR <- vroom(system.file("logR.example.txt", package = "Qploidy"), show_col_types = FALSE)
            logR <- NULL
            baf <- vroom(system.file("baf.example.txt", package = "Qploidy"), show_col_types = FALSE)
            temp <- load(system.file("mappoly.homoprob.ex.RData", package = "Qploidy"))
            haplo_mappoly <- get(temp)
            polyorigin  <- vroom(system.file("genofile_sub.csv", package = "Qploidy"), show_col_types = FALSE)
            f1.codes <- vroom(system.file("F1codes.polyorigin.txt", package = "Qploidy"), show_col_types = FALSE)
            ploidy <- 4
            haplo_polyorigin <- get_probs_polyorigin(polyorigin,
                                                     f1.codes = f1.codes,
                                                     ploidy = 4, n.cores = 2)
            return(list(logR, baf, haplo_mappoly, haplo_polyorigin))
          } else if(input$example_data == "roses_texas"){
            #logR <- vroom("C:/Users/Rose_Lab/Documents/Cris_temp/Qploidy_data/roses_texas/fitpoly/logR_roses_texas.txt", show_col_types = FALSE)
            logR <- NULL
            baf <- vroom("C:/Users/Rose_Lab/Documents/Cris_temp/Qploidy_data/roses_texas/fitpoly/baf_roses_texas.txt", show_col_types = FALSE)
            incProgress(0.3, detail = paste("Loading MAPpoly file..."))

            haplo_mappoly <- readRDS("C:/Users/Rose_Lab/Documents/Cris_temp/TAMU-RoseLab/standalone_apps/ploidy_estimation/data/count_breaks_poly/homoprob_normal.RDS")
            # polyorigin  <- read_polyOrigin("C:/Users/Rose_Lab/Documents/Cris_temp/Qploidy_data/roses_texas/5pop_compareprob2_fulldata_polyancestry.csv")
            # f1.codes <- vroom(system.file("F1codes.polyorigin.txt", package = "Qploidy"), show_col_types = FALSE)
            # ploidy <- 4
            # haplo_polyorigin <- Qploidy:::get_probs_polyorigin_sd(polyorigin,
            #                                          f1.codes = f1.codes,
            #                                          ploidy = 4, n.cores = 20)
            # saveRDS(haplo_polyorigin, file = "homoprob_polyorigin.rds")
            incProgress(0.3, detail = paste("Loading PolyOrigin file..."))
            haplo_polyorigin <- readRDS("C:/Users/Rose_Lab/Documents/Cris_temp/Qploidy_data/roses_texas/homoprob_polyorigin.rds")
            return(list(logR, baf, haplo_mappoly, haplo_polyorigin))
          } else if(input$example_data == "roses_france"){
            cat("Developing")
          } else if(input$example_data == "potatoes") {
            cat("Developing")
          }
        } else {
          return(NULL)
        }
      })
    })

    logR_baf <- reactive({
      if(is.null(input_logR_baf())){
        logR_baf <- loadExample()[1:2]
      } else {
        logR_baf <- input_logR_baf()[1:2]
      }
      return(logR_baf)
    })

    haplo <- reactive({
      if(is.null(input_logR_baf()$haplo_mappoly) & is.null(input_logR_baf()$haplo_polyorigin)){
        haplo <- loadExample()[3:4]
      } else {
        haplo <- input_logR_baf()[3:4]
      }
      return(haplo)
    })

    observe({
      choices_names <- as.list(unique(colnames(logR_baf()[[2]])[-c(1:3)]))
      names(choices_names) <- unique(colnames(logR_baf()[[2]])[-c(1:3)])

      updatePickerInput(session, "samples",
                        label = "Select samples for overall analysis",
                        choices = choices_names,
                        selected=unlist(choices_names)[1])
    })

    est.ploidy.chr_df <- eventReactive(input$run_overal,{
      withProgress(message = 'Working:', value = 0, {
        incProgress(0.5, detail = paste("Estimating ploidy..."))

        data_sample <- logR_baf()[[2]][,c(2,3,which(colnames(logR_baf()[[2]]) %in% c(input$samples)))]

        data_sample <- data_sample[order(data_sample$Chr, data_sample$Position),]

        est.ploidy.chr_df <- area_estimate_ploidy_by_chr(data_sample, ploidys = input$ploidys, area = input$area)
        est.ploidy.chr_df
      })
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
      withProgress(message = 'Working:', value = 0, {
        incProgress(0.5, detail = paste("Generating plots..."))

        # Build graphics to overview the estimations
        ps <- plots_overall(est.ploidy.chr_df(),  input$filter_diff, input$filter_corr)
        return(ps)
      })
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

    graphics_breaks <- reactive({
      withProgress(message = 'Working:', value = 0, {
        incProgress(0.5, detail = paste("Counting breaks..."))

        aneuploids <- data.frame(X = rownames(est.ploidy.chr_df()[[1]]), est.ploidy.chr_df()[[1]])

        if(any(unique(haplo()[[1]]$homoprob$individual) %in% input$samples)){
          breaks_mappoly_df <- count_breaks_df(homoprob =haplo()[[1]]$homoprob,
                                               aneuploids = aneuploids,
                                               inds = input$samples)

          breaks_mappoly_plot <- count_breaks_plot(breaks_mappoly_df,
                                                   by_LG = FALSE)
          breaks_mappoly <- list(df=breaks_mappoly_df, plot = breaks_mappoly_plot)
        } else breaks_mappoly <- NULL

        if(any(unique(haplo()[[2]]$homoprob$individual) %in% input$samples)){
          breaks_polyorigin_df <- count_breaks_df(homoprob =haplo()[[2]]$homoprob,
                                                  aneuploids = aneuploids,
                                                  inds = input$samples)

          breaks_polyorigin_plot <- count_breaks_plot(breaks_polyorigin_df,
                                                      by_LG = FALSE)

          breaks_polyorigin <- list(df = breaks_polyorigin_df, plot = breaks_polyorigin_plot)
        } else breaks_polyorigin <- NULL

        return(list(breaks_mappoly, breaks_polyorigin))
      })
    })

    output$breaks_mappoly_plot <- renderPlot({
      graphics_breaks()[[1]]$plot
    })

    output$breaks_polyorigin_plot <- renderPlot({
      graphics_breaks()[[2]]$plot
    })

    output$breaks_mappoly_df <- DT::renderDataTable({
      DT::datatable(graphics_breaks()[[1]]$df, extensions = 'Buttons',
                    options = list(
                      scrollX = TRUE,
                      dom = 'Bfrtlp',
                      buttons = c('copy', 'csv', 'excel', 'pdf')
                    ),
                    class = "display")
    })

    output$breaks_mappoly_df_download <- downloadHandler(
      filename = "table.csv",
      content = function(file) {
        write.csv(graphics_breaks()[[1]]$df, file = file)
      }
    )

    output$breaks_polyorigin_df <- DT::renderDataTable({
      DT::datatable(graphics_breaks()[[2]]$df, extensions = 'Buttons',
                    options = list(
                      scrollX = TRUE,
                      dom = 'Bfrtlp',
                      buttons = c('copy', 'csv', 'excel', 'pdf')
                    ),
                    class = "display")
    })

    output$breaks_polyorigin_df_download <- downloadHandler(
      filename = "table.csv",
      content = function(file) {
        write.csv(graphics_breaks()[[2]]$df, file = file)
      }
    )

    observe({
      choices_names <- as.list(unique(colnames(logR_baf()[[2]])[-c(1:3)]))
      names(choices_names) <- unique(colnames(logR_baf()[[2]])[-c(1:3)])

      updatePickerInput(session, "graphics",
                        label = "Select sample for the graphic",
                        choices = choices_names,
                        selected=unlist(choices_names)[1])
    })

    # Single individual analysis
    graphics_baf <- eventReactive(input$run_individual,{
      withProgress(message = 'Working:', value = 0, {
        incProgress(0.5, detail = paste("Generating individual BAF plots..."))

        data_sample <- logR_baf()[[2]][,c(2,3,which(colnames(logR_baf()[[2]]) %in% c(input$graphics)))]
        colnames(data_sample)[3] <- "sample"

        p_baf <- plot_baf(data_sample, input$area_single, input$ploidy, input$dot.size, input$add_lines, input$colors)
        p_hist <- plot_baf_hist(data_sample, input$area_single, input$ploidy, input$colors, input$add_lines)

        return(list(p_baf, p_hist))
      })
    })

    individual_haplo <- reactive({
      withProgress(message = 'Working:', value = 0, {
        incProgress(0.5, detail = paste("Generating individual haplotypes plots..."))

        mappoly <- haplo()[[1]]
        if(all(grep("^X", unique(mappoly$homoprob$individual))))
          mappoly$homoprob$individual <- gsub("^X", "", mappoly$homoprob$individual)

        print(unique(mappoly$homoprob$individual))
        if(any(unique(mappoly$homoprob$individual) %in% input$graphics)){
          haplo_lst <- list()
          for(i in 1:length(unique(mappoly$homoprob$homoprob$LG))){
            haplo_lst[[i]] <- plot(mappoly, lg = unique(mappoly$homoprob$LG)[i],
                                   ind = input$graphics,
                                   use.plotly = FALSE)
          }

          all_haplo_mappoly <- ggarrange(plotlist = haplo_lst, common.legend = TRUE)
        } else all_haplo_mappoly <- NULL

        polyorigin <- haplo()[[2]]
        if(all(grep("^X", unique(polyorigin$homoprob$individual))))
          polyorigin$homoprob$individual <- gsub("^X", "",polyorigin$homoprob$individual)

        print(unique(polyorigin$homoprob$individual))

        if(any(polyorigin$homoprob$individual %in% input$graphics)){
          haplo_lst <- list()
          for(i in 1:length(unique(polyorigin$homoprob$LG))){
            haplo_lst[[i]] <- plot(polyorigin, lg = unique(polyorigin$homoprob$LG)[i],
                                   ind = input$graphics,
                                   use.plotly = FALSE)
          }

          all_haplo_polyorigin <- ggarrange(plotlist = haplo_lst, common.legend = TRUE)
        } else all_haplo_polyorigin <- NULL

        return(list(all_haplo_mappoly, all_haplo_polyorigin))
      })
    })

    output$plot_lines <- renderPlot({
      graphics_baf()[[1]]
    })

    output$plot_hist <- renderPlot({
      graphics_baf()[[2]]
    })

    output$plot_haplo_mappoly <- renderPlot({
      validate(
        need(!is.null(individual_haplo()[[1]]), "MAPpoly haplotype information not provided for this individual."),
      )
      individual_haplo()[[1]]
    })

    output$plot_haplo_polyorigin <- renderPlot({
      validate(
        need(!is.null(individual_haplo()[[2]]), "PolyOrigin haplotype information not provided for this individual."),
      )
      individual_haplo()[[2]]
    })
  })
}

## To be copied in the UI
# mod_all_ui("all_1")

## To be copied in the server
# mod_all_server("all_1")
