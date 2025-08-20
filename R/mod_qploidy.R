#' qploidy UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @importFrom shinyWidgets virtualSelectInput
#' @import shinydisconnect
#' @importFrom future availableCores
#' @importFrom DT DTOutput
#' @importFrom bs4Dash valueBoxOutput
#'
mod_qploidy_ui <- function(id){
  ns <- NS(id)
  tagList(
    # Add qploidy content here
    fluidRow(
      disconnectMessage(
        text = "An unexpected error occurred, please reload the application and check the input file(s).",
        refresh = "Reload now",
        background = "white",
        colour = "grey",
        overlayColour = "grey",
        overlayOpacity = 0.3,
        refreshColour = "purple"
      ),
      column(width = 3,
             box(title="Inputs", width = 12, collapsible = TRUE, collapsed = FALSE, status = "info", solidHeader = TRUE,
                 "* Required",
                 selectInput(ns('file_type'),
                             label = 'Select File Format',
                             choices = c("VCF","Axiom Array", "Illumina Array", "Qploidy Standardized Dataset"),
                             selected = "VCF"),
                 fileInput(ns("input_file"), "Choose File*", accept = c(".csv",".vcf",".gz")),
                 conditionalPanel(condition = "input.file_type == 'Axiom Array' || input.file_type == 'Illumina Array'",
                                  ns = ns,
                                  fileInput(ns("marker_pos"), "Input Marker Position file")
                 ),
                 conditionalPanel(condition = "input.file_type != 'Qploidy Standardized Dataset'",
                                  ns=ns,
                                  radioButtons(ns("known_ploidy"),
                                               label = "Do you have ≥60 samples with the same, known ploidy?",
                                               choices = list("Yes"= "TRUE", "No" = "FALSE"),
                                               selected = "TRUE"),
                                  conditionalPanel(condition = "input.file_type == 'VCF' && input.known_ploidy == 'TRUE'",
                                                   ns=ns,
                                                   fileInput(ns("reference_samples"), "Upload a CSV with a single Sample_ID column of reference IDs*")),
                                  conditionalPanel(condition = "(input.file_type == 'Axiom Array' || input.file_type == 'Illumina Array') && input.known_ploidy == 'TRUE'",
                                                   ns=ns,
                                                   fileInput(ns("geno_ref"), "Provide fitpoly result _scores file for the reference samples*")),
                                  conditionalPanel(condition = "input.known_ploidy == 'FALSE' && (input.file_type == 'Axiom Array' || input.file_type == 'Illumina Array')",
                                                   ns = ns,
                                                   fileInput(ns("geno_all"), "Provide fitpoly result _scores file for all samples*")
                                  ),
                                  sliderInput(ns("cores"), "Number of CPU Cores*", min = 1, max = (availableCores() - 1), value = 1, step = 1),
                                  conditionalPanel(condition = "input.known_ploidy == 'TRUE'",
                                                   ns=ns,
                                                   numericInput(ns("ref_ploidy"), "Reference Samples Ploidy*", min = 1, value = NULL)),
                                  conditionalPanel(condition = "input.known_ploidy == 'FALSE'",
                                                   ns=ns,
                                                   numericInput(ns("ref_ploidy"), "Most Common Ploidy*", min = 1, value = NULL)),
                                  div(style="text-align: left; margin-top: 10px;",
                                      actionButton(ns("advanced_options"),
                                                   label = HTML(paste(icon("cog", style = "color: #007bff;"), "Advanced Options")),
                                                   style = "background-color: transparent; border: none; color: #007bff; font-size: smaller; text-decoration: underline; padding: 0;"
                                      )
                                  ),
                 ),
                 actionButton(ns("run_stand"), "Run Analysis"),

                 div(style="display:inline-block; float:right",dropdownButton(
                   HTML("<b>Input files</b>"),
                   p(downloadButton(ns('download_vcf'),""), "VCF Example File"),
                   p(downloadButton(ns('download_pheno'),""), "Trait Example File"), hr(),
                   p(HTML("<b>Parameters description:</b>"), actionButton(ns("goQploidypar"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
                   p(HTML("<b>Results description:</b>"), actionButton(ns("goQploidygraph"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
                   p(HTML("<b>How to cite:</b>"), actionButton(ns("goQploidycite"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
                   p(HTML("<b>Qploidy tutorial:</b>"), actionButton(ns("goQploidypoly"), icon("arrow-up-right-from-square", verify_fa = FALSE), onclick ="window.open('file:///Users/cht47/Documents/github/Qploidy/doc/Qploidy.html', '_blank')" )),hr(),
                   actionButton(ns("qploidy_summary"), "Summary"),
                   circle = FALSE,
                   status = "warning",
                   icon = icon("info"), width = "300px",
                   tooltip = tooltipOptions(title = "Click to see info!")
                 ))
             )
      ),
      column(width = 6,
             box(
               title = "Plots", status = "info", solidHeader = FALSE, width = 12,
               bs4Dash::tabsetPanel(
                 tabPanel("BIC Plot", plotOutput(ns("bic_plot"), height = "500px")),
                 tabPanel("Manhattan Plot", plotOutput(ns("manhattan_plot"), height = "500px")),
                 tabPanel("QQ Plot", plotOutput(ns("qq_plot"), height = "500px")),
                 tabPanel("BIC Table", DTOutput(ns("bic_table")),style = "overflow-y: auto; height: 500px"),
                 tabPanel("QTL - significant markers", DTOutput(ns("all_qtl")),style = "overflow-y: auto; height: 500px"),
                 tabPanel("Filter QTL by LD window",
                          br(),
                          box(
                            title = "LD plot",solidHeader = FALSE, width = 12, br(),
                            sliderInput(ns("bp_window_after"), label = "Adjust base pair window here to filter QTLs", min = 0,
                                        max = 100, value = 5, step = 1)
                          ),
                          box(
                            title = "Filtered QTL", solidHeader = FALSE, width = 12
                          )
                 ),
                 tabPanel("Multiple QTL model",
                          br(),
                          pickerInput(
                            inputId = ns("sele_models"),
                            label = "Select model",
                            choices = "will be updated",
                            options = list(
                              `actions-box` = TRUE),
                            multiple = FALSE
                          ), hr(),
                          pickerInput(
                            inputId = ns("sele_qtl"),
                            label = "Select QTL",
                            choices = "will be updated",
                            options = list(
                              `actions-box` = TRUE),
                            multiple = TRUE
                          ), hr(),
                          style = "overflow-y: auto; height: 500px")
               )
             )
      ),
      column(width = 3,
             box(title = "Status", width = 12, collapsible = TRUE, status = "info",
                 valueBoxOutput(ns("n.markers.start"), width=12),

                 progressBar(id = ns("pb_qploidy"), value = 0, status = "info", display_pct = TRUE, striped = TRUE, title = " ")
             ),
             box(title = "Plot Controls", status = "warning", solidHeader = TRUE, collapsible = TRUE, width = 12,
                 virtualSelectInput(
                   inputId = ns("plots"),
                   label = "Select Plot(s)*:",
                   choices = c("All",
                               "Heterozygosity",
                               "Standardized Ratio (BAF) Scatterplot",
                               "Sum of Allele Counts/Intensities Zscore",
                               "Standardized Ratio (BAF) Histogram",
                               "Non-standardized Ratio Scatterplot",
                               "Non-standardized Ratio Histogram"),
                   showValueAsTags = TRUE,
                   search = TRUE,
                   multiple = TRUE
                 ),
                 virtualSelectInput(
                   inputId = ns("sele_chr"),
                   label = "Select Chromosome(s)*:",
                   choices = NULL,
                   showValueAsTags = TRUE,
                   search = TRUE,
                   multiple = TRUE
                 ),

                 radioButtons(
                   inputId = ns("add_estimated_peaks"),
                   label = "Add Estimated Peaks:",
                   choiceNames  = c("True", "False"),
                   choiceValues = list(TRUE, FALSE),
                   selected = FALSE,
                   inline = TRUE
                 ),

                 radioButtons(
                   inputId = ns("add_expected_peaks"),
                   label = "Add Expected Peaks:",
                   choiceNames  = c("True", "False"),
                   choiceValues = list(TRUE, FALSE),
                   selected = FALSE,
                   inline = TRUE
                 ),

                 conditionalPanel(condition = "input.add_expected_peaks == 'TRUE'",
                                  ns = ns,
                                  numericInput(ns('ploidy'), label = 'Ploidy', value = 2),
                 ),

                 radioButtons(
                   inputId = ns("add_centromeres"),
                   label = "Add Centromere Position:",
                   choiceNames  = c("True", "False"),
                   choiceValues = list(TRUE, FALSE),
                   selected = FALSE,
                   inline = TRUE
                 ),

                 conditionalPanel(condition = "input.add_centromeres == 'TRUE'",
                                  ns = ns,
                                  fileInput(ns("centromeres_file"),
                                            "Centromere positions file (CSV format with columns: Chromosome, Centromere Start Position)", accept = ".csv")
                 ),

                 radioButtons(
                   inputId = ns("rm_homozygous"),
                   label = "Remove Homozygous Peaks:",
                   choiceNames  = c("True", "False"),
                   choiceValues = list(TRUE, FALSE),
                   selected = FALSE,
                   inline = TRUE
                 ),

                 numericInput(ns('dot.size'), label = 'Dot size', value = 0.05),
                 numericInput(ns('font'), label = 'Font size', value = 2),

                 div(style="display:inline-block; float:left",dropdownButton(
                   tags$h3("Save Image"),
                   selectInput(inputId = ns('image_type'), label = 'File Type', choices = c("jpeg","tiff","png","svg"), selected = "jpeg"),
                   sliderInput(inputId = ns('image_res'), label = 'Resolution', value = 300, min = 50, max = 1000, step=50),
                   sliderInput(inputId = ns('image_width'), label = 'Width', value = 9, min = 1, max = 20, step=0.5),
                   sliderInput(inputId = ns('image_height'), label = 'Height', value = 5, min = 1, max = 20, step = 0.5),
                   fluidRow(
                     downloadButton(ns("download_figure"), "Save Image"),
                     downloadButton(ns("download_file"), "Save Files"),
                     downloadButton(ns("download_stand"), "Save Standardized Data")),
                   circle = FALSE,
                   status = "danger",
                   icon = icon("floppy-disk"), width = "300px", label = "Save",
                   tooltip = tooltipOptions(title = "Click to see inputs!")
                 ))
             )
      )
    )
  )
}

#' gwas Server Functions
#'
#' @importFrom DT renderDT
#' @importFrom vcfR read.vcfR
#' @importFrom Matrix nearPD
#' @importFrom stats BIC as.formula lm logLik median model.matrix na.omit prcomp qbeta quantile runif sd setNames
#' @importFrom bs4Dash updatebs4TabItems updateBox
#' @importFrom shiny updateTabsetPanel
#' @importFrom plotly ggplotly
#' @import dplyr
#' @noRd
mod_qploidy_server <- function(input, output, session, parent_session){

  ns <- session$ns

  # Help links
  observeEvent(input$goQploidypar, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "GWAS_tabset",
                      selected = "Qploidy_par")
    # expand specific box
    updateBox(id = "Qploidy_box", action = "toggle", session = parent_session)
  })

  observeEvent(input$goQploidygraph, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "Qploidy_tabset",
                      selected = "Qploidy_results")
    # expand specific box
    updateBox(id = "Qploidy_box", action = "toggle", session = parent_session)
  })

  observeEvent(input$goQploidycite, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "Qploidy_tabset",
                      selected = "Qploidy_cite")
    # expand specific box
    updateBox(id = "Qploidy_box", action = "toggle", session = parent_session)
  })

  #UI popup window for input
  observeEvent(input$advanced_options, {
    showModal(modalDialog(
      title = "Advanced Options",
      numericInput(ns("miss"), "Remove SNPs with >= % missing data", min = 1, value = 90),
      numericInput(ns("prob"), "Remove genotypes with <= probability", min = 0, value = 0.8),
      numericInput(ns("n.clusters"), "Remove genotypes with <= # of dosage clusters",
                   min = 2, value = as.numeric(input$ref_ploidy) + 1),
      footer = tagList(
        modalButton("Close"),
        actionButton(ns("save_advanced_options"), "Save")
      )
    ))
  })

  n.markers.start <- geno.prob.rm <- miss.rm <- clusters.rm <- no.geno.info.rm <- n.markers.end <- reactiveVal(0)

  #SNP counts value box
  output$n.markers.start <- renderValueBox({
    valueBox(n.markers.start(), "Markers in uploaded file", icon = icon("dna"), color = "info")

    valueBox(
      value = tagList(
        div(style="font-size:1.2em; line-height:1.1;", span("Total markers: "), formatC(n.markers.start(), big.mark=",")),
        div(style="font-size:0.95em; opacity:0.9;", "Markers removed due:"),
        div(style="font-size:0.95em; opacity:0.9;", span("Low geno prob: "), strong(geno.prob.rm())),
        div(style="font-size:0.95em; opacity:0.9;", span("High missing data: "), strong(miss.rm())),
        div(style="font-size:0.95em; opacity:0.9;", span("Low # of dosage clusters: "), strong(clusters.rm())),
        div(style="font-size:0.95em; opacity:0.9;", span("Missing genomic position: "), strong(no.geno.info.rm())),
        div(style="font-size:1.2em; line-height:1.1;", span("Remaining markers: "), strong(n.markers.end()))
      ),
      subtitle = "",
      icon = icon("dna"),
      color = if(n.markers.end() != 0)
        if (n.markers.end()/n.markers.start() < 0.25) "darkgreen" else "darkred"
    )
  })

  ## Process inputs



}

## To be copied in the UI
# mod_qploidy_ui("qploidy_1")

## To be copied in the server
# mod_qploidy_server("qploidy_1")
