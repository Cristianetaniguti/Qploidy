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
#' @importFrom DT DTOutput
#' @importFrom bs4Dash valueBoxOutput
#' @importFrom shinydisconnect disconnectMessage
#'
mod_qploidy_ui <- function(id){
  ns <- NS(id)
  tagList(
    tags$head(tags$script(HTML("
    Shiny.addCustomMessageHandler('scroll-top', function(_) {
      window.scrollTo({ top: 0, behavior: 'smooth' });
    });
  "))),
    tags$style(HTML("#to_top { position: fixed; right: 20px; bottom: 20px; z-index: 9999; }",
                    "    /* make each tab's content white with some padding and a soft border */
    .tab-content > .tab-pane {
      background: #ffffff;
      padding: 16px;
      border: 1px solid #e6e6e6;
      border-radius: 8px;
    }")),
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
                             choices = c("VCF","Axiom Array", "Illumina Array", "Qploidy input formats","Qploidy Standardized Dataset"),
                             selected = "VCF"),
                 fileInput(ns("input_file"), "Choose File*", accept = c(".csv",".vcf",".gz", ".txt",".tsv")),
                 conditionalPanel(condition = "input.file_type == 'Axiom Array' || input.file_type == 'Illumina Array'",
                                  ns = ns,
                                  fileInput(ns("marker_pos"), "Input Marker Position file")
                 ),
                 conditionalPanel(condition = "input.file_type == 'Qploidy input formats'",
                                  ns = ns,
                                  fileInput(ns("data"), "CSV or TSV containing the following columns: MarkerName, SampleName, X, Y, R, ratio"),
                                  fileInput(ns("genos"), "CSV or TSV containing the following columns: MarkerName, SampleName, geno, prob"),
                                  fileInput(ns("geno_pos"), "CSV or TSV containing the following columns: MarkerName, Chromosome, Position")
                 ),
                 conditionalPanel(condition = "input.file_type != 'Qploidy Standardized Dataset'",
                                  ns=ns,
                                  radioButtons(ns("known_ploidy"),
                                               label = "Do you have >=60 samples with the same, known ploidy?",
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
                                  numericInput(ns("cores"), "Number of CPU Cores*", min = 1, value = 1),
                                  conditionalPanel(condition = "input.known_ploidy == 'TRUE'",
                                                   ns=ns,
                                                   numericInput(ns("ref_ploidy"), "Reference Samples Ploidy*", min = 1, value = NULL)),
                                  conditionalPanel(condition = "input.known_ploidy == 'FALSE'",
                                                   ns=ns,
                                                   numericInput(ns("common_ploidy"), "Most Common Ploidy*", min = 1, value = NULL)),
                                  div(style="text-align: left; margin-top: 10px;",
                                      actionButton(ns("advanced_options"),
                                                   label = HTML(paste(icon("cog", style = "color: #007bff;"), "Advanced Options")),
                                                   style = "background-color: transparent; border: none; color: #007bff; font-size: smaller; text-decoration: underline; padding: 0;"
                                      )
                                  ),
                 ),br(),
                 actionButton(ns("run_stand"), "Run Analysis"),

                 div(style="display:inline-block; float:right",dropdownButton(
                   HTML("<b>Input files</b>"),
                   p(downloadButton(ns('download_vcf'),""), "VCF Example File"),
                   p(downloadButton(ns('download_axiom'),""), "Axiom Summary Example File"), hr(),
                   p(downloadButton(ns('download_fitpoly'),""), "fitpoly Result Example File"), hr(),
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
             ),
             box(title = "Status", width = 12, collapsible = TRUE, status = "info",
                 valueBoxOutput(ns("markers_status"), width=12),
                 progressBar(id = ns("pb_qploidy"), value = 0, status = "info", display_pct = TRUE, striped = TRUE, title = " ")

             ),
      ),
      column(width = 9,
             bs4Dash::tabsetPanel(
               id = "Qploidy_results",
               tabPanel("Standardization \n Plots by Sample", value = "stand_by_sample", br(),
                        tags$div(
                          role = "note", `aria-live` = "polite",
                          style = "margin-bottom:10px;",
                          tags$p(
                            tags$strong(icon("chart-line"), " Explore Standardized Plots. "),
                            "Interactively visualize standardized allele ratios, BAF scatterplots, z-scores, and more for each sample. Use the controls below to select samples, chromosomes, and plot types."
                          ),
                          tags$ul(
                            tags$li(" Select a sample and chromosome/s to view its plots."),
                            tags$li(" Adjust plot options to improve visualization."),
                            tags$li(" Review diagnostic plots before further analysis.")
                          )
                        ),
                        box(title = "Plot Controls", status = "info", solidHeader = FALSE, collapsible = TRUE, width = 12,
                            fluidRow(
                              column(width = 6,
                                     virtualSelectInput(
                                       inputId = ns("sample"),
                                       label = "Select Sample*:",
                                       choices = NULL,
                                       showValueAsTags = TRUE,
                                       search = TRUE,
                                       multiple = FALSE
                                     ),
                                     virtualSelectInput(
                                       inputId = ns("plots"),
                                       label = "Select Plot(s)*:",
                                       choices = c("Heterozygosity" = "het",
                                                   "Standardized Ratio (BAF) Scatterplot" = "BAF",
                                                   "Sum of Allele Counts/Intensities Zscore" = "zscore",
                                                   "Standardized Ratio (BAF) Histogram" = "BAF_hist",
                                                   "Non-standardized Ratio Scatterplot" = "ratio",
                                                   "Sample-level Non-standardized Ratio Histogram" = "Ratio_hist_overall",
                                                   "Sample-level Standardized Ratio (BAF) Histogram" = "BAF_hist_overall"),
                                       showValueAsTags = TRUE,
                                       search = TRUE,
                                       multiple = TRUE
                                     ),
                                     virtualSelectInput(
                                       inputId = ns("sele_chr_standardization"),
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
                                                      numericInput(ns('ploidy_standardization'), label = 'Ploidy', value = 2),
                                     )
                              ),
                              column(width = 6,
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
                                                      fileInput(ns("centromeres_file_standardization"),
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
                                     actionButton(ns("build_plot"), "Build plot"),
                                     div(style="float:right",dropdownButton(
                                       tags$h3("Save Image"),
                                       selectInput(inputId = ns('image_type'), label = 'File Type', choices = c("jpeg","tiff","png","svg"), selected = "jpeg"),
                                       sliderInput(inputId = ns('image_res'), label = 'Resolution', value = 300, min = 50, max = 1000, step=50),
                                       sliderInput(inputId = ns('image_width'), label = 'Width', value = 9, min = 1, max = 20, step=0.5),
                                       sliderInput(inputId = ns('image_height'), label = 'Height', value = 5, min = 1, max = 20, step = 0.5),
                                       fluidRow(
                                         downloadButton(ns("download_figure"), "Save Image"), br(),
                                         downloadButton(ns("download_stand"), "Save Standardized Data")),
                                       circle = FALSE,
                                       status = "danger",
                                       icon = icon("floppy-disk"), width = "300px", label = "Save",
                                       tooltip = tooltipOptions(title = "Click to see inputs!")
                                     )
                                     )
                              )
                            ), hr(),

                            uiOutput(ns("plot_ui"))
                        )
               ),
               tabPanel(title = "HMM by Sample (BETA)", value = "hmm_by_sample", br(),
                        tags$div(
                          role = "note", `aria-live` = "polite",
                          style = "margin-bottom:10px;",
                          tags$p(
                            tags$strong(icon("flask"), " This feature is in beta. "),
                            "The multipoint Hidden Markov Model (HMM) copy-number estimation is under active development and testing.",
                            "Results may change between versions and should be treated as provisional."
                          ),
                          tags$ul(
                            tags$li("The algorithm operates on standardized values."),
                            tags$li("Before running, inspect the diagnostic plots and ensure they look reasonable."),
                            tags$li("Adjust the parameters below as needed, then click ",
                                    tags$code("Run HMM"), " to start the estimation.")
                          )
                        ),
                        hr(),
                        fluidRow(
                          column(width = 6,
                                 virtualSelectInput(
                                   inputId = ns("sample_hmm"),
                                   label = "Select Sample*:",
                                   choices = NULL,
                                   showValueAsTags = TRUE,
                                   search = TRUE,
                                   multiple = FALSE
                                 ),
                                 virtualSelectInput(
                                   inputId = ns("sele_chr_hmm"),
                                   label = "Select Chromosome(s)*:",
                                   choices = NULL,
                                   showValueAsTags = TRUE,
                                   search = TRUE,
                                   multiple = TRUE
                                 ),

                                 textInput(ns("ploidy_range_hmm"), "Ploidies to be tested", placeholder = "2,3,4", value = "2,3,4"), br(),
                                 actionButton(ns("est_ploidy_hmm"), "Run HMM"), br(),
                          ),
                          column(width = 6,
                                 numericInput(ns("snps_per_window_hmm"), "Number of SNPs per window", min = 5, value = 20),
                                 numericInput(ns("exp_ploidy_hmm"), "Expected overall ploidy", min = 1, value = 2),
                                 actionButton(ns("advanced_options_hmm"),
                                              label = HTML(paste(icon("cog", style = "color: #007bff;"), "Advanced Options")),
                                              style = "background-color: transparent; border: none; color: #007bff; font-size: smaller; text-decoration: underline; padding: 0;"
                                 )
                          )
                        ),
                        br(), hr(),
                        plotOutput(ns("plot_hmm"), height = 800), br(),
                        box(title="HMM results", collapsed = TRUE, width = 12,
                            DTOutput(ns("ploidy_table_hmm")),br(),
                            downloadButton(ns("download_ploidy_table_hmm"), "Download HMM CN Estimations Table")
                        )
               ),
               tabPanel(title = "All Samples Ploidy", value = "all_samples_ploidy", br(),
                        tags$div(
                          role = "note", `aria-live` = "polite",
                          style = "margin-bottom:10px;",
                          tags$p(
                            tags$strong(icon("users"), " Estimate Ploidy for All Samples. "),
                            "Review the plots in the 'Standardization Plots by Sample' and 'HMM by Sample (BETA)' tabs before estimating ploidies for all samples. Once satisfied, adjust the parameters below and click ",
                            tags$code("Estimate All Samples Ploidy"), " to run the analysis."
                          ),
                          tags$ul(
                            tags$li(" Inspect standardization and HMM plots for quality control."),
                            tags$li(" Set analysis parameters and resolution level."),
                            tags$li(" Review and download the results table below.")
                          )
                        ),
                        hr(),
                        textInput(ns("ploidy_range"), "Ploidies to be tested", placeholder = "2,3,4", value = "2,3,4"),
                        selectInput(ns('res_lvl'),
                                    label = 'Select Resolution Level',
                                    choices = c("Sample" = "sample",
                                                "Chromosome" = "chromosome",
                                                "Chromosome-arm" = "chromosome-arm"),
                                    selected = "Chromosome"),
                        conditionalPanel(condition = "input.res_lvl == 'chromosome-arm'",
                                         ns = ns,
                                         fileInput(ns("centromeres_file2"),
                                                   "Centromere positions file (CSV format with columns: Chromosome, Centromere Start Position)", accept = ".csv")
                        ),
                        checkboxInput(ns("all_hmm"), "Run HMM CN estimations (beta)", value = FALSE),
                        conditionalPanel(condition = "input.all_hmm",
                                         ns = ns,
                                         fluidRow(
                                           column(width = 6,
                                                  virtualSelectInput(
                                                    inputId = ns("sele_chr_all"),
                                                    label = "Select Chromosome(s)*:",
                                                    choices = NULL,
                                                    showValueAsTags = TRUE,
                                                    search = TRUE,
                                                    multiple = TRUE
                                                  ),
                                                  numericInput(ns("exp_ploidy_all"), "Expected overall ploidy", min = 1, value = 2),
                                           ),
                                           column(width = 6,
                                                  numericInput(ns("snps_per_window_all"), "Number of SNPs per window", min = 5, value = 20),
                                                  numericInput(ns("n_cores_hmm"), "The process is paralellized by sampl. Define the number of CPU cores to be used:",
                                                               min = 1, value = 1),
                                                  actionButton(ns("advanced_options_hmm_all"),
                                                               label = HTML(paste(icon("cog", style = "color: #007bff;"), "Advanced Options")),
                                                               style = "background-color: transparent; border: none; color: #007bff; font-size: smaller; text-decoration: underline; padding: 0;"
                                                  )
                                           )

                                         )
                        ),
                        br(),
                        actionButton(ns("est_ploidy"), "Estimate All Samples Ploidy"), hr(), br(),
                        uiOutput(ns("ploidy_table_note")),
                        DTOutput(ns("ploidy_table")),br(),
                        downloadButton(ns("download_ploidy_table"), "Download Table")
               )
             )
      )
    )
  )
}

#' qploidy Server Functions
#'
#' @importFrom DT renderDT datatable
#' @importFrom bs4Dash updatebs4TabItems updateBox
#' @importFrom shiny updateTabsetPanel
#' @import dplyr
#' @importFrom shinyalert shinyalert
#'
#' @noRd
mod_qploidy_server <- function(input, output, session, parent_session){

  ns <- session$ns

  # Helper to check file input existence and non-empty
  check_file_input <- function(x, name) {
    if (is.null(x) || is.null(x$datapath) || !file.exists(x$datapath)) {
      shinyalert(
        title = paste(name, "Missing"),
        text = paste("Please upload a valid", name, "file."),
        type = "error"
      )
      return(FALSE)
    }
    TRUE
  }

  # Helper to check required columns in a data.frame
  check_required_columns <- function(df, required, name) {
    if (!all(required %in% colnames(df))) {
      shinyalert(
        title = paste(name, "Format Error"),
        text = paste("The", name, "file must contain the following columns:", paste(required, collapse = ", ")),
        type = "error"
      )
      return(FALSE)
    }
    TRUE
  }

  # Help links
  observeEvent(input$goQploidypar, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "Qploidy_tabset",
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

  observeEvent(list(input$ref_ploidy, input$common_ploidy), {
    # Update n.clusters reactively using input values
    ploidy <- c(input$ref_ploidy, input$common_ploidy)
    ploidy <- ploidy[!is.na(ploidy)]
    if(length(ploidy) > 0) {
      advanced_options$n.clusters <- as.numeric(ploidy) + 1
    } else {
      advanced_options$n.clusters <- input$n.clusters
    }
  })

  # Initialize with static defaults only
  advanced_options <- reactiveValues(
    miss = 90,
    prob = 0.8,
    n.clusters = 3 # default value, will be updated reactively
  )

  advanced_options_hmm <- reactiveValues(
    min_snps_per_window = 20,
    n_bins = 100,
    bw = 0.03,
    het_lims = c(0,1),
    het_weight = 0.8,
    z_range = 0.2,
    transition_jump =0.995,
    z_only = FALSE
  )

  #UI popup window for input
  observeEvent(input$advanced_options, {
    showModal(modalDialog(
      title = "Standardization Advanced Options",
      numericInput(ns("miss"), "Remove SNPs with >= % missing data", min = 1, max = 100,value = advanced_options$miss),
      numericInput(ns("prob"), "Remove genotypes with <= probability", min = 0, value = advanced_options$prob),
      numericInput(ns("n.clusters"), "Remove genotypes with <= # of dosage clusters",
                   min = 2, value = advanced_options$n.clusters),
      footer = tagList(
        modalButton("Close"),
        actionButton(ns("save_advanced_options"), "Save")
      )
    ))
  })

  #UI popup window for input
  observeEvent(input$advanced_options_hmm, {
    showModal(modalDialog(
      title = "Advanced Options for HMM CN Estimation",
      numericInput(ns("min_snps_per_window_hmm"), "Minimum number of SNPs per window",
                   min = 1, value = advanced_options_hmm$min_snps_per_window),
      numericInput(ns("n_bins_hmm"), "Number of bins to split BAF",
                   min = 5, value = advanced_options_hmm$n_bins),
      numericInput(ns("bw_hmm"), "Bandwidth (SD) of the Gaussian kernels used to generate the BAF comb templates",
                   min = 0, value = advanced_options_hmm$bw),
      sliderInput(ns("het_lims_hmm"), "BAF limits to consider a SNP heterozygous",
                  min = 0, max = 1, value = advanced_options_hmm$het_lims, step = 0.005),
      numericInput(ns("het_weight_hmm"), "Quantile used to scale BAF emission weight based on heterozygote count",
                   min = 0, max = 1, value = advanced_options_hmm$het_weight),
      numericInput(ns("z_range_hmm"), "Padding added to min/max z for initial mean estimation",
                   min = 0, value = advanced_options_hmm$z_range),
      numericInput(ns("transition_jump_hmm"), "Diagonal value for transition matrix (probability to stay in same CN state)",
                   min = 0, max = 1,value = advanced_options_hmm$transition_jump),
      checkboxInput(ns("z_only_hmm"), "Fit the HMM using the z-emission only (ignores BAF)", value = advanced_options_hmm$z_only),
      footer = tagList(
        modalButton("Close"),
        actionButton(ns("save_advanced_options_hmm"), "Save")
      )
    ))
  })

  observeEvent(input$advanced_options_hmm_all, {
    showModal(modalDialog(
      title = "Advanced Options for HMM CN Estimation",
      numericInput(ns("min_snps_per_window_all"), "Minimum number of SNPs per window",
                   min = 1, value = advanced_options_hmm$min_snps_per_window),
      numericInput(ns("n_bins_all"), "Number of bins to split BAF",
                   min = 5, value = advanced_options_hmm$n_bins),
      numericInput(ns("bw_all"), "Bandwidth (SD) of the Gaussian kernels used to generate the BAF comb templates",
                   min = 0, value = advanced_options_hmm$bw),
      sliderInput(ns("het_lims_all"), "BAF limits to consider a SNP heterozygous",
                  min = 0, max = 1, value = advanced_options_hmm$het_lims , step = 0.005),
      numericInput(ns("het_weight_all"), "Quantile used to scale BAF emission weight based on heterozygote count",
                   min = 0, max = 1, value = advanced_options_hmm$het_weight),
      numericInput(ns("z_range_all"), "Padding added to min/max z for initial mean estimation",
                   min = 0, value = advanced_options_hmm$z_range),
      numericInput(ns("transition_jump_all"), "Diagonal value for transition matrix (probability to stay in same CN state)",
                   min = 0, max = 1,value = advanced_options_hmm$transition_jump),
      checkboxInput(ns("z_only_all"), "Fit the HMM using the z-emission only (ignores BAF)", value = advanced_options_hmm$z_only),
      footer = tagList(
        modalButton("Close"),
        actionButton(ns("save_advanced_options_hmm_all"), "Save")
      )
    ))
  })

  # Close popup window when user "saves options"
  observeEvent(input$save_advanced_options, {
    advanced_options$miss <- input$miss
    advanced_options$prob <- input$prob
    advanced_options$n.clusters <- input$n.clusters

    removeModal() # Close the modal after saving
  })

  observeEvent(input$save_advanced_options_hmm, {
    advanced_options_hmm$min_snps_per_window <- input$min_snps_per_window_hmm
    advanced_options_hmm$n_bins <- input$n_bins_hmm
    advanced_options_hmm$bw <- input$bw_hmm
    advanced_options_hmm$het_lims <- input$het_lims_hmm
    advanced_options_hmm$het_weight <- input$het_weight_hmm
    advanced_options_hmm$z_range <- input$z_range_hmm
    advanced_options_hmm$transition_jump <- input$transition_jump_hmm
    advanced_options_hmm$z_only <- input$z_only_hmm

    removeModal() # Close the modal after saving
  })

  observeEvent(input$save_advanced_options_hmm_all, {
    advanced_options_hmm$min_snps_per_window <- input$min_snps_per_window_hmm_all
    advanced_options_hmm$n_bins <- input$n_bins_hmm_all
    advanced_options_hmm$bw <- input$bw_hmm_all
    advanced_options_hmm$het_lims <- input$het_lims_hmm_all
    advanced_options_hmm$het_weight <- input$het_weight_hmm_all
    advanced_options_hmm$z_range <- input$z_range_hmm_all
    advanced_options_hmm$transition_jump <- input$transition_jump_hmm_all
    advanced_options_hmm$z_only <- input$z_only_hmm_all

    removeModal() # Close the modal after saving
  })

  n.markers.start <- geno.prob.rm <- miss.rm <- clusters.rm <- no.geno.info.rm <- n.markers.end <- reactiveVal(0)

  #SNP counts value box
  output$markers_status <- renderValueBox({
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

  stand_file <- reactiveVal(NULL)

  ## Process inputs
  data_standardized <- eventReactive(input$run_stand, {
    req(input$file_type)

    # Reset progress bar
    updateProgressBar(session = session, id = "pb_qploidy", value = 5, total = 100)

    # File existence checks
    if (!check_file_input(input$input_file, "Input")) return()
    if (input$file_type == "Axiom Array" || input$file_type == "Illumina Array") {
      if (!check_file_input(input$marker_pos, "Marker Position")) return()
    }
    if (input$file_type == "Qploidy input formats") {
      if (!check_file_input(input$data, "Data")) return()
      if (!check_file_input(input$genos, "Genotypes")) return()
      if (!check_file_input(input$geno_pos, "Genotype Position")) return()
    }
    if (input$file_type != "Qploidy Standardized Dataset") {
      if (input$file_type == "VCF" && identical(input$known_ploidy, TRUE)) {
        if (!check_file_input(input$reference_samples, "Reference Samples")) return()
      }
      if ((input$file_type == "Axiom Array" || input$file_type == "Illumina Array") && identical(input$known_ploidy, TRUE)) {
        if (!check_file_input(input$geno_ref, "fitpoly Reference Scores")) return()
      }
      if ((input$file_type == "Axiom Array" || input$file_type == "Illumina Array") && identical(input$known_ploidy, FALSE)) {
        if (!check_file_input(input$geno_all, "fitpoly All Scores")) return()
      }
    }

    # Update progress bar
    updateProgressBar(session = session, id = "pb_qploidy", value = 20)

    # Read input file based on type
    if (input$file_type == "VCF") {
      input_data <- qploidy_read_vcf(input$input_file$datapath)
    } else if (input$file_type == "Axiom Array") {
      input_data <- read_axiom(input$input_file$datapath)
    } else if(input$file_type == "Illumina Array") {
      input_data <- read_illumina_array(input$input_file$datapath)
    } else if (input$file_type == "Qploidy Standardized Dataset") {
      input_data <- "skip"
      data_standardized <- read_qploidy_standardization(input$input_file$datapath)
    } else if (input$file_type == "Qploidy input formats") {
      input_data <- read.csv(input$data$datapath, header = TRUE, sep = if(grepl(".tsv$", input$data$name)) "\t" else ",")
      genos <- read.csv(input$genos$datapath, header = TRUE, sep = if(grepl(".tsv$", input$genos$name)) "\t" else ",")
      geno.pos <- read.csv(input$geno_pos$datapath, header = TRUE, sep = if(grepl(".tsv$", input$geno_pos$name)) "\t" else ",")
      # Check required columns
      req_cols_data <- c("MarkerName", "SampleName", "X", "Y", "R", "ratio")
      req_cols_genos <- c("MarkerName", "SampleName", "geno", "prob")
      req_cols_geno_pos <- c("MarkerName", "Chromosome", "Position")
      if (!check_required_columns(input_data, req_cols_data, "Data")) return()
      if (!check_required_columns(genos, req_cols_genos, "Genotypes")) return()
      if (!check_required_columns(geno.pos, req_cols_geno_pos, "Genotype Position")) return()
    }

    if(input$file_type != "Qploidy Standardized Dataset"){

      if(input$file_type == "VCF"){
        genos <- qploidy_read_vcf(input$input_file$datapath, geno = TRUE)
        genos.pos <- qploidy_read_vcf(input$input_file$datapath, geno.pos = TRUE)

        if(all(is.na(genos.pos[,1]))) genos.pos$MarkerName <- paste0(genos.pos$Chromosome, "_", genos.pos$Position)

        if(input$known_ploidy){
          reference_samples <- read.csv(input$reference_samples$datapath, header = TRUE, stringsAsFactors = FALSE)$Sample_ID
          genos <- genos[which(genos$SampleName %in% reference_samples), ]
        }

      } else if (input$file_type == "Axiom Array" || input$file_type == "Illumina Array") {
        if(input$known_ploidy)
          fitpoly_scores <- read.table(input$geno_ref$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        else
          fitpoly_scores <- read.table(input$geno_all$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

        genos <- data.frame(
          MarkerName = fitpoly_scores$MarkerName,
          SampleName = fitpoly_scores$SampleName,
          geno = fitpoly_scores$maxgeno,
          prob = fitpoly_scores$maxP
        )

        genos.pos_file <- read.table(input$marker_pos, header = T)

        genos.pos <- data.frame(
          MarkerName = genos.pos_file$probes,
          Chromosome = genos.pos_file$chr,
          Position = genos.pos_file$pos
        )
      }

      # Verify if informed ploidy is the same in file
      ploidy <- c(input$ref_ploidy, input$common_ploidy)
      ploidy <- ploidy[-is.na(ploidy)]

      if (length(ploidy) > 1) {
        if (!all(max(genos$geno, na.rm = TRUE) %in% ploidy)) {
          shinyalert(
            title = "Ploidy Mismatch",
            text = "Genotype on File doesn't match any informed Ploidy.",
            type = "error"
          )
          return()
        }
      } else {
        if (max(genos$geno, na.rm = TRUE) != ploidy) {
          shinyalert(
            title = "Ploidy Mismatch",
            text = "Genotype on File doesn't match Ploidy informed",
            type = "error"
          )
          return()
        }
      }

      # Updated progress bar
      updateProgressBar(session = session, id = "pb_qploidy", value = 40)


      # Create reactive value for temp path output file
      tmp <- tempfile(fileext = ".tsv.gz")

      data_standardized <- standardize(
        data = input_data,
        genos = genos,
        geno.pos = genos.pos,
        ploidy.standardization = ploidy,
        threshold.geno.prob = if(is.null(input$prob)) 0.8 else input$prob,
        threshold.missing.geno = if(is.null(input$miss)) 0.9 else input$miss/100,
        threshold.n.clusters = if(is.null(input$n.clusters)) ploidy + 1 else input$n.clusters,
        n.cores = input$cores,
        out_filename = tmp,
        verbose = FALSE
      )
      stand_file(tmp)   # remember the path we wrote

    } else {
      data_standardized <- read_qploidy_standardization(input$input_file$datapath)
    }

    # Update progress bar
    updateProgressBar(session = session, id = "pb_qploidy", value = 100)

    data_standardized
  })

  output$markers_status <- renderValueBox({
    req(data_standardized())
    total <- data_standardized()$filters[1]
    remain <- data_standardized()$filters[6]
    color <- "success"
    if (is.null(total) || total == 0) color <- "danger"
    else if (remain/total < 0.25) color <- "danger"
    valueBox(
      value = tagList(
        div(style="font-size:1.2em; line-height:1.1;", span("Total markers: "), formatC(total, big.mark=",")),
        div(style="font-size:0.95em; opacity:0.9;", "Markers removed due:"),
        div(style="font-size:0.95em; opacity:0.9;", span("Low geno prob: "), strong(data_standardized()$filters[2])),
        div(style="font-size:0.95em; opacity:0.9;", span("High missing data: "), strong(data_standardized()$filters[3])),
        div(style="font-size:0.95em; opacity:0.9;", span("Low # of dosage clusters: "), strong(data_standardized()$filters[4])),
        div(style="font-size:0.95em; opacity:0.9;", span("Missing genomic position: "), strong(data_standardized()$filters[5])),
        div(style="font-size:1.2em; line-height:1.1;", span("Remaining markers: "), strong(remain))
      ),
      subtitle = "",
      icon = icon("dna"),
      color = color
    )
  })


  observe({
    req(data_standardized())
    samples <- unique(data_standardized()$data$SampleName)
    updateVirtualSelect(
      session = session,
      inputId = "sample",
      choices = samples,
      selected = samples[1]
    )

    updateVirtualSelect(
      session = session,
      inputId = "sample_hmm",
      choices = samples,
      selected = samples[1]
    )

    chrs <- unique(data_standardized()$data$Chr)
    if(any(is.na(chrs))) chrs <- chrs[-which(is.na(chrs))]
    chrs <- sort(chrs)

    sel <- if(length(chrs) <= 5) chrs else chrs[1:5]
    names(chrs) <- chrs

    updateVirtualSelect(
      session = session,
      inputId = "sele_chr_standardization",
      choices = chrs,
      selected = sel
    )

    updateVirtualSelect(
      session = session,
      inputId = "sele_chr_hmm",
      choices = chrs,
      selected = sel
    )

    updateVirtualSelect(
      session = session,
      inputId = "sele_chr_all",
      choices = chrs,
      selected = sel
    )
  })

  ploidies <- eventReactive(input$est_ploidy, {
    req(data_standardized())

    updateProgressBar(session = session, id = "pb_qploidy", value = 5)

    if(!is.null(input$centromeres_file2$datapath)){
      centromeres <- read.csv(input$centromeres_file2$datapath, header = TRUE, stringsAsFactors = FALSE)
      centromeres_vec <- centromeres[,2]
      names(centromeres_vec) <- centromeres[,1]
    } else {
      centromeres_vec <- NULL
    }

    ploidies <- as.numeric(unlist(strsplit(input$ploidy_range, ",")))
    estimated_ploidies <- area_estimate_ploidy(
      qploidy_standardization = data_standardized(),
      samples = "all",
      level = input$res_lvl,
      ploidies = ploidies,
      centromeres = if(!is.null(centromeres_vec)) centromeres_vec
    )

    updateProgressBar(session = session, id = "pb_qploidy", value = 20)

    if(input$all_hmm){

      multi_esti <- hmm_estimate_CN_multi(qploidy_standarize_result = data_standardized(),
                                          sample_ids = "all",
                                          n_cores = if(is.null(input$n_cores_hmm)) 1 else input$n_cores_hmm,
                                          chr = input$sele_chr_all,
                                          snps_per_window = if(is.null(input$snps_per_window_all)) 20 else input$snps_per_window_all,
                                          min_snps_per_window = if(is.null(advanced_options_hmm$min_snps_per_window)) 20 else advanced_options_hmm$min_snps_per_window,
                                          cn_grid = if(is.null(ploidies)) c(2,3,4) else ploidies,
                                          M = if(is.null(advanced_options_hmm$n_bins)) 100 else advanced_options_hmm$n_bins,
                                          bw = if(is.null(advanced_options_hmm$bw)) 0.03 else advanced_options_hmm$bw,
                                          exp_ploidy = if(is.null(input$exp_ploidy_all)) 2 else input$exp_ploidy_all,
                                          het_lims = if(is.null(advanced_options_hmm$het_lims)) c(0,1) else advanced_options_hmm$het_lims,
                                          het_weight = if(is.null(advanced_options_hmm$het_weight)) 0.8 else advanced_options_hmm$het_weight,
                                          z_range = if(is.null(advanced_options_hmm$z_range)) 0.2 else advanced_options_hmm$z_range,
                                          transition_jump = if(is.null(advanced_options_hmm$transition_jump)) 0.995 else advanced_options_hmm$transition_jump,
                                          max_iter = 60,
                                          z_only = if(is.null(input$z_only)) FALSE else input$z_only)

      updateProgressBar(session = session, id = "pb_qploidy", value = 60)

      summ_hmm <- summarize_cn_mode(multi_esti, level = input$res_lvl)
      merged_table <- merge_cn_summary_with_estimates(summ_hmm, estimated_ploidies, level=input$res_lvl)

    } else merged_table <- NULL

    updateProgressBar(session = session, id = "pb_qploidy", value = 100)

    list(area = estimated_ploidies, hmm = merged_table)
  })

  output$ploidy_table <- renderDT({
    req(ploidies())
    if(input$all_hmm){
      tab <- ploidies()$hmm
      num_cols <- sapply(tab, is.numeric)
      tab[num_cols] <- lapply(tab[num_cols], function(x) round(x, 4))
      datatable(tab,
                selection = "single",
                options = list(pageLength = 10, scrollX = TRUE))
    } else {
      ploidy <- ploidies()$area$ploidy
      ploidy[which(ploidies()$area$diff_first_second < 0.01)] <- NA

      if(ncol(ploidy) > 1){
        num_cols <- sapply(ploidy, is.numeric)
        ploidy[num_cols] <- lapply(ploidy[num_cols], function(x) round(x, 4))
      } else {
        ploidy <- as.data.frame(round(ploidy, 4))
        rownames(ploidy) <- rownames(ploidies()$area$ploidy)
        colnames(ploidy) <- colnames(ploidies()$area$ploidy)
      }

      datatable(ploidy,
                selection = "single",
                options = list(pageLength = 10, scrollX = TRUE))
    }
  })

  output$ploidy_table_note <- renderUI({
    req(ploidies())
    tab <- if (input$all_hmm) ploidies()$hmm else ploidies()$area$ploidy
    if (is.null(tab) || nrow(tab) == 0) return(NULL) else {
      tags$div(
        style = "margin-bottom:10px; background-color:#f8f9fa; border-radius:6px; padding:10px; border:1px solid #e2e3e5;",
        tags$strong(icon("mouse-pointer"), " Tip: "),
        "Click a row in the table below to update the selected sample in the plot and HMM estimation panels."
      )
    }
  })

  picked_samples <- eventReactive(input$ploidy_table_rows_selected, {
    req(ploidies())
    s <- input$ploidy_table_rows_selected
    req(length(s) == 1)

    if(input$all_hmm) sample_picked <- ploidies()$hmm[s,1] else sample_picked <- rownames(ploidies()$area$ploidy)[s]
    sample_picked
  })

  observeEvent(picked_samples(),{
    updateVirtualSelect(
      session = session,
      inputId = "sample",
      selected = picked_samples()
    )

    updateVirtualSelect(
      session = session,
      inputId = "sample_hmm",
      selected = picked_samples()
    )

      # select specific tab
      updateTabsetPanel(session = parent_session, inputId = "Qploidy_results",
                        selected = "stand_by_sample")

  })

  output$download_ploidy_table <- downloadHandler(
    filename = function() {
      paste0("ploidy_estimation_", Sys.Date(), ".tsv")
    },
    content = function(file) {
      req(ploidies())
      if(input$all_hmm){
        write.table(ploidies()$hmm, file, sep = "\t", row.names = FALSE, quote = FALSE)
      } else {
        ploidy <- ploidies()$area$ploidy
        write.table(ploidy, file, sep = "\t", row.names = TRUE, quote = FALSE)
      }
    }
  )

  built_plot <- eventReactive(input$build_plot, {
    req(data_standardized())
    req(input$sample)
    req(input$plots)
    req(input$sele_chr_standardization)

    # Update progress bar
    updateProgressBar(session = session, id = "pb_qploidy", value = 0, total = 100)

    if(!is.null(input$centromeres_file_standardization$datapath)){
      centromeres <- read.csv(input$centromeres_file_standardization$datapath, header = TRUE, stringsAsFactors = FALSE)
      centromeres_vec <- centromeres[,2]
      names(centromeres_vec) <- centromeres[,1]
    } else {
      centromeres_vec <- NULL
    }

    updateProgressBar(session = session, id = "pb_qploidy", value = 50)

    p <- plot_qploidy_standardization(
      x = data_standardized(),
      sample = input$sample,
      type = input$plots,
      chr = input$sele_chr_standardization,
      ploidy = input$ploidy_standardization,
      add_expected_peaks = as.logical(input$add_expected_peaks),
      add_estimated_peaks = as.logical(input$add_estimated_peaks),
      rm_homozygous = as.logical(input$rm_homozygous),
      add_centromeres = as.logical(input$add_centromeres),
      centromeres = centromeres_vec
    )


    # Update progress bar
    updateProgressBar(session = session, id = "pb_qploidy", value = 100)
    p
  })

  output$plot <- renderPlot({
    built_plot()
  }, res = 96)

  output$plot_ui <- renderUI({
    n <- max(1L, length(input$plots))        # treat 0 as 1
    base_h <- 350                             # px per plot “unit”
    plotOutput(session$ns("plot"), height = paste0(base_h * n, "px"))
  })

  ## Add HMM
  ploidies_hmm <- eventReactive(input$est_ploidy_hmm, {
    req(data_standardized())
    updateProgressBar(session = session, id = "pb_qploidy", value = 10)

    ploidies <- as.numeric(unlist(strsplit(input$ploidy_range_hmm, ",")))

    esti <- hmm_estimate_CN(
      qploidy_standarize_result = data_standardized(),
      sample_id = input$sample_hmm,
      chr = input$sele_chr_hmm,
      snps_per_window = if(is.null(input$snps_per_window_hmm)) 20 else input$snps_per_window_hmm,
      min_snps_per_window = if(is.null(advanced_options_hmm$min_snps_per_window)) 20 else advanced_options_hmm$min_snps_per_window,
      cn_grid = if(is.null(ploidies)) c(2,3,4) else ploidies,
      M = if(is.null(advanced_options_hmm$n_bins)) 100 else advanced_options_hmm$n_bins,
      bw = if(is.null(advanced_options_hmm$bw)) 0.03 else advanced_options_hmm$bw,
      exp_ploidy = if(is.null(input$exp_ploidy_hmm)) 2 else input$exp_ploidy_hmm,
      het_lims = if(is.null(advanced_options_hmm$het_lims)) c(0,1) else advanced_options_hmm$het_lims,
      het_weight = if(is.null(advanced_options_hmm$het_weight)) 0.8 else advanced_options_hmm$het_weight,
      z_range = if(is.null(advanced_options_hmm$z_range)) 0.2 else advanced_options_hmm$z_range,
      transition_jump = if(is.null(advanced_options_hmm$transition_jump)) 0.995 else advanced_options_hmm$transition_jump,
      max_iter = 60,
      z_only = if(is.null(advanced_options_hmm$z_only)) FALSE else advanced_options_hmm$z_only
    )

    updateProgressBar(session = session, id = "pb_qploidy", value = 75)

    esti
  })

  output$ploidy_table_hmm <- renderDT({
    req(ploidies_hmm())
    updateProgressBar(session = session, id = "pb_qploidy", value = 85)

    tab <- ploidies_hmm()[[1]]
    # Round all numeric columns to 4 decimals
    num_cols <- sapply(tab, is.numeric)
    tab[num_cols] <- lapply(tab[num_cols], function(x) round(x, 4))

    datatable(
      tab,
      options = list(
        scrollX   = TRUE,
        autoWidth = TRUE
      ),
      rownames = FALSE
    )
  })

  output$plot_hmm <- renderPlot({
    req(ploidies_hmm())

    p <- plot_cn_track(hmm_CN = ploidies_hmm(), 
                      qploidy_standarize_result= data_standardized(), 
                      sample_id = input$sample_hmm, 
                      show_window_lines = TRUE)

    updateProgressBar(session = session, id = "pb_qploidy", value = 100)
    p
  })

  output$download_stand <- downloadHandler(
    filename = function() paste0("qploidy_standardized_data_", Sys.Date(), ".tsv.gz"),
    content  = function(file) {
      req(data_standardized())                                  # ensure it ran
      src <- stand_file()
      validate(need(!is.null(src) && file.exists(src), "File not generated yet."))
      ok <- file.copy(from = src, to = file, overwrite = TRUE)
      validate(need(ok, "Failed to copy the standardized file."))
    }
  )

  output$download_figure <- downloadHandler(
    filename = function() {
      ext <- switch(input$image_type,
                    "jpeg" = "jpg",
                    "png"  = "png",
                    "tiff" = "tiff",
                    "svg"  = "svg")  # assume UI sends one of these
      paste0("Qploidy-", Sys.Date(), ".", ext)
    },
    content = function(file) {
      # open the right device
      switch(input$image_type,
             "jpeg" = jpeg(file,
                           width  = as.numeric(input$image_width),
                           height = as.numeric(input$image_height),
                           units  = "in", res = as.numeric(input$image_res)),
             "png"  = png(file,
                          width  = as.numeric(input$image_width),
                          height = as.numeric(input$image_height),
                          units  = "in", res = as.numeric(input$image_res)),
             "tiff" = tiff(file,
                           width  = as.numeric(input$image_width),
                           height = as.numeric(input$image_height),
                           units  = "in", res = as.numeric(input$image_res)),
             "svg"  = svg(file,
                          width  = as.numeric(input$image_width),
                          height = as.numeric(input$image_height))
      )
      on.exit(dev.off(), add = TRUE)

      # draw!
      p <- built_plot()         # ggplot object (or similar)
      print(p)                  # crucial for ggplot2
    }
  )

  output$download_ploidy_table_hmm <- downloadHandler(
    filename = function() {
      paste0("ploidy_table_hmm_", Sys.Date(), ".tsv")
    },
    content = function(file) {
      req(ploidies_hmm())
      write.table(ploidies_hmm()[[1]], file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )

}

## To be copied in the UI
# mod_qploidy_ui("qploidy_1")

## To be copied in the server
# mod_qploidy_server("qploidy_1")
