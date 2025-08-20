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
             box(
               title = "Plots", status = "info", solidHeader = FALSE, width = 12,
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
                                          "Non-standardized Ratio Histogram" = "Ratio_hist_overall"),
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
                            actionButton(ns("build_plot"), "Build plot"),

                     )
                   )
               ), hr(),

               plotOutput(ns("plot"), height = "500px"),

               div(style="display:inline-block; float:left",dropdownButton(
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
               ))
             ),
             box(
               title = "Ploidy Estimation", status = "info", solidHeader = FALSE, width = 12,
               "Check the plots above before obtaining the estimated ploidies for all samples.",
               "If the plots look good, adjust the following paramenters and click the 'Estimate Ploidies' button below to run the estimation.", br(),
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
               br(),
               actionButton(ns("est_ploidy"), "Estimate Ploidies"), hr(), br(),
               DTOutput(ns("ploidy_table")),br(),
               downloadButton(ns("download_ploidy_table"), "Download Ploidy Table")
             )
      )
    )
  )
}

#' qploidy Server Functions
#'
#' @importFrom DT renderDT
#' @importFrom bs4Dash updatebs4TabItems updateBox
#' @importFrom shiny updateTabsetPanel
#' @import dplyr
#' @importFrom shinyalert shinyalert
#'
#' @noRd
mod_qploidy_server <- function(input, output, session, parent_session){

  ns <- session$ns

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

  #UI popup window for input
  observeEvent(input$advanced_options, {
    showModal(modalDialog(
      title = "Advanced Options",
      numericInput(ns("miss"), "Remove SNPs with >= % missing data", min = 1, max = 100,value = 90),
      numericInput(ns("prob"), "Remove genotypes with <= probability", min = 0, value = 0.8),
      numericInput(ns("n.clusters"), "Remove genotypes with <= # of dosage clusters",
                   min = 2, value = {
                     ploidy <- c(input$ref_ploidy, input$common_ploidy)
                     ploidy <- ploidy[-is.na(ploidy)]
                     as.numeric(ploidy) + 1
                   }),
      footer = tagList(
        modalButton("Close"),
        actionButton(ns("save_advanced_options"), "Save")
      )
    ))
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
    req(input$input_file)
    req(input$file_type)

    # Reset progress bar
    updateProgressBar(session = session, id = "pb_qploidy", value = 0, total = 100)

    # Read input file based on type
    if (input$file_type == "VCF") {
      # Add VCF format checks

      input_data <- qploidy_read_vcf(input$input_file$datapath)

    } else if (input$file_type == "Axiom Array") {
      # Add axiom array summary file format checks

      input_data <- read_axiom(input$input_file$datapath)

    } else if(input$file_type == "Illumina Array") {
      # Add Illumina array summary file format checks

      input_data <- read_illumina_array(input$input_file$datapath)

    } else if (input$file_type == "Qploidy Standardized Dataset") {
      # Add Qploidy standardized dataset format checks
      input_data <- "skip"
      data_standardized <- read_qploidy_standardization(input$input_file$datapath)
    }

    # Update progress bar
    updateProgressBar(session = session, id = "pb_qploidy", value = 20)

    if(input$file_type != "Qploidy Standardized Dataset"){

      if(input$file_type == "VCF"){
        genos <- qploidy_read_vcf(input$input_file$datapath, geno = TRUE)
        genos.pos <- qploidy_read_vcf(input$input_file$datapath, geno.pos = TRUE)

        if(input$known_ploidy){
          reference_samples <- read.csv(input$reference_samples$datapath, header = TRUE, stringsAsFactors = FALSE)$Sample_ID
          genos <- genos[which(genos$SampleName %in% reference_samples), ]
        }

      } else {
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

      if(max(genos$geno) != ploidy){
        shinyalert(
          title = "Ploidy Mismatch",
          text = "Genotype on File doesn't match Ploidy informed",
          size = "s",
          closeOnEsc = TRUE,
          closeOnClickOutside = FALSE,
          html = TRUE,
          type = "error",
          showConfirmButton = TRUE,
          confirmButtonText = "OK",
          confirmButtonCol = "#004192",
          showCancelButton = FALSE,
          animation = TRUE
        )
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

    valueBox(
      value = tagList(
        div(style="font-size:1.2em; line-height:1.1;", span("Total markers: "), formatC(data_standardized()$filters[1], big.mark=",")),
        div(style="font-size:0.95em; opacity:0.9;", "Markers removed due:"),
        div(style="font-size:0.95em; opacity:0.9;", span("Low geno prob: "), strong(data_standardized()$filters[2])),
        div(style="font-size:0.95em; opacity:0.9;", span("High missing data: "), strong(data_standardized()$filters[3])),
        div(style="font-size:0.95em; opacity:0.9;", span("Low # of dosage clusters: "), strong(data_standardized()$filters[4])),
        div(style="font-size:0.95em; opacity:0.9;", span("Missing genomic position: "), strong(data_standardized()$filters[5])),
        div(style="font-size:1.2em; line-height:1.1;", span("Remaining markers: "), strong(data_standardized()$filters[6]))
      ),
      subtitle = "",
      icon = icon("dna"),
      color = if(data_standardized()$filters[6] != 0)
        if (data_standardized()$filters[6]/data_standardized()$filters[1] > 0.25) "success" else "danger"
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

    chrs <- unique(data_standardized()$data$Chr)
    updateVirtualSelect(
      session = session,
      inputId = "sele_chr",
      choices = chrs,
      selected = if(length(chrs) <= 5) chrs else chrs[1:5]
    )
  })

  ploidies <- eventReactive(input$est_ploidy, {
    req(data_standardized())

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
  })

  output$ploidy_table <- renderDT({
    req(ploidies())
    ploidies()$ploidy
  }, options = list(pageLength = 10, scrollX = TRUE))

  output$download_ploidy_table <- downloadHandler(
    filename = function() {
      paste0("ploidy_estimation_", Sys.Date(), ".tsv")
    },
    content = function(file) {
      req(ploidies())
      write.table(ploidies()$ploidy, file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )

  built_plot <- eventReactive(input$build_plot, {
    req(data_standardized())
    req(input$sample)
    req(input$plots)
    req(input$sele_chr)

    # Update progress bar
    updateProgressBar(session = session, id = "pb_qploidy", value = 0, total = 100)

    chrs <- unique(data_standardized()$data$Chr)
    chrs_idx <- which(chrs %in% input$sele_chr)

    if(!is.null(input$centromeres_file$datapath)){
      centromeres <- read.csv(input$centromeres_file$datapath, header = TRUE, stringsAsFactors = FALSE)
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
      chr = chrs_idx,
      ploidy = input$ploidy,
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

}

## To be copied in the UI
# mod_qploidy_ui("qploidy_1")

## To be copied in the server
# mod_qploidy_server("qploidy_1")
