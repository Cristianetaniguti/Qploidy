#' interpolation UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_interpolation_ui <- function(id){
  ns <- NS(id)
  tagList(
    fixedRow(
      column(12,
             p("Upload allele intensities or counts"),
             box(width= 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,  status="primary", title = "Axiom array files", label = tags$b("Axiom array files"),
                 p("Axiom array summary file:"),
                 fileInput(ns("load_summary"), label = "File input"),
                 p("File with sample names:"),
                 fileInput(ns("load_ind_names"), label = "File input"),
                 p("File with genomic position of each marker:"),
                 fileInput(ns("load_geno_pos"), label = "File input")
             ),
             p("Illumina array intensities file:"),
             fileInput(ns("load_illumina"), label = "File input"),
             p("VCF file with read counts:"),
             fileInput(ns("load_vcf"), label = "File input"),
             p("Or choose an example dataset:"),
             radioButtons(ns("example_data"), label = "Choose example data set",
                          choices = c("Roses Texas" = "roses_texas",
                                      "Roses France" = "roses_france",
                                      "Potatoes Texas" = "potatoes"),
                          selected = "roses_texas")
      ),

      column(12,
             p("Choose samples to be used as reference for interpolation (minimum 50).
        We suggest they to be the most known abundant ploidy in your population.
        Example: If you 100 tetraploids, 30 diploids, and 1000 unknown,
               select the tetraploids and use `interpolation model ploidy = 4`."),
             pickerInput(ns("refs"),
                         label = "Select reference samples",
                         choices = "This will be updated after samples are selected above",
                         selected = "This will be updated after samples are selected above",
                         options = pickerOptions(
                           size = 12,
                           `selected-text-format` = "count > 3",
                           `live-search`=TRUE,
                           actionsBox = TRUE,
                           dropupAuto = FALSE
                         ),
                         multiple = TRUE, width = '100%'),
             numericInput(ns("ploidy"), label = "Interpolation model ploidy", value = 4),
             hr()
      ),
      column(12,
             h3("Array data clusterization"), br(),
             p("If you uploaded array data, the clusterization will be made using fitpoly.
        There are many parameters in fitpoly that can be adjusted according with your population structure.
               Download the fitpoly input data, run the model according to your data set population structure.
               fitPoly will take some time to run depending on the number of individuals, markers, and cores.
               Upload your `_scores.dat` results below."),

             downloadButton(ns("down_fitpoly"), "Download fitpoly input"),
             hr(),
             p("Example code for Hardy-Weinberg model:"),
             code("saveMarkerModels(ploidy= 4, \n
                                    data=refs_fitpoly_inputs, \n
                                    p.threshold=0.9, \n
                                    filePrefix= 'fitpoly_out_', \n
                                    ncores=30)"),
             hr(),
             fileInput(ns("load_scores"), label = "Upload fitpoly scores output file:"),
             hr(),
      ), br(),
    ), hr(),
    mainPanel(
      tabsetPanel(
        tabPanel("Tables"),
        tabPanel("Graphics")
      )
    )
  )
}

#' interpolation Server Functions
#'
#' @noRd
mod_interpolation_server <- function(id){
  moduleServer(id, function(input, output, session){

    ns <- session$ns

    input_summary <- reactive({
      if(!is.null(input$load_summary)){
        summary <- vroom(input$load_summary$datapath)
        cleaned_summary <- clean_summary(summary_df = summary)

        ind.names <- vroom(input$load_ind_names$datapath)
        geno.pos <- vroom(input$load_geno_pos$datapath)

        fitpoly_input <- summary_to_fitpoly(cleaned_summary = cleaned_summary, ind.names, geno.pos)

        fitpoly_input
      } else {
        NULL
      }
    })

    loadExample <- reactive({
      if(is.null(input_summary())){
        if(input$example_data == "roses_texas"){
          summary <- vroom(system.file("fitpoly_input.txt", package = "Qploidy"))
        } else if(input$example_data == "roses_france"){
          summary <- vroom(system.file("fitpoly_input.txt", package = "Qploidy"))
        } else if(input$example_data == "potatoes") {
          summary <- vroom(system.file("fitpoly_input.txt", package = "Qploidy"))
        }
        summary
      } else NULL
    })

    summary <- reactive({
      if(!is.null(input_summary())) return(input_summary()) else return(loadExample())
    })

    observe({
      choices_names <- as.list(unique(summary()$SampleName))
      names(choices_names) <- as.list(unique(summary()$SampleName))

      updatePickerInput(session, "refs",
                        label = "Select reference samples",
                        choices = choices_names,
                        selected=unlist(choices_names)[1:3])
    })

    refs_fitpoly <- reactive({
        df <- summary()
        df <- df %>% filter(SampleName %in% input$refs)
        df
    })

    output$down_fitpoly <- downloadHandler(
      filename = function() {
        # Use the selected dataset as the suggested file name
        paste0("reference_samples.txt")
      },
      content = function(file) {
        # Write the dataset to the `file` that will be downloaded
        write.table(refs_fitpoly(), file)
      }
    )

  })
}

## To be copied in the UI
# mod_interpolation_ui("interpolation_1")

## To be copied in the server
# mod_interpolation_server("interpolation_1")
