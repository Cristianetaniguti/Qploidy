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
             box(width= 12, solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,  status="primary", title = "Upload allele intensities or counts", label = tags$b("Upload allele intensities or counts"),

                 box(width= 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,  status="primary", title = "Axiom array files", label = tags$b("Axiom array files"),
                     p("Axiom array summary file:"),
                     fileInput(ns("load_summary"), label = "File input"),
                     p("File with sample names:"),
                     fileInput(ns("load_ind_names"), label = "File input"),
                     p("File with genomic position of each marker:"),
                     fileInput(ns("load_geno_pos"), label = "File input")
                 ),
                 box(width= 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,  status="primary", title = "Illumina array intensities file", label = tags$b("Illumina array intensities file"),
                     p("Illumina array intensities file:"),
                     fileInput(ns("load_illumina"), label = "File input")
                 ),
                 box(width= 12, solidHeader = FALSE, collapsible = TRUE, collapsed = TRUE,  status="primary", title = "VCF file with read counts", label = tags$b("VCF file with read counts"),
                     p("VCF file with read counts:"),
                     fileInput(ns("load_vcf"), label = "File input")
                 ),
                 box(width= 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,  status="primary", title = "Or choose an example dataset", label = tags$b("Or choose an example dataset"),
                     radioButtons(ns("example_data"), label = "Choose example data set",
                                  choices = c("Example data" = "example_data"),
                                  selected = "example_data")
                 )
             )
      ),

      column(12,
             box(width= 12, solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,  status="primary", title = "Array data clusterization", label = tags$b("Array data clusterization"),
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
                 hr(),
                 column(12,
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
                        numericInput(ns("n.cores"), label = "Number of cores to be used", value = 1),
                        numericInput(ns("ploidy"), label = "Interpolation model ploidy", value = 4),
                        downloadButton(ns("down_baf"), "Download BAF"),
                        downloadButton(ns("down_z"), "Download z scores"),
                        downloadButton(ns("data_interpolation"), "Download interpolated data"),
                        hr()
                 )
             )
      ),
      box(width= 12, solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,  status="primary", title = "Interpolation results", label = tags$b("Interpolation results"),

          box(width= 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,  status="primary", title = "Upload interpolated data", label = tags$b("Upload interpolated data"),
              fileInput(ns("data_interpolation_in"), label = "Upload interpolated data:"),
          ),
          box(width= 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,  status="primary", title = "Filtered markers", label = tags$b("Filtered markers"),

              p("Filtered markers:"),
              DT::dataTableOutput(ns("filtered_markers")), br(), hr(),
          ),
          box(width= 12, solidHeader = FALSE, collapsible = TRUE, collapsed = FALSE,  status="primary", title = "Single marker interpolation plot", label = tags$b("Single marker interpolation plot"),

              pickerInput(ns("int_marker"),
                          label = "Select marker",
                          choices = "This will be updated",
                          selected = "This will be updated",
                          options = pickerOptions(
                            size = 12,
                            `selected-text-format` = "count > 3",
                            `live-search`=TRUE,
                            actionsBox = TRUE,
                            dropupAuto = FALSE
                          ),
                          multiple = FALSE, width = '100%'), hr(),
              plotOutput(ns("interpolation"))
          )
      )
    )
  )
}

#' interpolation Server Functions
#'
#'
#' @import parallel
#' @import vroom
#'
#' @noRd
mod_interpolation_server <- function(id){
  moduleServer(id, function(input, output, session){

    ns <- session$ns

    input_summary <- reactive({
      if(!is.null(input$load_summary)){
        summary <- vroom(input$load_summary$datapath, show_col_types = FALSE)
        ind.names <- vroom(input$load_ind_names$datapath, show_col_types = FALSE)

        cleaned_summary <- clean_summary(summary_df = summary)

        R_theta <- get_R_theta(cleaned_summary, ind.names)
        R_all <- R_theta[[1]]
        theta_all <- R_theta[[2]]
        fitpoly_input <- summary_to_fitpoly(R_all, theta_all)

        result <- list(fitpoly_input=fitpoly_input, R_all=R_all, theta_all=theta_all)
        return(result)
      } else {
        return(NULL)
      }
    })

    loadExample <- reactive({
      withProgress(message = 'Working:', value = 0, {
        if(is.null(input_summary())){
          if(input$example_data == "example_data"){
            incProgress(0.5, detail = paste("Loading example data..."))
            result <- readRDS(system.file("summary_result_example.rds", package = "Qploidy"))
            return(result)
          }
        } else {
          return(NULL)
        }
      })
    })

    summary <- reactive({
      if(is.null(input_summary())){
        summary_df <- loadExample()
      } else {
        summary_df <- input_summary()
      }
      return(summary_df)
    })

    observe({
      choices_names <- as.list(unique(summary()[[1]]$SampleName))
      names(choices_names) <- as.list(unique(summary()[[1]]$SampleName))

      updatePickerInput(session, "refs",
                        label = "Select reference samples",
                        choices = choices_names,
                        selected=unlist(choices_names)[1:3])
    })

    refs_fitpoly <- reactive({
      withProgress(message = 'Working:', value = 0, {
        incProgress(0.5, detail = paste("Filtering scores..."))

        refs_fitpoly <- summary()[[1]]
        refs_fitpoly <- refs_fitpoly %>% filter(.data$SampleName %in% input$refs)
        return(refs_fitpoly)
      })
    })

    output$down_fitpoly <- downloadHandler(
      filename = function() {
        # Use the selected dataset as the suggested file name
        temp <- sample(1:1000, 1)
        paste0("reference_samples_",temp, ".txt")
      },
      content = function(file) {
        # Write the dataset to the `file` that will be downloaded
        write.table(refs_fitpoly(), file, row.names = FALSE)
      }
    )

    fitpoly_scores <- reactive({
      withProgress(message = 'Working:', value = 0, {

        if(is.null(input$load_scores)){
          scores <- NULL
        } else {
          incProgress(0.5, detail = paste("Loading scores data..."))

          scores <- vroom(input$load_scores$datapath, show_col_types = FALSE)
        }
      })

      scores
    })

    lst_interpolation <- reactive({
      withProgress(message = 'Working:', value = 0, {
        incProgress(0.5, detail = paste("Preparing for clusterization..."))
        if(!is.null(fitpoly_scores())){
          n.na <- fitpoly_scores() %>% group_by(.data$MarkerName) %>% summarize(n.na = (sum(is.na(.data$geno))/length(.data$geno))*100)
          rm.mks <- n.na$MarkerName[which(n.na$n.na > 25)]
          missing.data <- length(rm.mks)
          scores_filt <- fitpoly_scores()[-which(fitpoly_scores()$MarkerName %in% rm.mks),]

          keep.mks <- match(paste0(scores_filt$MarkerName, "_",scores_filt$SampleName),
                            paste0(summary()[[1]]$MarkerName, "_", summary()[[1]]$SampleName))

          fitpoly_input_filt <- summary()[[1]][keep.mks,]
          theta <- fitpoly_input_filt$ratio
          R <- fitpoly_input_filt$R

          rm.na <- which(is.na(scores_filt$geno))
          if(length(rm.na) > 0){
            R[rm.na] <- NA
            theta[rm.na] <- NA
          }

          data_interpolation <- data.frame(mks = fitpoly_input_filt$MarkerName,
                                           ind = fitpoly_input_filt$SampleName,
                                           R = R,
                                           theta = theta,
                                           geno = scores_filt$geno)
          lst_inter <- split(data_interpolation, data_interpolation$mks)

          return(list(data_interpolation, missing.data))
        } else {
          return(NULL)
        }
      })
    })

    clusters <- reactive({
      withProgress(message = 'Working:', value = 0, {
        incProgress(0.5, detail = paste("Defining clusters centroids..."))

        # Generate clusters
        if(!is.null(fitpoly_scores())){
          ploidy_r <- input$ploidy
          clust <- makeCluster(input$n.cores)
          clusterExport(clust, c("par_fitpoly_interpolation"))
          clusters <- parLapply(clust, lst_interpolation()[[1]], function(x) {
            par_fitpoly_interpolation(x, ploidy= ploidy_r)
          })
          stopCluster(clust)
          return(clusters)
        } else {
          return(NULL)
        }
      })
    })

    filters <- reactive({
      withProgress(message = 'Working:', value = 0, {
        incProgress(0.5, detail = paste("Counting filtered markers..."))

        if(!is.null(clusters())){
          # Filter by number of clusters
          rm.mks <- sapply(clusters(), function(x) is.null(x$mod))
          wrong_n_clusters <- sum(unlist(rm.mks))
          clusters_filt <- clusters()[-which(rm.mks)]

          # Filtered markers table
          filtered.markers <- data.frame(nrow(summary()[[2]]),
                                         lst_interpolation()[[2]],
                                         wrong_n_clusters,
                                         length(clusters_filt))

          colnames(filtered.markers) <- c("n.mk.start", "missing.data", "wrong_n_clusters", "n.mk.selected")
          return(list(clusters_filt, filtered.markers, rm.mks))
        } else if(!is.null(input$filters_table)){
          filtered.markers <- vroom(input$filter_table)
          return(list(NULL, filtered.markers, NULL))
        } else {
          filtered.markers <- vroom(system.file("filtered.markers.txt", package = "Qploidy"))
          return(list(NULL, filtered.markers, NULL))
        }
      })
    })

    data_interpolation_disc.mks <- reactive({
      withProgress(message = 'Working:', value = 0, {
        incProgress(0.5, detail = paste("Counting filtered markers..."))

        if(is.null(input$data_interpolation_in) & is.null(lst_interpolation())){
          data_interpolation_disc.mks <- readRDS(system.file("data_interpolation.rds", package = "Qploidy"))
        } else if(is.null(input$data_interpolation_in)) {
          # Add discarted
          rm.mks <- sapply(clusters, function(x) is.null(x$mod))
          mks <- names(clusters)
          mks[which(rm.mks)] <- paste(mks[which(rm.mks)], "(discarded)")
          data_interpolation_disc.mks <- do.call(rbind, lst_interpolation()[[1]])
          data_interpolation_disc.mks$mks.disc <- mks[match(data_interpolation_disc.mks$mks, names(clusters()))]
        } else {
          data_interpolation_disc.mks <- vroom(input$data_interpolation_in$datapath)
        }
        return(data_interpolation_disc.mks)
      })
    })

    output$data_interpolation <- downloadHandler(
      filename = function() {
        # Use the selected dataset as the suggested file name
        temp <- sample(1:1000, 1)
        paste0("data_interpolation_",temp, ".txt")
      },
      content = function(file) {
        vroom_write(data_interpolation_disc.mks(), file = file)
      }
    )

    observe({
      withProgress(message = 'Working:', value = 0, {
        incProgress(0.5, detail = paste("Counting filtered markers..."))

        choices <- as.list(unique(data_interpolation_disc.mks()$mks))
        names(choices) <- unique(data_interpolation_disc.mks()$mks.disc)

        updatePickerInput(session, "int_marker",
                          label = "Selected marker",
                          choices = choices,
                          selected=unlist(choices)[1])
      })

    })

    output$interpolation <- renderPlot({
      if(input$int_marker != "This will be updated"){
        lst_interpolation <- split(data_interpolation_disc.mks(), data_interpolation_disc.mks()$mks)
        plot_one_marker(lst_interpolation[[input$int_marker]], ploidy = input$ploidy)
      } else NULL
    })

    output$filtered_markers <- DT::renderDataTable(server = FALSE, {
      DT::datatable(filters()[[2]], extensions = 'Buttons',
                    options = list(
                      dom = 'Bfrtlp',
                      buttons = c('copy', 'csv', 'excel', 'pdf')
                    ),
                    class = "display")
    })

    logR_BAF <- reactive({
      withProgress(message = 'Working:', value = 0, {
        incProgress(0.5, detail = paste("Interpolating BAF..."))

        if(!is.null(filters()[[1]])){
          keep.mks <- names(filters()[[1]])
          # Getting logR for entire dataset
          R_filt <- summary()[[2]][match(keep.mks, summary()[[2]]$MarkerName),]
          theta_filt <- summary()[[3]][match(keep.mks, summary()[[3]]$MarkerName),]

          par <- rep(1:input$n.cores, each=round((nrow(R_filt)/input$n.cores)+1,0))[1:nrow(R_filt)]

          par_R <- split.data.frame(R_filt[,-1], par)
          par_theta <- split.data.frame(theta_filt[,-1], par)
          par_clusters_filt <- split(filters()[[1]], par)

          par_all <- list()
          for(i in 1:input$n.cores){
            par_all[[i]] <- list()
            par_all[[i]][[1]] <- par_R[[i]]
            par_all[[i]][[2]] <- par_theta[[i]]
            par_all[[i]][[3]] <- par_clusters_filt[[i]]
          }

          clust <- makeCluster(input$n.cores)
          ploidy_r <- input$ploidy

          # Get BAF
          clust <- makeCluster(input$n.cores)
          clusterExport(clust, c("get_baf", "get_baf_par"))
          bafs_diplo <- parLapply(clust, par_all, function(x) {
            get_baf_par(x, ploidy = ploidy_r)
          })
          stopCluster(clust)

          bafs_diplo_lt <- unlist(bafs_diplo, recursive = F)
          bafs_diplo_m <- do.call(rbind, bafs_diplo_lt)
          rownames(bafs_diplo_m) <- keep.mks
          bafs_diplo_df <- as.data.frame(bafs_diplo_m)
          bafs_diplo_df <- cbind(mks=rownames(bafs_diplo_df), bafs_diplo_df)

          return(bafs_diplo_df)
        } else {
          return(NULL)
        }
      })
    })

    output$down_baf <- downloadHandler(
      filename = function() {
        # Use the selected dataset as the suggested file name
        temp <- sample(1:1000, 1)
        paste0("BAF_",temp, ".txt")
      },
      content = function(file) {
        vroom_write(logR_BAF(), file = file)
      }
    )
  })
}
## To be copied in the UI
# mod_interpolation_ui("interpolation_1")

## To be copied in the server
# mod_interpolation_server("interpolation_1")
