#' single UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_single_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
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
    mainPanel(
      tabsetPanel(
        tabPanel("Tables"),
        tabPanel("Graphics",
                 column(12,
                        p("Segmented logR plot"),
                        plotOutput(ns("plot_logR")), br(),
                 ),
                 column(12,
                        p("BAF plot"),
                        plotOutput(ns("plot_lines")), br(),
                        plotOutput(ns("plot_hist")))
        )
      )
    )
  )
}

#' single Server Functions
#'
#' @noRd
mod_single_server <- function(input, output, session, data, parent_session){
    ns <- session$ns

    observe({
      choices_names <- as.list(unique(colnames(data()[[1]])[-c(1:3)]))
      names(choices_names) <- unique(colnames(data()[[1]])[-c(1:3)])

      updatePickerInput(session, "graphics",
                        label = "Select sample for the graphic",
                        choices = choices_names,
                        selected=unlist(choices_names)[1])

      updatePickerInput(session, "graphics_parents",
                        label = "Select parents samples for the graphic",
                        choices = choices_names,
                        selected=NULL)
    })

    # input <- list()
    # input$graphics <- "Old_Blush"
    #
    # input$graphics <- "Banana"
    # input$graphics <- "10038-N001"
    # input$graphics <- "George_Vancouver"
    # input$graphics <- "16405_N122" # tetraploid Chr1 aneuploid 5
    # input$area_single <- 0.75
    # input$ploidy <- 2

    graphics_baf <- eventReactive(input$run_individual,{
      #data_sample <- data[[2]][,c(2,3,which(colnames(data[[2]]) %in% c(input$graphics)))]
      data_sample <- data()[[2]][,c(2,3,which(colnames(data()[[2]]) %in% c(input$graphics)))]
      colnames(data_sample)[3] <- "sample"

      # sample.chr <- split(data_sample$sample, data_sample$Chr)
      # modes <- lapply(sample.chr, function(x) {
      #   modes <- Modes(x)
      #   modes.observed <- sort(modes$modes)
      #   modes.expected <- seq(0,1,1/(input$ploidy))
      #   joint <- c(modes.expected, modes.observed)
      #   joint.name <- c(rep("modes expected", length(modes.expected)), rep("modes observed", length(modes.observed)))
      #   alpha <- c(rep(1, length(modes.expected)), rep(0.5, length(modes.observed)))
      #
      #   df <- data.frame(value=joint, name = joint.name, alpha)
      #   return(df)
      # })
      #
      # modes.df <- data.frame()
      # for(i in 1:length(modes)){
      #   modes.temp <- cbind(modes[[i]], Chr = names(modes)[i])
      #   modes.df <- rbind(modes.df, modes.temp)
      # }

      ymin <- seq(0, 1, 1/input$ploidy) - (input$area_single/(input$ploidy*2))
      ymax <- seq(0, 1, 1/input$ploidy) + (input$area_single/(input$ploidy*2))

      ymin[which(ymin < 0)] <- 0
      ymax[which(ymax > 1)] <- 1
      rets <- data.frame(ymin, ymax, xmax = Inf, xmin = -Inf)

      idx_tot <- FALSE
      idx <- list()
      for(i in 1:nrow(rets)){
        idx <- data_sample$sample >= rets$ymin[i] & data_sample$sample <= rets$ymax[i]
        idx_tot <- idx_tot | idx
      }

      data_sample$color <- NA
      data_sample$color[which(!idx_tot)] <- "red"
      data_sample$color[which(idx_tot)] <- "black"

      p_baf <- data_sample %>% ggplot(aes(x=Position, y=sample)) +
        {if(input$colors) geom_point(aes(color = color), alpha=0.7, size=input$dot.size) else geom_point(alpha=0.7, size=input$dot.size)} +
        scale_color_manual(values = c("red", "black")) +
        {if(input$add_lines) geom_rect(data = rets, inherit.aes=FALSE,
                                       aes(ymin=ymin,ymax=ymax,
                                           xmax = xmax, xmin = xmin,
                                           alpha=0.001, color = "red"))} +
        ylab("BAF") +
        facet_grid(~ Chr, scales = "free_x") + theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              legend.position="none")

      data_sample$color[which(data_sample$color == "black")] <- "inside area"
      data_sample$color[which(data_sample$color == "red")] <- "outside area"

      modes.df2 <- list()
      for(i in 1:nrow(rets)) {
        split_data <- split(data_sample$sample, data_sample$Chr)
        x <- split_data[[1]]
        modes.df2[[i]] <- matrix(unlist(sapply(split_data, function(x){
          estimated <- median(x[which(x > rets$ymin[i] & x < rets$ymax[i])])
          expected <- seq(0, 1, 1/input$ploidy)[i]
          return(data.frame(estimated, expected))
        })), nrow = 2)
        colnames(modes.df2[[i]]) <- names(split_data)
        rownames(modes.df2[[i]]) <- c("modes estimated", "modes expected")
      }

      modes.df2 <- lapply(modes.df2, function(x) {
        y <- t(x)
        y <- data.frame(Chr=rownames(y), y)
        pivot_longer(y, cols = 2:3)
      })

      modes.df2 <- do.call(rbind, modes.df2)
      modes.df2$alpha <- modes.df2$name
      modes.df2$alpha <- gsub("modes.expected", 1, modes.df2$alpha)
      modes.df2$alpha <- as.numeric(gsub("modes.estimated", 0.5, modes.df2$alpha))

      p_hist <- data_sample %>% ggplot(aes(x=sample)) +
        {if(input$colors) geom_histogram(aes(fill = color)) else  geom_histogram()} +
        scale_x_continuous(breaks = round(seq(0, 1, 1/input$ploidy),2)) +
        {if(input$add_lines) geom_vline(data = modes.df2, aes(xintercept= value,
                                                              color = name,
                                                              linetype = name,
                                                              alpha = alpha),
                                        linewidth = 0.8)}+
        {if(input$colors) scale_fill_manual(values = c("red", "black"))} +
        scale_color_manual(values = c("blue", "purple")) +
        scale_linetype_manual(values = c("dashed", "solid"), guide="none") +
        scale_alpha(range = c(0.7, 1), guide="none") +
        facet_grid(~ Chr, scales = "free_x") + theme_bw() +  xlab("BAF") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              legend.position="bottom") +
        labs(color= "Peaks", fill = "Area")

      p_hist_all <- data_sample %>% ggplot(aes(x=sample)) +
        {if(input$colors) geom_histogram(aes(fill = color)) else  geom_histogram()} +
        scale_x_continuous(breaks = round(seq(0, 1, 1/input$ploidy),2)) +
        {if(input$colors) scale_fill_manual(values = c("red", "black"))} +
        scale_color_manual(values = c("blue", "purple")) +
        scale_linetype_manual(values = c("dashed", "solid"), guide="none") +
        scale_alpha(range = c(0.7, 1), guide="none") +
        theme_bw() +  xlab("BAF") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              legend.position="bottom") +
        labs(color= "Peaks", fill = "Area")

      p_hist <- ggarrange(p_hist_all, p_hist, widths = c(1,max(data_sample$Chr)-2))

      list(p_baf, p_hist)
    })


    output$plot_lines <- renderPlot({
      graphics_baf()[[1]]
    })

    output$plot_hist <- renderPlot({
      graphics_baf()[[2]]
    })

    graphics_logR <- eventReactive(input$run_individual,{
      # logR_sample <- data[[1]][,c(2,3,which(colnames(data[[1]]) %in% c(input$graphics)))]
      logR_sample <- data()[[1]][,c(2,3,which(colnames(data()[[1]]) %in% c(input$graphics)))]
      colnames(logR_sample)[3] <- "sample"
      # Segmentation of logR
      logR_list <- split(logR_sample, logR_sample[,1])
      logR_list <- lapply(logR_list, function(x) x[order(x[,2]$Position),])
      segs_all <- lapply(logR_list, function(x) {
        if(length(which(is.na(x$sample))) > 0) {
          x1 <- x$Position[-which(is.na(x$sample))]
          y1 <- x$sample[-which(is.na(x$sample))]
        } else {
          x1 <- x$Position
          y1 <- x$sample
        }

        dpseg(x= x1,
              y= y1,
              jumps=FALSE,
              P=0.001,
              type="var",
              store.matrix=TRUE)
      })

      logR_list_int <- logR_list
      for(j in 1:length(segs_all)){
        intercept <- vector()
        for(i in 1:length(segs_all[[j]]$segments$x1)){
          idx_start <- which(logR_list[[j]]$Position == segs_all[[j]]$segments$x1[i])
          idx_end <- which(logR_list[[j]]$Position == segs_all[[j]]$segments$x2[i])
          intercept[idx_start[1]:idx_end[length(idx_end)]] <- segs_all[[j]]$segments$intercept[i]
        }
        if(length(intercept) != nrow(logR_list[[j]])) intercept[length(intercept):nrow(logR_list[[j]])] <- NA # bugfix
        logR_list_int[[j]] <- cbind(logR_list[[j]],intercept)
      }

      logR_df_int <- do.call(rbind, logR_list_int)

      #lim <- max(c(abs(max(logR_df_int$sample, na.rm = T)), abs(min(logR_df_int$sample, na.rm = T))))
      segs_plot <- logR_df_int %>% ggplot(aes(x=Position)) +
        geom_point(aes(y = sample), alpha = 0.3) +
        geom_line(aes(y=intercept), color = "red") + ylab("logR") +
        facet_grid(.~Chr , scales = "free_x") + theme_bw() + ylim(-1,1)+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

      segs_plot
    })

    output$plot_logR <- renderPlot({
      graphics_logR()
    })
}

## To be copied in the UI
# mod_single_ui("single_1")

## To be copied in the server
# mod_single_server("single_1")
