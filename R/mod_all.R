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
        tabPanel(ns("Tables"),
                 p("The table include the ploidy estimated by chromosome and overall using the area peak method.
            The measures of quality are the difference between first and second place in the area method and the correlation coefficient (Pearson) between estimated and expected peak position."),
                 column(12,
                        p("Download ploidy estimations by chromosome table:"),
                        downloadButton(ns('result.ploidy_download'), "Download"), br(), hr(),
                        DT::dataTableOutput(ns("result.ploidy")), br(),

                        p("Download Proportion of dots inside selected area:"),
                        downloadButton(ns('dots.int_tot_mt_download'), "Download"), br(), hr(),
                        DT::dataTableOutput(ns("dots.int_tot_mt")), br(),

                        p("Download Difference between first and second place in area method table:"),
                        downloadButton(ns('diff.first.second_download'), "Download"), br(), hr(),
                        DT::dataTableOutput(ns("diff.first.second")), br(),

                        p("Download standard deviation inside area table:"),
                        downloadButton(ns('sd_tot_mt_download'), "Download"), br(), hr(),
                        DT::dataTableOutput(ns("sd_tot_mt")), br(),

                        p("Download Highest correlation table:"),
                        downloadButton(ns('corr_tot_mt_download'), "Download"), br(), hr(),
                        DT::dataTableOutput(ns("corr_tot_mt")), br(),

                        p("Download Modes inside areas table:"),
                        downloadButton(ns('modes_paste_tot_mt_download'), "Download"), br(), hr(),
                        DT::dataTableOutput(ns("modes_paste_tot_mt")), br()
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
mod_all_server <- function(input, output, session, data, parent_session){
    ns <- session$ns

    observe({
      choices_names <- as.list(unique(colnames(data()[[1]])[-c(1:3)]))
      names(choices_names) <- unique(colnames(data()[[1]])[-c(1:3)])

      updatePickerInput(session, "samples_chr",
                        label = "Select sample for analysis by chr",
                        choices = choices_names,
                        selected=unlist(choices_names)[1])

      updatePickerInput(session, "samples",
                        label = "Select samples for overall analysis",
                        choices = choices_names,
                        selected=unlist(choices_names)[1])
    })

    est.ploidy.chr_df <- eventReactive(input$run_overal,{
      data_sample <- data()[[2]][,c(2,3,which(colnames(data()[[2]]) %in% c(input$samples)))]
      #data_sample <- data[[2]][,c(2,3,which(colnames(data[[2]]) %in% c(input$samples)))]
      #data_sample <- data[[2]][,-1]

      data_sample <- data_sample[order(data_sample$Chr, data_sample$Position),]
      by_chr <- split(data_sample, data_sample$Chr)

      freq <- pascalTriangle(input$ploidys[2])
      freq <- freq[-1]
      ploidys <- input$ploidys[1]:input$ploidys[2]
      p.values_tot <- vector()
      dots.int_tot <- corr_tot <- max_sd_tot <- modes_paste_tot <- as.list(rep(NA, length(by_chr)))
      means <- list()
      for(j in 1:length(ploidys)){
        # Area method
        ymin <- seq(0, 1, 1/ploidys[j]) - (input$area/(ploidys[j]*2))
        ymax <- seq(0, 1, 1/ploidys[j]) + (input$area/(ploidys[j]*2))

        ymin[which(ymin < 0)] <- 0
        ymax[which(ymax > 1)] <- 1
        rets <- data.frame(ymin, ymax)

        prop_tot <- modes_tot  <- sd_areas_tot <- as.list(rep(NA, length(by_chr)))
        for(z in 1:length(by_chr)){
          for(i in 1:nrow(rets)){
            prop <- apply(by_chr[[z]][,-c(1,2)], 2, function(x) sum(x >= rets$ymin[i] & x <= rets$ymax[i]))
            modes <- apply(by_chr[[z]][,-c(1,2)], 2, function(x) median(x[x >= rets$ymin[i] & x <= rets$ymax[i]]))
            sd_areas <- apply(by_chr[[z]][,-c(1,2)], 2, function(x) sd(x[x >= rets$ymin[i] & x <= rets$ymax[i]]))

            prop_tot[[z]] <- rbind(prop_tot[[z]], prop)
            modes_tot[[z]] <- rbind(modes_tot[[z]], modes)
            sd_areas_tot[[z]] <- rbind(sd_areas_tot[[z]], sd_areas)
          }
        }

        # remove NAs
        prop_tot <- lapply(prop_tot, function(x) x[-1,])
        for(z in 1:length(prop_tot)){
          if(is.null(ncol(prop_tot[[z]]))) dots.int <- sum(prop_tot[[z]])/dim(by_chr[[z]])[1] else
            dots.int <- apply(prop_tot[[z]], 2, function(x) sum(x)/dim(by_chr[[z]])[1])
          dots.int_tot[[z]] <- rbind(dots.int_tot[[z]], dots.int)

          modes_paste <- apply(modes_tot[[z]], 2, function(x) paste0(x[-1], collapse = "/"))
          modes_paste_tot[[z]] <- rbind(modes_paste_tot[[z]], modes_paste)

          corr <- apply(modes_tot[[z]], 2, function(x) cor(x = x[-1], y = seq(0,1,1/ploidys[j])))
          corr_tot[[z]] <- rbind(corr_tot[[z]], corr)

          max_sd <- apply(sd_areas_tot[[z]], 2, function(x) max(x[-1]))
          max_sd_tot[[z]] <- rbind(max_sd_tot[[z]], max_sd)
        }
      }

      dots.int_tot <- lapply(dots.int_tot, function(x) x[-1,])
      modes_paste_tot <- lapply(modes_paste_tot, function(x) x[-1,])
      corr_tot <- lapply(corr_tot, function(x) x[-1,])
      max_sd_tot <- lapply(max_sd_tot, function(x) x[-1,])

      # Area method
      result.ploidy  <- diff.count <- second <- diff.second <- diff.first.second <- list()
      for(z in 1:length(dots.int_tot)){
        if(is.null(rownames(dots.int_tot[[z]]))) {
          names(dots.int_tot[[z]]) <- ploidys
          result.ploidy[[z]] <- ploidys[order(dots.int_tot[[z]], decreasing =T)][1]

          diff.count[[z]] <- dots.int_tot[[z]][order(dots.int_tot[[z]], decreasing =T)][1]
          second[[z]] <- ploidys[order(dots.int_tot[[z]], decreasing =T)][2]
          diff.second[[z]] <- dots.int_tot[[z]][order(dots.int_tot[[z]], decreasing =T)][2]
        } else {
          rownames(dots.int_tot[[z]]) <- ploidys
          result.ploidy[[z]] <- apply(dots.int_tot[[z]], 2, function(x) ploidys[order(x, decreasing =T)][1])
          diff.count[[z]] <- apply(dots.int_tot[[z]], 2, function(x) x[order(x, decreasing =T)][1])
          second[[z]] <- apply(dots.int_tot[[z]], 2, function(x) ploidys[order(x, decreasing =T)][2])
          diff.second[[z]] <- apply(dots.int_tot[[z]], 2, function(x) x[order(x, decreasing =T)][2])
        }
        diff.first.second[[z]] <- diff.count[[z]] - diff.second[[z]]
        # filt <- which((diff.count[[z]] - diff.second[[z]]) < input$filter)
        # result.ploidy[[z]][filt] <- NA
        # second[[z]][filt] <- NA
      }

      result.ploidy <- do.call(rbind, result.ploidy)
      diff.first.second <- do.call(rbind, diff.first.second)

      if(is.vector(max_sd_tot[[1]])){ # only one individual selected
        max_sd_tot <- lapply(max_sd_tot, as.matrix)
        corr_tot <- lapply(corr_tot, as.matrix)
        modes_paste_tot <- lapply(modes_paste_tot, as.matrix)
        dots.int_tot <- lapply(dots.int_tot, as.matrix)
      } else if(dim(max_sd_tot[[1]])[1] == 1){
        max_sd_tot <- lapply(max_sd_tot, t)
        corr_tot <- lapply(corr_tot, t)
        modes_paste_tot <- lapply(modes_paste_tot, t)
        dots.int_tot <- lapply(dots.int_tot, t)
      }

      sd_tot_mt <- corr_tot_mt <- dots.int_tot_mt <- modes_paste_tot_mt <- matrix(NA, nrow = dim(result.ploidy)[1], ncol = dim(result.ploidy)[2])
      for(i in 1:dim(result.ploidy)[2]){ # Ind
        for(j in 1:dim(result.ploidy)[1]){ # Chr
          sd_tot_mt[j, i] <- max_sd_tot[[j]][which(ploidys == result.ploidy[j,i]),i]
          corr_tot_mt[j,i] <- corr_tot[[j]][which(ploidys == result.ploidy[j,i]),i]
          modes_paste_tot_mt[j,i] <- modes_paste_tot[[j]][which(ploidys == result.ploidy[j,i]),i]
          dots.int_tot_mt[j,i] <- dots.int_tot[[j]][which(ploidys == result.ploidy[j,i]),i]
        }
      }

      result.ploidy <- t(result.ploidy)
      diff.first.second <- t(diff.first.second)
      sd_tot_mt <- t(sd_tot_mt)
      corr_tot_mt <- t(corr_tot_mt)
      modes_paste_tot_mt <- t(modes_paste_tot_mt)
      dots.int_tot_mt <- t(dots.int_tot_mt)

      colnames(result.ploidy) <- colnames(diff.first.second) <- paste0("Chr",names(by_chr))
      colnames(sd_tot_mt) <- colnames(corr_tot_mt) <- colnames(modes_paste_tot_mt) <- colnames(dots.int_tot_mt) <- paste0("Chr",names(by_chr))

      rownames(diff.first.second) <- rownames(result.ploidy)
      rownames(sd_tot_mt) <- rownames(corr_tot_mt) <- rownames(modes_paste_tot_mt) <- rownames(dots.int_tot_mt) <- rownames(result.ploidy)

      est.ploidy.chr_df <- list(result.ploidy,         # Estimated ploidy by area method
                                dots.int_tot_mt,       # Proportion of dots inside selected area
                                diff.first.second,     # Difference between first and second place in area method
                                sd_tot_mt,             # standard deviation inside area
                                corr_tot_mt,           # Highest correlation
                                modes_paste_tot_mt)    # Modes inside areas
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

      idx <- which(apply(est.ploidy.chr_df()[[1]], 1, function(x) length(unique(x)) > 1))
      if(length(idx) > 0){
        ploidies <- est.ploidy.chr_df()[[1]][-idx,]
        ploidies.prop <- est.ploidy.chr_df()[[2]][-idx,]
        ploidies.sd <- est.ploidy.chr_df()[[4]][-idx,]
        ploidies.corr <- est.ploidy.chr_df()[[5]][-idx,]
        ploidies.diff <- est.ploidy.chr_df()[[3]][-idx,]

        # Graphics using only area method and aneuploidy individuals
        aneuploids <- est.ploidy.chr_df()[[1]][idx,]
        aneuploids.prop <- est.ploidy.chr_df()[[2]][idx,]
        aneuploids.sd <- est.ploidy.chr_df()[[4]][idx,]
        aneuploids.corr <- est.ploidy.chr_df()[[5]][idx,]
        aneuploids.diff <- est.ploidy.chr_df()[[3]][idx,]

        if(is.vector(aneuploids)){
          aneuploids <- data.frame(Var1 = rownames(est.ploidy.chr_df()[[1]])[idx], Var2=names(aneuploids),melt(aneuploids, value.name = "ploidy"))
          aneuploids.prop <- data.frame(Var1 = rownames(est.ploidy.chr_df()[[1]])[idx], Var2=names(aneuploids.prop),melt(aneuploids.prop, value.name = "prop"))
          aneuploids.sd <- data.frame(Var1 = rownames(est.ploidy.chr_df()[[1]])[idx], Var2=names(aneuploids.sd),melt(aneuploids.sd, value.name = "sd"))
          aneuploids.corr <- data.frame(Var1 = rownames(est.ploidy.chr_df()[[1]])[idx], Var2=names(aneuploids.corr),melt(aneuploids.corr, value.name = "corr"))
          aneuploids.diff <- data.frame(Var1 = rownames(est.ploidy.chr_df()[[1]])[idx], Var2=names(aneuploids.diff),melt(aneuploids.diff, value.name = "diff"))
        } else {
          aneuploids <- melt(aneuploids, value.name = "ploidy")
          aneuploids.prop <- melt(aneuploids.prop, value.name = "prop")
          aneuploids.sd <- melt(aneuploids.sd, value.name = "sd")
          aneuploids.corr <- melt(aneuploids.corr, value.name = "corr")
          aneuploids.diff <- melt(aneuploids.diff, value.name = "diff")
        }

        # Aneuploid table and samples list
        main_ploidy <- aneuploids %>% group_by(Var1) %>% mutate(main_ploidy = median(ploidy))

        aneuploidy_df <- main_ploidy[-which(main_ploidy$ploidy == main_ploidy$main_ploidy),]
        colnames(aneuploidy_df) <- c("sample", "chromosome", "aneuploidy_ploidy", "main_ploidy")
        aneuploidy_list <- aneuploidy_df %>% group_by(chromosome, aneuploidy_ploidy, main_ploidy) %>% summarise(samples = paste0(sample, collapse = " "))
        aneuploidy_list <- aneuploidy_list[order(aneuploidy_list$main_ploidy, aneuploidy_list$aneuploidy_ploidy, aneuploidy_list$chromosome),]
        aneuploid_text <- paste0("<br/>Main ploidy - ", aneuploidy_list$main_ploidy, "; ploidy - ", aneuploidy_list$aneuploidy_ploidy, "; in chromosome - ", aneuploidy_list$chromosome,":<br/>", aneuploidy_list$samples)

        aneuploidy_all_weird <- unique(c(as.character(aneuploids.corr$Var1)[which(aneuploids.corr$corr < input$filter_corr)],
                                         as.character(aneuploids.diff$Var1)[which(aneuploids.diff$diff < input$filter_diff)]))

        if(length(aneuploidy_all_weird) == 0) aneuploidy_all_weird <- "Set filters"

        # Graphics
        df_n <- melt(table(aneuploids$ploidy), varnames = "ploidy", value.name = "#individuals*#chrom")
        df_n$ploidy <- as.factor(df_n$ploidy)
        p1 <- ggplot(df_n, aes(x = ploidy, y=`#individuals*#chrom`, fill =ploidy)) +
          geom_bar(stat="identity") +
          theme_bw()

        df <- merge(aneuploids, aneuploids.prop, by = c(1,2))
        df <- merge(df, aneuploids.sd, by = c(1,2))
        df <- merge(df, aneuploids.corr, by = c(1,2))
        df <- merge(df, aneuploids.diff, by = c(1,2))

        df <- melt(df, id.vars = c("ploidy", "Var1", "Var2"))

        df$variable <- gsub("sd", "standard deviation", df$variable)
        df$variable <- gsub("prop", "proportion in area", df$variable)
        df$variable <- gsub("corr", "correlation (Pearson)", df$variable)
        df$variable <- gsub("diff", "difference in proportion", df$variable)

        df$ploidy <- as.factor(df$ploidy)
        p2 <- ggplot(df, aes(x = ploidy, y=value, group=ploidy, fill=ploidy)) +
          geom_boxplot() +
          facet_grid(variable~., scales = "free") +
          theme_bw()

        p_aneuploids <- ggarrange(p1, p2, common.legend = T)

        df <- as.data.frame(table(data.frame(chr=aneuploids$Var2, ploidy = aneuploids$ploidy)))
        df <- df %>% group_by(chr, ploidy) %>% summarise(Freq_all = sum(Freq))
        p_aneuploids2 <- ggplot(df, aes(chr,ploidy)) +
          geom_tile(aes(fill = Freq_all)) +
          geom_text(aes(label=Freq_all)) +
          scale_fill_gradient(low="white", high = "red")

      } else {
        ploidies <- est.ploidy.chr_df()[[1]]
        ploidies.prop <- est.ploidy.chr_df()[[2]]
        ploidies.sd <- est.ploidy.chr_df()[[4]]
        ploidies.corr <- est.ploidy.chr_df()[[5]]
        ploidies.diff <- est.ploidy.chr_df()[[3]]

        p_aneuploids <- paste0("No aneuploid samples")
        p_aneuploids2 <- paste0("No aneuploid samples")
        aneuploid_text <- paste0("No aneuploid samples")
        aneuploidy_all_weird <- paste0("No aneuploid samples")
        aneuploidy_df <-  paste0("No aneuploid samples")
      }

      # Graphics using only area method and non-aneuploidy individuals
      ploidy_all_tex <- split(names(ploidies[,1]), ploidies[,1])
      text_ploidy <- vector()
      for(i in 1:length(ploidy_all_tex)){
        line1 <- paste0("<br/>ploidy ",names(ploidy_all_tex)[i], ":<br/>")
        line2 <- paste0(ploidy_all_tex[[i]], collapse = " ")
        text_ploidy <- paste0(text_ploidy, line1, line2)
      }

      df_n <- melt(table(ploidies[,1]), varnames = "ploidy", value.name = "#individuals")
      df_n$ploidy <- as.factor(df_n$ploidy)
      p1 <- ggplot(df_n, aes(x = ploidy, y=`#individuals`, fill =ploidy)) +
        geom_bar(stat="identity") +
        theme_bw()

      if(is.vector(ploidies)){
        ploidies <- data.frame(Var1 = rownames(est.ploidy.chr_df()[[1]])[idx], Var2=names(ploidies),melt(ploidies, value.name = "ploidy"))
        ploidies.prop <- data.frame(Var1 = rownames(est.ploidy.chr_df()[[1]])[idx], Var2=names(ploidies.prop),melt(ploidies.prop, value.name = "prop"))
        ploidies.sd <- data.frame(Var1 = rownames(est.ploidy.chr_df()[[1]])[idx], Var2=names(ploidies.sd),melt(ploidies.sd, value.name = "sd"))
        ploidies.corr <- data.frame(Var1 = rownames(est.ploidy.chr_df()[[1]])[idx], Var2=names(ploidies.corr),melt(ploidies.corr, value.name = "corr"))
        ploidies.diff <- data.frame(Var1 = rownames(est.ploidy.chr_df()[[1]])[idx], Var2=names(ploidies.diff),melt(ploidies.diff, value.name = "diff"))
      } else {
        ploidies <- melt(ploidies, value.name = "ploidy")
        ploidies.prop <- melt(ploidies.prop, value.name = "prop")
        ploidies.sd <- melt(ploidies.sd, value.name = "sd")
        ploidies.corr <- melt(ploidies.corr, value.name = "corr")
        ploidies.diff <- melt(ploidies.diff, value.name = "diff")
      }

      df <- merge(ploidies, ploidies.prop, by = c(1,2))
      df <- merge(df, ploidies.sd, by = c(1,2))
      df <- merge(df, ploidies.corr, by = c(1,2))
      df <- merge(df, ploidies.diff, by = c(1,2))

      df <- melt(df, id.vars = c("ploidy", "Var1", "Var2"))

      df$variable <- gsub("sd", "standard deviation", df$variable)
      df$variable <- gsub("prop", "proportion in area", df$variable)
      df$variable <- gsub("corr", "correlation (Pearson)", df$variable)
      df$variable <- gsub("diff", "difference in proportion", df$variable)

      df$ploidy <- as.factor(df$ploidy)
      p2 <- ggplot(df, aes(x = ploidy, y=value, group=ploidy, fill=ploidy)) +
        geom_boxplot() +
        facet_grid(variable~., scales = "free") +
        theme_bw()

      p <- ggarrange(p1, p2, common.legend = T)

      # Weird individuals
      # set filters
      # All
      ploidy_all_weird <- unique(c(as.character(ploidies.corr$Var1)[which(ploidies.corr$corr < input$filter_corr)],
                                   as.character(ploidies.diff$Var1)[which(ploidies.diff$diff < input$filter_diff)]))

      if(length(ploidy_all_weird) == 0) ploidy_all_weird <- "Set filters"

      list(p, p_aneuploids, p_aneuploids2, text_ploidy, aneuploid_text, aneuploidy_df, ploidy_all_weird, aneuploidy_all_weird)
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
}

## To be copied in the UI
# mod_all_ui("all_1")

## To be copied in the server
# mod_all_server("all_1")
