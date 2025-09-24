#' Home UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @importFrom bs4Dash renderValueBox valueBox
#'
#'
mod_Home_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      fluidRow(
        column(width = 4,
               box(
                 title = "Qploidy + BIGapp", status = "info", solidHeader = FALSE, width = 12, collapsible = FALSE,
                 HTML(
                   "<p><b>About Qploidy.</b> Qploidy package provides a standardization method for allele counts or intensities that allows ploidy and aneuploidy estimation. It was first developed by Cristiane (Cris) Taniguti while working with Oscar Riera-Lizarazu’s group at Texas AM University. After Cris moved to Breeding Insight (BI), Qploidy’s maintenance and new features have continued under BI—this interface is a direct result of that ongoing effort. Qploidy is now integrated into BIGapp.</p>

     <p><b>About BIGapp</b> BIGapp is a user-friendly R Shiny application that streamlines low- to mid-density genotyping workflows for diploid and polyploid species. It provides a web-based interface so users can analyze genomic data without command-line tools. New analyses will be added over time, with an initial focus on features that support breeding decisions.</p>

     <p><b>Supported Analyses</b></p>
     <p>Initial supported analyses draw from mature genomics/bioinformatics pipelines developed within Breeding Insight:</p>
     <ul>
       <li>Genotype Processing (incl. Qploidy standardization)</li>
       <li>Summary Metrics</li>
       <li>Population Structure</li>
       <li>GWAS</li>
       <li>Genomic Selection</li>
       <li>Ploidy Estimation</li>
     </ul>"
                 ),
                 style = "overflow-y: auto; height: 500px"
               )
        ),
        column(width = 4,
               box(
                 title = "About Breeding Insight", status = "success", solidHeader = FALSE, width = 12, collapsible = FALSE,
                 HTML(
                   "We provide scientific consultation and data management software to the specialty crop and animal breeding communities.
            <ul>
              <li>Genomics</li>
              <li>Phenomics</li>
              <li>Data Management</li>
              <li>Software Tools</li>
              <li>Analysis</li>
            </ul>
            Breeding Insight is funded by the U.S. Department of Agriculture (USDA) Agricultural Research Service (ARS) through Cornell University.
            <div style='text-align: center; margin-top: 20px;'>
              <img src='www/BreedingInsight.png' alt='Breeding Insight' style='width: 85px; height: 85px;'>
            </div>"
                 ),
                 style = "overflow-y: auto; height: 500px"
               )
        ),
        column(width = 4,
               a(
                 href = "https://www.breedinginsight.org",  # Replace with your desired URL
                 target = "_blank",  # Optional: opens the link in a new tab
                 valueBox(
                   value = NULL,
                   subtitle = "Learn More About Breeding Insight",
                   icon = icon("link"),
                   color = "purple",
                   gradient = TRUE,
                   width = 11
                 ),
                 style = "text-decoration: none; color: inherit;"  # Optional: removes underline and retains original color
               ),
               a(
                 href = "https://breedinginsight.org/contact-us/",  # Replace with your desired URL
                 target = "_blank",  # Optional: opens the link in a new tab
                 valueBox(
                   value = NULL,
                   subtitle = "Contact Us",
                   icon = icon("envelope"),
                   color = "danger",
                   gradient = TRUE,
                   width = 11
                 ),
                 style = "text-decoration: none; color: inherit;"  # Optional: removes underline and retains original color
               ),
               a(
                 href = "file:///Users/cht47/Documents/github/Qploidy/doc/Qploidy.html",  # Replace with your desired URL
                 target = "_blank",  # Optional: opens the link in a new tab
                 valueBox(
                   value = NULL,
                   subtitle = "Qploidy Tutorial",
                   icon = icon("compass"),
                   color = "info",
                   gradient = TRUE,
                   width = 11
                 ),
                 style = "text-decoration: none; color: inherit;"  # Optional: removes underline and retains original color
               ),
               a(
                 href = "https://scribehow.com/page/BIGapp_Tutorials__FdLsY9ZxQsi6kgT9p-U2Zg",  # Replace with your desired URL
                 target = "_blank",  # Optional: opens the link in a new tab
                 valueBox(
                   value = NULL,
                   subtitle = "BIGapp Tutorials",
                   icon = icon("compass"),
                   color = "warning",
                   gradient = TRUE,
                   width = 11
                 ),
                 style = "text-decoration: none; color: inherit;"  # Optional: removes underline and retains original color
               )
        )
      )
    )
  )
}

#' Home Server Functions
#'
#'
#' @noRd
mod_Home_server <- function(input, output, session, parent_session){

  ns <- session$ns

}

## To be copied in the UI
# mod_Home_ui("Home_1")

## To be copied in the server
# mod_Home_server("Home_1")
