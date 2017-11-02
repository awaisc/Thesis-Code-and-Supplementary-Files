
library(shiny)
library(shinycssloaders)
library(shinydashboard)
library(gridExtra)


dashboardPage(
  dashboardHeader(title = "ARX Genome Browser"),
  dashboardSidebar(
    numericInput("fromM", "Starting Base",value = 25016813),
    numericInput("toM", "Finishing Base", value = 25038065),
    numericInput("counts", "Contact Probability", value = 1),
    selectInput("chrM", label = h3("Select box"), 
                choices = list("Chromosome 1" = "chr1",
                               "Chromosome 2" = "chr2", 
                               "Chromosome 3" = "chr3",
                               "Chromosome 4" = "chr4",
                               "Chromosome 5" = "chr5",
                               "Chromosome 6" = "chr6",
                               "Chromosome 7" = "chr7",
                               "Chromosome 8" = "chr8",
                               "Chromosome 9" = "chr9",
                               "Chromosome 10" = "chr10",
                               "Chromosome 11" = "chr11",
                               "Chromosome 12" = "chr12",
                               "Chromosome 13" = "chr13",
                               "Chromosome 14" = "chr14",
                               "Chromosome 15" = "chr15",
                               "Chromosome 16" = "chr16",
                               "Chromosome 17" = "chr17",
                               "Chromosome 18" = "chr18",
                               "Chromosome 19" = "chr19",
                               "Chromosome 20" = "chr20",
                               "Chromosome 21" = "chr21",
                               "Chromosome 22" = "chr22",
                               "Chromosome X" = "chrX",
                               "Chromosome Y" = "chrY"), selected = "chrX"),
    hr(),
    fluidRow(column(3, verbatimTextOutput("value"))
    ),
    
    # Copy the line below to make a checkbox
    checkboxInput("contactProbabilities", label = "Raw Interactions", value = TRUE)
    
  ),
  dashboardBody(tabsetPanel(
    tabPanel("Human (hg19) Genome Browser",      
             withSpinner(plotOutput("HumangvizPlot"),
                         type = getOption("spinner.type", default = 3),
                         color = getOption("spinner.color", default = "#0275D8"),
                         color.background = getOption("spinner.color.background", default = "#FFFFFF"))),
    tabPanel("Mouse (mm9) Genome Browser", 
             withSpinner(plotOutput("MousegvizPlot"),
                         type = getOption("spinner.type", default = 3),
                         color = getOption("spinner.color", default = "#0275D8"),
                         color.background = getOption("spinner.color.background", default = "#FFFFFF"))),
    tabPanel("Legends", 
             fluidRow(
               column(6,plotOutput(outputId="LegendsPlot", width="300px",height="700px")),  
               column(6,plotOutput(outputId="LegendsPlotMouse", width="300px",height="700px"))
               
             # withSpinner(plotOutput(c("LegendsPlot")),
             #             type = getOption("spinner.type", default = 3),
             #             color = getOption("spinner.color", default = "#0275D8"),
             #             color.background = getOption("spinner.color.background", default = "#FFFFFF")),
             # withSpinner(plotOutput(c("LegendsPlotMouse")),
             #             type = getOption("spinner.type", default = 3),
             #             color = getOption("spinner.color", default = "#0275D8"),
             #             color.background = getOption("spinner.color.background", default = "#FFFFFF"))
             ))
    )
    )
    )