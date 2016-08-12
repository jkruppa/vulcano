library(shiny)

shinyUI(fluidPage(
  titlePanel("Visualisation of different normalization methods"),
  sidebarLayout(
      ## sidebar panel
      sidebarPanel(
          ## Controls
          h3("Controls"),
          sliderInput("fc.thresh", 
                      "Absolute log fold change", 
                      min = 0,
                      max = 2,
                      value = 0.5,
                      step = 0.25),
          sliderInput("alpha.error", 
                      "Type I error", 
                      min = 0.001,
                      max = 0.1,
                      value = 0.05,
                      step = 0.01),
          helpText("Data was generated at random. The false discovery rate (FDR) on the y axis is calculated by -log10 of the Type I error."),
          helpText("To avoid overplotting some gene symbols are not displayed in the vulcano plot but can be found in the corresponding table."),
          hr(),
          ## Dependencies
          h3("Dependencies"),
          p("The following R packages are needed for the shiny app."),
          code('install.packages("shiny")'),
          br(),
          code('install.packages("ggplot2")'), 
          hr(),
          ## References
          h3("References"),
          p("Cui X and Churchill GA (2003). 'Statistical tests for differential expression in cDNA microarray experiments'. Genome Biology, 4(4), 210.")
      ),
      ## main panel
      mainPanel(
          h4("Vulcano plot"),
          plotOutput("plots", width = "100%", height = "600px"),
          h4("Significant hits ordered by log fold change"),
          tableOutput("view")
      ),
  )
))
