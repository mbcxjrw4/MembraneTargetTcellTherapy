# app.R

library(shiny)
library(ggplot2)
library(DT)

# Source your modularized functions
source("scripts/01_select_membrane_targets.R")
source("scripts/02_extract_expression_data.R")
source("scripts/03_percentage_above_cutoff.R")
source("scripts/04_apply_safety_filters.R")
source("scripts/05_visualize_overview.R")

# UI ----
ui <- fluidPage(
  titlePanel("Membrane Target Discovery for Cell Therapy"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("compartments_score",
                  "COMPARTMENTS Knowledge Score Threshold:",
                  min = 1, max = 5, value = 3, step = 1),
      
      sliderInput("tpm_cutoff", 
                  "TPM Cutoff (Tumor Expression Threshold):", 
                  min = 1, max = 50, value = 10, step = 1),
      
      sliderInput("tumor_prevalence", 
                  "Tumor Prevalence Threshold (%):", 
                  min = 0, max = 100, value = 20, step = 1),
      
      actionButton("run_button", "Run Analysis", class = "btn-primary")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Overview Plot", 
                 plotOutput("overviewPlot", height = "500px")),
        tabPanel("Target Table", 
                 DTOutput("targetTable")),
        tabPanel("Download",
                 downloadButton("downloadData", "Download Filtered Targets"))
      )
    )
  )
)

# Server ----
server <- function(input, output, session) {
  
  results <- reactiveVal(NULL)
  
  observeEvent(input$run_button, {
    # Step 1: Select membrane targets based on knowledge score
    membrane_targets <- select_membrane_targets(score_threshold = input$compartments_score)
    
    # Step 2: Extract RNA expression data
    expression_data <- extract_expression_data(membrane_targets)
    
    # Step 3: Apply safety filters
    filtered <- apply_safety_filters(
      expression_data = expression_data,
      tpm_cutoff = input$tpm_cutoff,
      prevalence_threshold = input$tumor_prevalence
    )
    
    # Save results
    results(filtered)
  })
  
  # Overview plot
  output$overviewPlot <- renderPlot({
    req(results())
    plot_overview(results())
  })
  
  # Target result table
  output$targetTable <- renderDT({
    req(results())
    datatable(results(), options = list(pageLength = 10), rownames = FALSE)
  })
  
  # Download
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("filtered_targets_score", input$compartments_score,
             "_tpm", input$tpm_cutoff,
             "_prev", input$tumor_prevalence, ".csv")
    },
    content = function(file) {
      write.csv(results(), file, row.names = FALSE)
    }
  )
}

# Run app
shinyApp(ui, server)
