


library(shiny)

# Define UI for application
ui <- fluidPage(
  titlePanel("Azenta Analyses App"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("raw_counts_input", "Provide the raw_counts.csv file found in the hit-counts folder",
                accept = ".csv"),
      fileInput("gene_lengths_input", "Provide the gene length file (expression file from any sample)",
                accept = ".txt"),
      radioButtons("normalization_option", "Normalization method:",
                   choices = list("RPKM" = "RPKM", "TPM" = "TPM"),
                   selected = "RPKM"), # Default option
      numericInput("exp_cutoff", "Expression Cutoff:", value = 1, min = 0, step = 0.1),
      checkboxGroupInput("analyses", "Select Analyses to Run:",
                         choices = list("Pathway" = "pathway_analysis",
                                        "Germline" = "germline_analysis",
                                        "Differential Gene Expression" = "DGE_analysis")),
      conditionalPanel(
        condition = "input.analyses.includes('DGE_analysis')",
        selectInput("DGE_option1", "Choose first set of values for DGE analysis:",
                    choices = NULL, # Will be populated by server
                    selected = NULL,
                    multiple = TRUE),
        selectInput("DGE_option2", "Choose second set of values for DGE analysis:",
                    choices = NULL, # Will be populated by server
                    selected = NULL,
                    multiple = TRUE)
      ),
      actionButton("run", "Run Analyses")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Data", tableOutput("data")),
        tabPanel("Progress", textOutput("progress")) # Changed to display progress messages
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  options(shiny.maxRequestSize = 50 * 1024^2)
  
  # Reactive expression to read the uploaded raw counts file
  raw_counts <- reactive({
    req(input$raw_counts_input)
    inFile <- input$raw_counts_input
    df <- read.csv(inFile$datapath)
    return(df)
  })
  
  # Reactive expression to read the uploaded gene length file
  gene_length <- reactive({
    req(input$gene_lengths_input)
    inFile <- input$gene_lengths_input
    df <- read.delim(inFile$datapath)
    return(df)
  })
  
  
  
  # Display the uploaded raw counts and gene length data
  output$data <- renderTable({
    req(raw_counts(), gene_length())
    list(RawCounts = raw_counts(), GeneLength = gene_length())
  })
  
  # Reactive expression to store results of analyses
  results <- reactiveValues(data = NULL)
  progress_message <- reactiveVal("") # Store progress messages
  
  # Observe event when the Run button is pressed
  observeEvent(input$run, {
    req(raw_counts(), gene_length())
    
    # Show progress message
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Running normalization...", value = 0.1)
    
    # Always run normalization
    source('./Scripts/Normalization.R')
    if (input$normalization_option == "RPKM") {
      results$data <- cbind(results$data, normalization_function(raw_counts = raw_counts(), 
                                                                 gene_length = gene_length(), 
                                                                 normalization_type = 'RPKM',
                                                                 exp_cutoff = input$exp_cutoff))
    } else if (input$normalization_option == "TPM") {
      results$data <- cbind(results$data, normalization_function(raw_counts = raw_counts(), 
                                                                 gene_length = gene_length(), 
                                                                 normalization_type = 'TPM',
                                                                 exp_cutoff = input$exp_cutoff))
    }
    progress$set(message = "Normalization complete", value = 0.3)
    progress_message("Normalization complete")
    
    # Run other selected analyses
    if ("pathway_analysis" %in% input$analyses) {
      source('./Scripts/Pathway_Analysis.R')
      results$data <- pathway_function(normalization_type = input$normalization_option)
      progress$set(message = "Running pathway analysis...", value = 0.5)
      progress_message("Running pathway analysis...")
    }
    if ("germline_analysis" %in% input$analyses) {
      source('./Scripts/Germline_Analysis.R')
      results$data <- germline_function(normalization_type = input$normalization_option)
      progress$set(message = "Running germline analysis...", value = 0.7)
      progress_message("Running germline analysis...")
    }
    if ("DGE_analysis" %in% input$analyses) {
      req(input$DGE_option1, input$DGE_option2)
      
      source('./Scripts/DGE_Analysis.R')
      
      results$data <- DGE_function(selected_group_1 = input$DGE_option1,
                                   selected_group_2 = input$DGE_option2)
      progress$set(message = "Running DGE analysis...", value = 0.9)
      progress_message("Running DGE analysis...")
    }
    
    # Finalize progress
    progress$set(message = "All analyses complete", value = 1)
    progress_message("All analyses complete")
  })
  
  # Display progress messages
  output$progress <- renderText({
    progress_message()
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)
