library(shiny)
library(DT)
library(zip)
library(shinyalert)

# Define UI for application
ui <- fluidPage(
  # Add custom CSS to center the download button
  tags$head(
    tags$style(HTML("
      #download_link {
        display: flex;
        justify-content: center;
        align-items: center;
        height: 20vh; /* Full viewport height */
        margin: 0;
      }
      .btn-large {
        font-size: 1.5em; /* Larger font size */
        padding: 15px 30px; /* Larger padding */
      }
    "))
  ),
  
  titlePanel("Azenta Analyses App"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("raw_counts_input", "raw_counts.csv file",
                accept = ".csv"),
      fileInput("gene_lengths_input", "[sample_name].counts.txt file",
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
        selectInput("DGE_group1", "Choose first group:",
                    choices = NULL, # Will be populated by server
                    selected = NULL,
                    multiple = FALSE),  # Only one selection allowed
        selectInput("DGE_group2", "Choose group to compare to:",
                    choices = NULL, # Will be populated by server
                    selected = NULL,
                    multiple = FALSE)  # Only one selection allowed
      ),
      actionButton("run", "Run Analyses")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("All Results", 
                 textOutput("progress"),  # Display progress messages here
                 uiOutput("download_link")  # Display download link here
        ),
        tabPanel("Input Descriptions",
                 p("1. **raw_counts.csv**: This file can be found in the hit-counts folder of the Azenta results."),
                 p("2. **[sample_name].counts.txt file**: This file is found in the hit-counts folder of the Azenta results. Choose a file with any sample_name as long as the file name is in the format of {sample_name}.counts.txt, where {sample_name} is replaced with the name of the sample. It doesn't matter which file because this is only to get gene lengths."),
                 p("3. **Normalization method**: Choose how to normalize the raw counts. The options are RPKM (Reads Per Kilobase of transcript per Million mapped reads) and TPM (Transcripts Per Million)."),
                 p("4. **Expression Cutoff**: Only expression values higher than this cutoff are considered in the analysis. A cutoff of 1 is the default. Cutoff values lower than 1 are not recommended."),
                 p("5. **Select Analyses to Run**: Choose which analyses you want to perform, including Pathway analysis, Germline analysis, and Differential Gene Expression (DGE) analysis."),
                 p("6. **Choose first group for DGE analysis**: Select the first group for Differential Gene Expression analysis. The analysis will run comparing all of the replicates of this group to all of the replicates of the second group chosen"),
                 p("7. **Choose group to compare to**: Select a DIFFERENT group to compare against the first group for Differential Gene Expression analysis.")
        ),
        tabPanel("Pathway Table", 
                 DTOutput("pathway_table")  # Display Pathway analysis results here
        ),
        tabPanel("Germline Table", 
                 DTOutput("germline_table")  # Display Germline analysis results here
        )
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
  
  # Update the choices for DGE options based on the raw counts input
  observe({
    req(raw_counts())
    group_names <- colnames(raw_counts())
    replicate_names <- group_names[!group_names %in% c('ID', 'Length', 'Gene.name')]
    group_names <- unique(substr(replicate_names, 1, nchar(replicate_names) - 1))
    updateSelectInput(session, "DGE_group1", choices = group_names)
    updateSelectInput(session, "DGE_group2", choices = group_names)
  })
  
  # Reactive expression to store results of analyses
  results <- reactiveValues(data = NULL, pathway_data = NULL, germline_data = NULL)
  progress_message <- reactiveVal("")  # Store progress messages
  
  # Observe event when the Run button is pressed
  observeEvent(input$run, {
    req(raw_counts(), gene_length())
    
    # simple unit tests
    if (!(colnames(gene_length())[ncol(gene_length())] %in% colnames(raw_counts()))) {
      shinyalert::shinyalert("Error", "The raw_counts.csv and sample_name.counts.txt files are not from the same Azenta result", type = "error")
      return()
    }
    
    # Show progress message
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Running normalization...", value = 0.1)
    progress_message("Running normalization...")
    
    # Always run normalization
    source('./Scripts/Normalization.R')
    if (input$normalization_option == "RPKM") {
      
      results$data <- normalization_function(raw_counts = raw_counts(), 
                                             gene_length = gene_length(), 
                                             normalization_type = 'RPKM',
                                             exp_cutoff = input$exp_cutoff)
      
      
      
    } else if (input$normalization_option == "TPM") {
      results$data <- normalization_function(raw_counts = raw_counts(), 
                                             gene_length = gene_length(), 
                                             normalization_type = 'TPM',
                                             exp_cutoff = input$exp_cutoff)
    }
    
    # check that each group has the same number of replicates
    replicate_groups <- read.table("./Data/replicate_groups.txt", quote="\"", comment.char="")
    t <- table(replicate_groups)
    
    if (min(t) != max(t)) {
      shinyalert::shinyalert("Error", "The number of replicates varies across groups.", type = "error")
      return()
    }
    
    progress$set(message = "Normalization complete", value = 0.3)
    progress_message("Normalization complete")
    
    # Run other selected analyses
    if ("pathway_analysis" %in% input$analyses) {
      source('./Scripts/Pathway_Analysis.R')
      progress$set(message = "Running pathway analysis...", value = 0.5)
      progress_message("Running pathway analysis...")
      results$pathway_data <- pathway_function(normalization_type = input$normalization_option)
      
      output$pathway_table <- renderDT({
        datatable(results$pathway_data, options = list(pageLength = 25, searchHighlight = TRUE))
      })
    }
    if ("germline_analysis" %in% input$analyses) {
      source('./Scripts/Germline_Analysis.R')
      progress$set(message = "Running germline analysis...", value = 0.7)
      progress_message("Running germline analysis...")
      results$germline_data <- germline_function(normalization_type = input$normalization_option)
      
      output$germline_table <- renderDT({
        datatable(results$germline_data, options = list(pageLength = 25, searchHighlight = TRUE))
      })
    }
    if ("DGE_analysis" %in% input$analyses) {
      req(input$DGE_group1, input$DGE_group2)
      
      source('./Scripts/DGE_Analysis.R')
      progress$set(message = "Running DGE analysis...", value = 0.9)
      progress_message("Running DGE analysis...")
      results$data$DGE <- DGE_function(selected_group_1 = input$DGE_group1,
                                       selected_group_2 = input$DGE_group2)
      
    }
    
    # Finalize progress
    progress$set(message = "All analyses complete", value = 1)
    progress_message("All analyses are complete.")
    
    # Create zip file of the ./Data folder
    output$download_link <- renderUI({
      tagList(
        # Create a large download button
        downloadButton("download_zip", "Download Results", class = "btn-large")
      )
    })
    
    output$download_zip <- downloadHandler(
      filename = function() {
        paste("Analysis_Results_", Sys.Date(), ".zip", sep = "")
      },
      content = function(file) {
        zip(zipfile = file, files = dir("./Data", full.names = TRUE))
      }
    )
  })
  
  # Show progress messages in mainPanel
  output$progress <- renderText({
    progress_message()
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
