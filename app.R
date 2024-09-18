library(shiny)
library(DT)
library(zip)
library(shinyalert)
library(dplyr)
library(reshape2)
library(tidyr)
library(fs)
library(ggplotify)
library(DESeq2)
library(ggplot2)
library(tibble)
library(ragg)
library(ComplexHeatmap)

#aws_prefix <- '/mnt/efs/fs1/destination_folder/Azenta_Analyses/'
aws_prefix <- 'C:/Users/17735/Downloads/Azenta_Analyses/'

# Update UI for application
ui <- fluidPage(
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
  
  sidebarPanel(
    selectInput("project", "Select Project Directory:",
                choices = NULL, # Will be populated by server
                selected = NULL),
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
               textOutput("progress"),
               uiOutput("download_link")
      ),
      tabPanel("Input Descriptions",
               p("**Select Project Directory**: Choose an Azenta project to run analyses for. Only projects inside of the Azenta_Projects directory in AWS are able to be chosen."),
               p("**Normalization method**: Choose how to normalize the raw counts. The options are RPKM (Reads Per Kilobase of transcript per Million mapped reads) and TPM (Transcripts Per Million)."),
               p("**Expression Cutoff**: Only expression values higher than this cutoff are considered in the analysis. A cutoff of 1 is the default. Cutoff values lower than 1 are not recommended."),
               p("**Select Analyses to Run**: Choose which analyses you want to perform, including Pathway analysis, Germline analysis, and Differential Gene Expression (DGE) analysis."),
               p("**Choose first group for DGE analysis**: Select the first group for Differential Gene Expression analysis. The analysis will run comparing all of the replicates of this group to all of the replicates of the second group chosen."),
               p("**Choose group to compare to**: Select a DIFFERENT group to compare against the first group for Differential Gene Expression analysis.")
      ),
      tabPanel("Pathway Table", 
               DTOutput("pathway_table")
      ),
      tabPanel("Germline Table", 
               DTOutput("germline_table")
      )
    )
  )
)
server <- function(input, output, session) {
  
  options(shiny.maxRequestSize = 50 * 1024^2)
  
  # Reactive expression to list directories in aws_prefix
  available_dirs <- function() {
    dirs <- dir(path = paste0(aws_prefix, 'Azenta_Projects'), full.names = TRUE, recursive = FALSE, include.dirs = TRUE)
    dir_names <- basename(dirs)
    return(dir_names)
  }
  
  # Update the directory selection inputs with available directories
  observe({
    updateSelectInput(session, "project", choices = available_dirs(), selected = input$project)
  })
  
  # Reactive expression to read the raw counts file based on selected project
  project_path <- reactive({
    req(input$project)
    # Construct the full path to the raw_counts.csv file
    path <- file.path(gsub('/$', '', aws_prefix), 'Azenta_Projects', input$project, '/')
    if (!file.exists(path)) {
      stop("File does not exist: ", path)
    }
    return(path)
  })
  
  # Update the choices for DGE options based on the raw counts input
  observe({
    req(project_path())
    raw_counts <- read.csv(paste0(project_path(), 'hit-counts/raw_counts.csv'))  # Use project_path() directly
    
    # Assuming the columns you want to exclude are 'ID', 'Length', and 'Gene.name'
    group_names <- colnames(raw_counts)
    replicate_names <- group_names[!group_names %in% c('ID', 'Length', 'Gene.name')]
    
    # make it work for the triple replicates labeled 1,2,nothing instead of 1,2,3
    if (all(grepl("(1|2|3)$", replicate_names))) {
      # all values end in 1, 2, or 3  -> remove last character 
      group_names <- lapply(replicate_names, function(x) substr(x, 1, nchar(x) - 1))
    } else {
      # Use sub to remove only the last dot and numbers that follow (only if only numbers follow)
      group_names <- sub("\\.[0-9]+$", "", replicate_names)
    }  
    
    # get list of group names 
    group_names <- as.character(unique(group_names))
    
    # Update the select input choices for DGE groups
    updateSelectInput(session, "DGE_group1", choices = group_names)
    updateSelectInput(session, "DGE_group2", choices = group_names)
    
  })
  
  # Reactive expression to store results of analyses
  results <- reactiveValues(data = NULL, pathway_data = NULL, germline_data = NULL)
  progress_message <- reactiveVal("")  # Store progress messages
  
  # Observe event when the Run button is pressed
  observeEvent(input$run, {
    req(project_path())
    
    
    if (("DGE_analysis" %in% input$analyses) & (input$DGE_group1 == input$DGE_group2)) {
      shinyalert::shinyalert("Error", "The same groups were chosen for DGE Analysis.", type = "error")
      return()
    }
    
    # remove old results 
    unlink(paste0(aws_prefix, "Data/*"))
    
    # Show progress message
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Running normalization...", value = 0.1)
    progress_message("Running normalization...")
    
    # Normalization logic remains the same
    source(file.path(aws_prefix, 'Scripts', 'Normalization.R'))
    if (input$normalization_option == "RPKM") {
      results$data <- normalization_function(project_path = project_path(), 
                                             normalization_type = 'RPKM',
                                             exp_cutoff = input$exp_cutoff)
    } else if (input$normalization_option == "TPM") {
      results$data <- normalization_function(project_path = project_path(), 
                                             normalization_type = 'TPM',
                                             exp_cutoff = input$exp_cutoff)
    }
    
    # Check that each group has the same number of replicates
    replicate_groups <- read.table(file = paste0(aws_prefix, "Data/replicate_groups.txt"), quote="\"", comment.char="")
    t <- table(replicate_groups)
    
    if (min(t) != max(t)) {
      shinyalert::shinyalert("Error", "The number of replicates varies across groups.", type = "error")
      return()
    }
    
    progress$set(message = "Normalization complete", value = 0.3)
    progress_message("Normalization complete")
    
    # Run other selected analyses
    if ("pathway_analysis" %in% input$analyses) {
      source(file.path(aws_prefix, 'Scripts', 'Pathway_Analysis.R'))
      progress$set(message = "Running pathway analysis...", value = 0.5)
      progress_message("Running pathway analysis...")
      results$pathway_data <- pathway_function(normalization_type = input$normalization_option)
      
      output$pathway_table <- renderDT({
        datatable(results$pathway_data, options = list(pageLength = 25, searchHighlight = TRUE))
      })
    }
    if ("germline_analysis" %in% input$analyses) {
      source(file.path(aws_prefix, 'Scripts', 'Germline_Analysis.R'))
      progress$set(message = "Running germline analysis...", value = 0.7)
      progress_message("Running germline analysis...")
      results$germline_data <- germline_function(normalization_type = input$normalization_option)
      
      output$germline_table <- renderDT({
        datatable(results$germline_data, options = list(pageLength = 25, searchHighlight = TRUE))
      })
    }
    if ("DGE_analysis" %in% input$analyses) {
      req(input$DGE_group1, input$DGE_group2)
      
      source(file.path(aws_prefix, 'Scripts', 'DGE_Analysis.R'))
      progress$set(message = "Running DGE analysis...", value = 0.9)
      progress_message("Running DGE analysis...")
      results$data$DGE <- DGE_function(selected_group_1 = input$DGE_group1,
                                       selected_group_2 = input$DGE_group2)
    }
    
    # Finalize progress
    progress$set(message = "All analyses complete", value = 1)
    progress_message("All analyses are complete.")
    
    # Create zip file of the Data folder
    output$download_link <- renderUI({
      tagList(
        # Create a large download button
        downloadButton("download_zip", "Download Results", class = "btn-large")
      )
    })
    
    output$download_zip <- downloadHandler(
      filename = function() {
        paste("Analysis_Results_", input$project, '_', Sys.Date(), ".zip", sep = "")
      },
      content = function(file) {
        zip::zipr(file, 
                 files = list.files(paste0(aws_prefix, "Data/"), full.names = T), 
                 recurse = F, 
                 include_directories = F)
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
