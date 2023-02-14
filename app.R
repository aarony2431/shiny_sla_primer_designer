# Import R packages needed for the UI
library(shiny)
library(dplyr)
library(tidyverse)
library(gsubfn)

# ------------------ App virtualenv setup (Do not edit) ------------------- #

virtualenv_dir = Sys.getenv('VIRTUALENV_NAME')
python_path = Sys.getenv('PYTHON_PATH')

# Define any Python packages needed for the app here:
PYTHON_DEPENDENCIES = c('pip', 'biopython', 'Cython', 'primer3-py')

# Create virtual env and install dependencies
reticulate::virtualenv_create(envname = virtualenv_dir, python = python_path)
reticulate::virtualenv_install(virtualenv_dir, packages = PYTHON_DEPENDENCIES, ignore_installed=TRUE)
reticulate::use_virtualenv(virtualenv_dir, required = T)

# --------------------------------- End ----------------------------------- #

# Source Python code
reticulate::source_python('Shiny_SLA_function.py')

# Begin UI for the R
ui <- fluidPage(
  
  titlePanel('SLA Primer Designer Tool'),
  
  sidebarLayout(
    
    # Sidebar with a text input
    sidebarPanel(
      helpText(
        tags$p(
          'The following tool was implemented and edited by Aaron Yu (2022) using code originally written by Sam Beppler (2020). 
          This tool picks the primers which yield the highest number of probes and the overall longest in 50 iterations of optimization.
          The Stem-Loop Assay is a qPCR method designed to detect the Guide Strand from the 3\' end.'
        ), 
        tags$p(
          'Please enter the full Guide Strand (GS; Antisense Strand) sequence including the ', tags$em('\'UN\''),  'overhangs.
          For example, an siRNA with a target antisense of ', tags$em('\'AATTGTGCATTCCCGTGCA\''),  'will have a final GS 
          sequence of ', tags$em('\'UAUUGUGCAUUCCCGUGCAUU\''), '.'
        ), 
        tags$p(
          tags$strong('Please note'), 'that it is possible to get different primers using the same sequence. 
          This is due to to a pseudo-random generator for GC content and melting temperature.'
        )
      ), 
      hr(), 
      textAreaInput('GS', 
                    'Enter full GS sequences in the area below:', 
                    placeholder = 'Enter each GS on a separate line...', 
                    cols = 30, 
                    rows = 20
                    ),
      actionButton('getPrimers', 
                   'Calculate Primers!', 
                   class = 'btn btn-primary', 
                   width = '100%'
                   ), 
      hr(), 
      verbatimTextOutput('inputText')
    ),
    
    # Show primer output
    mainPanel(
      downloadButton('downloadData', 'Download Output'),
      verbatimTextOutput('primerTable')
    )
  )
)

# Begin app server
server <- function(input, output) {
  
  # ------------------ App server logic (Edit anything below) --------------- #
  
  # Initialize
  outputs <- ''
  
  getDownloadOutput <- function(sequences.numbered, primers.output) {
    oligos <- 8
    Guide_Strand <- rep(sequences.numbered, each = oligos)
    primers.df <- lapply(primers.output, function(s) {
      lapply(s, function(p) {
        str_split(p, pattern = '\t') %>% unlist(.)
      }) %>% as.data.frame(.)
    }) %>% Reduce(cbind, .) %>% as.data.frame(.) %>% t(.) %>% as.data.frame(.)
    rownames(primers.df) <- 1:nrow(primers.df)
    colnames(primers.df) <- c('Primers\\Probes', 'Sequence', 'TM')
    out <- cbind(data.frame(Guide_Strand), primers.df)
    return(out)
  }
  
  # Add reactive event to action button
  primers <- eventReactive(input$getPrimers, {
    in.in <- input$GS
    sequences <- str_split(in.in, pattern = '(\\n|\\s)') %>% unlist(.)
    sequences.trimmed <- sequences[sequences != '']
    
    output$inputText <- renderText({
      paste('Your Input Text:', paste(sequences.trimmed, collapse = '\n'), sep = '\n')
    })
    
    primers.output <- lapply(sequences.trimmed, cross_check_sla)
    sequences.numbered <- paste(1:length(sequences.trimmed), sequences.trimmed, sep = '. ')
    outputs <<- getDownloadOutput(sequences.numbered, primers.output)
    return(list(sequences.numbered, primers.output))
  })
  
  # Assign output text
  output$primerTable <- renderText({
    list[sequences.numbered, primers.output] <- primers()
    table <- paste(sequences.numbered, lapply(primers.output, function (x) {
      spacer <- rep('\t', times = length(x))
      paste(spacer, x, sep = '', collapse = '\n')
    }), sep = '\n', collapse = '\n\n')
    return(table)
  })
  
  # Download Data
  output$downloadData <- downloadHandler(
    filename = 'sla_primers.csv',
    content = function(file) {
      write.csv(outputs, file)
    },
    contentType = 'text/csv'
  )
  
}

shinyApp(ui, server)