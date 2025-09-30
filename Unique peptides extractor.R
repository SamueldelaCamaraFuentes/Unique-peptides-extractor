library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(DT)
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(shinyjs)
library(writexl)

options(shiny.maxRequestSize = 30*1024^2)  # 30 MB max upload

# ---- UI ----
ui <- dashboardPage(
  dashboardHeader(title = "Unique Peptides Extractor"),
  
  dashboardSidebar(
    useShinyjs(),
    sidebarMenu(
      menuItem("Data Input", tabName = "input_tab", icon = icon("upload"))
    ),
    
    div(style = "padding: 10px;",
        prettyRadioButtons(
          inputId = "data_type",
          label = "Select Data Type:",
          choices = c("Proteome Discoverer" = "pd",
                      "MSFragger" = "msfragger"),
          icon = icon("check"),
          status = "info",
          animation = "smooth"
        )
    ),
    
    conditionalPanel(
      condition = "input.data_type == 'pd'",
      fileInput("pd_file1", "Upload Proteome Discoverer peptide file (.xlsx)"),
      fileInput("pd_file2", "Upload Proteome Discoverer protein file (.xlsx)")
    ),
    
    conditionalPanel(
      condition = "input.data_type == 'msfragger'",
      fileInput("msfragger_file", "Upload MSFragger protein file (.tsv)"),
      textInput("msfragger_dir", "Enter MSFragger Directory Path")
    ),
    
    br(),
    downloadButton("downloadData", "Download Dataset", class = "btn-primary")
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper {background-color: #f7f7f7;}
        .box {border-radius: 10px;}
        .box-header {font-weight: bold; font-size: 16px;}
      "))
    ),
    
    fluidRow(
      box(title = "Dataset Preview", width = 12, status = "primary",
          solidHeader = TRUE,
          DTOutput("table")
      )
    )
  )
)

# ---- SERVER ----
server <- function(input, output, session) {
  
  dataset <- reactive({
    req(input$data_type)
    
    # --- PROTEOME DISCOVERER ---
    if (input$data_type == "pd") {
      req(input$pd_file1, input$pd_file2)
      
      raw <- readxl::read_xlsx(input$pd_file1$datapath)
      proteins <- readxl::read_xlsx(input$pd_file2$datapath)
      
      abundance_cols <- grep("Abundance:", colnames(raw), value = TRUE)
      
      df_peptides <- raw %>%
        filter(!grepl("Oxidation", Modifications)) %>%
        distinct(Sequence, .keep_all = TRUE) %>%
        select(Sequence, `Master Protein Accessions`, all_of(abundance_cols)) %>%
        pivot_longer(cols = all_of(abundance_cols),
                     names_to = "Sample",
                     values_to = "Abundance") %>%
        filter(!is.na(Abundance)) %>%
        group_by(`Master Protein Accessions`, Sample) %>%
        summarise(Peptide_Count = n(), .groups = "drop") %>%
        pivot_wider(names_from = Sample, values_from = Peptide_Count, values_fill = 0) %>%
        rename_with(~ str_replace(., "^Abundance:", "Unique.peptides"), -`Master Protein Accessions`) %>%
        mutate(across(starts_with("Unique.peptides"), as.numeric))
      
      filtered_df <- df_peptides %>%
        inner_join(proteins, by = c("Master Protein Accessions" = "Accession")) %>%
        rename(Accession = `Master Protein Accessions`)
      
      return(filtered_df)
    }
    
    # --- MSFRAGGER ---
    if (input$data_type == "msfragger") {
      req(input$msfragger_file, input$msfragger_dir)
      
      folders <- list.dirs(path = input$msfragger_dir, recursive = FALSE)
      protein_data_list <- list()
      
      for (folder in folders) {
        file_path <- file.path(folder, "protein.tsv")
        
        if (file.exists(file_path)) {
          df <- tryCatch(
            read_tsv(file_path, show_col_types = FALSE),
            error = function(e) NULL
          )
          
          if (!is.null(df) && all(c("Protein", "Unique Peptides") %in% names(df))) {
            sample_name <- paste0(basename(folder), ".Unique.peptides")
            df <- df %>%
              select(Protein, `Unique Peptides`) %>%
              rename(!!sample_name := `Unique Peptides`)
            protein_data_list[[sample_name]] <- df
          }
        }
      }
      
      combined_df <- Reduce(function(x, y) full_join(x, y, by = "Protein"), protein_data_list)
      combined_df[is.na(combined_df)] <- 0
      combined_df <- combined_df %>%
        mutate(across(ends_with("Unique.peptides"), ~ as.numeric(.)))
      
      raw <- read.delim(input$msfragger_file$datapath,
                        sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
      
      merged <- merge(raw, combined_df, by = "Protein")
      return(merged)
    }
  })
  
  # ---- Dataset Preview ----
  output$table <- renderDT({
    req(dataset())
    datatable(
      head(dataset(), 20),
      options = list(
        scrollX = TRUE,
        pageLength = 10,
        lengthMenu = c(10, 20, 50),
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      rownames = FALSE
    )
  })
  
  # ---- Download Handler ----
  output$downloadData <- downloadHandler(
    filename = function() {
      if (input$data_type == "pd") {
        paste0("proteome_discoverer_dataset_", Sys.Date(), ".xlsx")
      } else if (input$data_type == "msfragger") {
        paste0("msfragger_dataset_", Sys.Date(), ".tsv")
      } else {
        paste0("dataset_", Sys.Date(), ".csv")
      }
    },
    content = function(file) {
      req(dataset())
      if (input$data_type == "pd") {
        writexl::write_xlsx(dataset(), file)
      } else if (input$data_type == "msfragger") {
        write.table(dataset(), file, sep = "\t", row.names = FALSE, quote = FALSE)
      } else {
        write.csv(dataset(), file, row.names = FALSE)
      }
    }
  )
}

# ---- RUN APP ----
shinyApp(ui, server)

