library(dplyr)
library(shiny)
library(shinythemes)
library(stringr)
library(plotly)
library(ggplot2)
library(DT)
library(shinyFeedback)

# Helper functions ----
clean_dna <- function(dna, clean) {
  if(clean) {
    toupper(str_replace_all(dna, "[^ATGC]", ""))
  } else {
    toupper(dna)
  }
}

calculate_stats <- function(dna, rna, cleaned) {
  data.frame(
    Metric = c("RNA Length", "A", "U", "G", "C", "GC Content (%)",
               "Original Gaps", "Original Ns", "Invalid in DNA", "Invalid in RNA"),
    Value = c(
      nchar(rna),
      str_count(rna, "A"),
      str_count(rna, "U"),
      str_count(rna, "G"),
      str_count(rna, "C"),
      ifelse(nchar(rna) > 0, round(100 * (str_count(rna, "G") + str_count(rna, "C")) / nchar(rna), 2), 0),
      str_count(dna, "-"),
      str_count(dna, "N"),
      paste(unique(str_extract_all(dna, "[^ATGCN-]")[[1]]), collapse = ", ") %>% na_if(""),
      paste(unique(str_extract_all(rna, "[^AUCG]")[[1]]), collapse = ", ") %>% na_if("")
    )
  ) %>% replace(is.na(.), "None")
}

generate_qc_report <- function(stats, checks, rna_seq) {
  report <- list()
  
  if("length" %in% checks) {
    len <- as.numeric(stats$Value[stats$Metric == "RNA Length"])
    status <- if(len > 0) "✅ Pass" else "❌ Fail (Empty sequence)"
    report <- c(report, paste("Length:", len, "-", status))
  }
  
  if("gc" %in% checks) {
    gc_val <- as.numeric(stats$Value[stats$Metric == "GC Content (%)"])
    status <- case_when(
      gc_val >= 60 ~ "⚠️ Warning (High GC)",
      gc_val <= 30 ~ "⚠️ Warning (Low GC)",
      TRUE ~ "✅ Pass"
    )
    report <- c(report, paste("GC Content:", gc_val, "% -", status))
  }
  
  if("gaps" %in% checks) {
    gaps <- as.numeric(stats$Value[stats$Metric == "Original Gaps"])
    status <- if(gaps > 0) "⚠️ Warning (Gaps present)" else "✅ Pass"
    report <- c(report, paste("Original Gaps:", gaps, "-", status))
  }
  
  if("homopolymer" %in% checks) {
    hp <- str_extract_all(rna_seq, "A{5,}|U{5,}|G{5,}|C{5,}")[[1]]
    status <- if(length(hp) > 0) paste("❌ Fail (Found", length(hp), ")") else "✅ Pass"
    report <- c(report, paste("Homopolymers:", length(hp), "-", status))
  }
  
  paste(report, collapse = "\n\n")
}

calculate_dinucleotide_freq <- function(sequence) {
  if(nchar(sequence) < 2) return(NULL)
  dinucleotides <- sapply(1:(nchar(sequence)-1), function(i) substr(sequence, i, i+1))
  freq_table <- table(dinucleotides)
  as.data.frame(freq_table) %>% 
    rename(Dinucleotide = dinucleotides, Frequency = Freq) %>%
    mutate(Percentage = round(Frequency / sum(Frequency) * 100, 2))
}

calculate_positional_distribution <- function(sequence) {
  bases <- strsplit(sequence, "")[[1]]
  positions <- 1:length(bases)
  data.frame(Position = positions, Base = bases)
}

calculate_sequence_complexity <- function(sequence, window_size = 20) {
  if(nchar(sequence) < window_size) return(NULL)
  positions <- seq(1, nchar(sequence) - window_size + 1)
  complexities <- sapply(positions, function(i) {
    window <- substr(sequence, i, i + window_size - 1)
    length(unique(strsplit(window, "")[[1]]))
  })
  data.frame(Position = positions, Complexity = complexities)
}

# UI ----
ui <- fluidPage(
  theme = shinytheme("cyborg"),
  useShinyFeedback(),
  
  tags$head(
    tags$style(HTML("
      .academic-footer {
        font-size: 0.8em;
        color: #888;
        border-top: 1px solid #444;
        margin-top: 20px;
        padding-top: 10px;
        line-height: 1.4;
      }
      .credential-highlight {
        color: #2CA02C;
        font-weight: 500;
      }
    "))
  ),
  
  titlePanel(
    tags$div(
      h1("NucleoConvert Analytics", style = "margin-bottom: 5px;"),
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      h3("🔄 Input Settings"),
      textAreaInput("dna_seq", "DNA Sequence:", rows = 6,
                    placeholder = "Enter raw DNA (A, T, G, C, N, -)"),
      checkboxInput("clean_seq", "Auto-clean (remove N's, gaps & invalid chars)", TRUE),
      actionButton("convert", "Process Sequence", class = "btn-primary"),
      hr(),
      
      h3("🧪 Quality Control"),
      checkboxGroupInput("qc_checks", "QC Checks:",
                         choices = c("Sequence Length" = "length",
                                     "GC Content" = "gc",
                                     "Original Gaps" = "gaps",
                                     "Homopolymers" = "homopolymer"),
                         selected = c("length", "gc")),
      hr(),
      
      h3("📥 Output Options"),
      downloadButton("download_fasta", "FASTA", class = "btn-success"),
      downloadButton("download_csv", "CSV Report", class = "btn-info"),
      hr(),
      
      tags$div(class = "academic-footer",
               tags$p(
                 class = "credential-highlight",
                 "Nihal Habib, PhD Candidate",
                 tags$br(),
                 tags$span(style = "color: #3C3B6E;", 
                           "Mohammed VI University of Sciences & Health (UM6SS)"),
                 tags$br(),
                 tags$span(style = "color: #CC9900;", 
                           "Fulbright Scholar at George Washington University"),
                 tags$br(),
                 tags$em("Casablanca, Morocco | Washington DC, USA")
               )
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("RNA Output", 
                 verbatimTextOutput("rna_seq"),
                 tags$hr(),
                 DTOutput("stats_table")),
        
        tabPanel("Analytics",
                 fluidRow(
                   column(6, plotlyOutput("gc_plot")),
                   column(6, plotlyOutput("nt_plot"))
                 ),
                 fluidRow(
                   column(6, plotlyOutput("dinucleotide_heatmap")),
                   column(6, plotlyOutput("positional_distribution"))
                 ),
                 fluidRow(
                   column(6, plotlyOutput("sequence_complexity")),
                   column(6, plotlyOutput("length_comparison"))
                 ),
                 fluidRow(
                   column(6, plotlyOutput("invalid_chars_plot")),
                   column(6, plotlyOutput("homopolymer_plot"))
                 )),
        
        tabPanel("RNA Analysis",
                 h3("Advanced RNA Analysis"),
                 tabsetPanel(
                   tabPanel("Codon Usage",
                            plotlyOutput("codon_usage_plot"),
                            DTOutput("codon_table")),
                   tabPanel("ORF Finder",
                            verbatimTextOutput("orf_results"),
                            plotlyOutput("orf_plot")),
                   tabPanel("Poly-A Tail",
                            verbatimTextOutput("polyA_analysis"))
                 )),
        
        tabPanel("QC Report", 
                 verbatimTextOutput("qc_report"),
                 tags$hr())
      )
    )
  )
)

# Server ----
server <- function(input, output, session) {
  rna <- reactiveVal("")
  stats_df <- reactiveVal()
  
  observeEvent(input$convert, {
    dna <- input$dna_seq
    
    if(nchar(dna) == 0) {
      showFeedbackDanger("dna_seq", text = "Input required!")
      return()
    } else {
      hideFeedback("dna_seq")
    }
    
    cleaned_dna <- clean_dna(dna, input$clean_seq)
    processed_rna <- str_replace_all(cleaned_dna, "T", "U")
    
    rna(processed_rna)
    stats_df(calculate_stats(dna, processed_rna, input$clean_seq))
  })
  
  # RNA Output
  output$rna_seq <- renderText({
    req(rna())
    validate(need(nchar(rna()) > 0, "No valid RNA produced - check input and cleaning settings"))
    rna()
  })
  
  # Statistics Table
  output$stats_table <- renderDT({
    req(stats_df())
    datatable(stats_df(), 
              options = list(dom = 't', pageLength = 10),
              rownames = FALSE) %>%
      formatStyle("Metric", fontWeight = "bold")
  })
  
  # QC Report
  output$qc_report <- renderText({
    req(stats_df(), input$qc_checks)
    generate_qc_report(stats_df(), input$qc_checks, rna())
  })
  
  # GC Content Plot
  output$gc_plot <- renderPlotly({
    req(rna())
    gc_pct <- stats_df()$Value[stats_df()$Metric == "GC Content (%)"] %>% as.numeric()
    
    plot_ly(
      values = c(gc_pct, 100 - gc_pct),
      labels = c("GC", "AU"),
      type = "pie",
      hole = 0.4,
      textinfo = "label+percent",
      marker = list(colors = c("#2ca02c", "#9467bd"))) %>%
      layout(title = "GC Content Distribution")
  })
  
  # Nucleotide Frequency Plot
  output$nt_plot <- renderPlotly({
    req(rna())
    nt_counts <- stats_df() %>%
      filter(Metric %in% c("A", "U", "G", "C")) %>%
      mutate(Value = as.numeric(Value))
    
    ggplotly(
      ggplot(nt_counts, aes(x = Metric, y = Value, fill = Metric)) +
        geom_col() +
        scale_fill_manual(values = c(A = "#1f77b4", U = "#ff7f0e", G = "#2ca02c", C = "#d62728")) +
        labs(title = "Nucleotide Frequency"))
  })
  
  # Dinucleotide Frequency Heatmap
  output$dinucleotide_heatmap <- renderPlotly({
    req(rna())
    df <- calculate_dinucleotide_freq(rna())
    validate(need(!is.null(df), "Sequence too short for dinucleotide analysis"))
    
    plot_ly(data = df, x = ~Dinucleotide, y = ~Percentage, type = "bar",
            color = ~Dinucleotide) %>%
      layout(title = "Dinucleotide Frequency",
             yaxis = list(title = "Percentage (%)"))
  })
  
  # Positional Nucleotide Distribution
  output$positional_distribution <- renderPlotly({
    req(rna())
    pos_df <- calculate_positional_distribution(rna())
    
    ggplotly(
      ggplot(pos_df, aes(x = Position, fill = Base)) +
        geom_bar(position = "fill") +
        scale_fill_manual(values = c(A = "#1f77b4", U = "#ff7f0e", G = "#2ca02c", C = "#d62728")) +
        labs(title = "Positional Nucleotide Distribution",
             y = "Proportion") +
        theme_minimal())
  })
  
  # Sequence Complexity Plot
  output$sequence_complexity <- renderPlotly({
    req(rna())
    complexity_df <- calculate_sequence_complexity(rna())
    validate(need(!is.null(complexity_df), "Sequence too short for complexity analysis"))
    
    ggplotly(
      ggplot(complexity_df, aes(x = Position, y = Complexity)) +
        geom_line(color = "#17becf") +
        labs(title = "Sequence Complexity (Unique Bases in 20bp Windows)") +
        theme_minimal())
  })
  
  # Length Comparison Plot
  output$length_comparison <- renderPlotly({
    req(stats_df())
    len_data <- data.frame(
      Type = c("Original DNA", "Processed RNA"),
      Length = c(nchar(input$dna_seq), nchar(rna())))
    
    plot_ly(len_data, x = ~Type, y = ~Length, type = "bar",
            color = ~Type, colors = c("#636efa", "#ef553b")) %>%
      layout(title = "Sequence Length Comparison",
             yaxis = list(title = "Length (bp)"))
  })
  
  # Invalid Characters Plot
  output$invalid_chars_plot <- renderPlotly({
    req(stats_df())
    invalid_chars <- stats_df()$Value[stats_df()$Metric == "Invalid in DNA"]
    if(invalid_chars == "None") return(plotly_empty())
    
    chars <- strsplit(invalid_chars, ", ")[[1]]
    freq <- table(chars)
    df <- data.frame(Character = names(freq), Count = as.numeric(freq))
    
    plot_ly(df, labels = ~Character, values = ~Count, type = "pie",
            hole = 0.4, textinfo = "label+percent") %>%
      layout(title = "Invalid Characters Distribution")
  })
  
  # Homopolymer Plot
  output$homopolymer_plot <- renderPlotly({
    req(rna())
    hp <- str_extract_all(rna(), "A{5,}|U{5,}|G{5,}|C{5,}")[[1]]
    
    if(length(hp) == 0) return(plotly_empty())
    
    hp_df <- data.frame(
      Base = str_sub(hp, 1, 1),
      Length = nchar(hp))
    
    ggplotly(
      ggplot(hp_df, aes(x = Base, y = Length, fill = Base)) +
        geom_boxplot() +
        labs(title = "Homopolymer Length Distribution"))
  })
  
  # Codon Usage Analysis
  output$codon_usage_plot <- renderPlotly({
    req(rna())
    df <- calculate_codon_usage(rna())
    validate(need(!is.null(df), "Sequence too short for codon analysis"))
    
    plot_ly(data = df, x = ~Codon, y = ~Percentage, type = "bar",
            color = ~Codon) %>%
      layout(title = "Codon Usage Frequency",
             yaxis = list(title = "Percentage (%)"))
  })
  
  output$codon_table <- renderDT({
    req(rna())
    df <- calculate_codon_usage(rna())
    validate(need(!is.null(df), "No codons found"))
    datatable(df, options = list(dom = 't', pageLength = 10))
  })
  
  # ORF Finder
  output$orf_results <- renderText({
    req(rna())
    orfs <- find_orfs(rna())
    if(is.null(orfs)) return("No ORFs found.")
    paste("Found", nrow(orfs), "ORFs:\n",
          paste("Start:", orfs$Start, "End:", orfs$End, "Length:", orfs$Length, collapse = "\n"))
  })
  
  output$orf_plot <- renderPlotly({
    req(rna())
    orfs <- find_orfs(rna())
    validate(need(!is.null(orfs), "No ORFs found"))
    
    plot_ly(data = orfs, x = ~Start, y = ~Length, type = "bar",
            color = ~Length) %>%
      layout(title = "ORF Length Distribution",
             xaxis = list(title = "Start Position"),
             yaxis = list(title = "Length (bp)"))
  })
  
  # Poly-A Tail Detection
  output$polyA_analysis <- renderText({
    req(rna())
    polyA <- detect_polyA_tail(rna())
    if(is.null(polyA)) return("No poly-A tail detected.")
    paste("Poly-A tail detected:", polyA$Tail_Length, "bases long.")
  })
  
  # Download Handlers
  output$download_fasta <- downloadHandler(
    filename = "rna_sequence.fasta",
    content = function(file) {
      writeLines(paste0(">Processed_RNA\n", rna()), file)
    }
  )
  
  output$download_csv <- downloadHandler(
    filename = "sequence_report.csv",
    content = function(file) {
      write.csv(stats_df(), file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
