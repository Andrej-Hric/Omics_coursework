library(shiny)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(DT)
library(dplyr)
library(tidyr)
library(biomaRt)
library(shinythemes)

set.seed(123)  # Ensure reproducibility

# Load data from GitHub
load_data <- function() {
  count_data_url <- "https://raw.githubusercontent.com/Andrej-Hric/Omics_coursework/main/input_data/p2/GSE64018_countlevel_12asd_12ctl_edited.txt"
  metadata_url <- "https://raw.githubusercontent.com/Andrej-Hric/Omics_coursework/main/input_data/p2/metadata.txt"
  
  count_data <- read.table(count_data_url, header = TRUE, row.names = 1)
  metadata <- read.table(metadata_url, header = TRUE, row.names = 1)
  
  metadata <- metadata %>% mutate(across(where(is.character), as.factor))
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~ diagnosis)
  dds$diagnosis <- relevel(dds$diagnosis, ref = "CTL")
  dds <- dds[rowSums(counts(dds)) > 10, ]
  dds <- DESeq(dds)
  
  res <- results(dds)
  df.plot <- as.data.frame(res) %>% drop_na(padj)
  df.plot <- df.plot[order(abs(df.plot$log2FoldChange), decreasing = TRUE), ]
  df.plot <- mutate(df.plot, sig = ifelse(df.plot$padj < 0.05, "FDR<0.05", "Not sig"))
  df.plot <- cbind(df.plot, top10 = "no")
  top10.index <- rownames(df.plot[df.plot$sig == "FDR<0.05", ][1:10, ])
  df.plot[top10.index, "top10"] <- "yes"
  
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  g <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), filters = "ensembl_gene_id", mart = ensembl, values = rownames(df.plot))
  df.plot$ensembl_id <- rownames(df.plot)
  df.plot.gene <- merge(df.plot, g, by.x = "ensembl_id", by.y = "ensembl_gene_id")
  colnames(df.plot.gene)[10] <- "gene_name"
  
  df.plot.gene$gene_label <- ifelse(df.plot.gene$gene_name == "", df.plot.gene$ensembl_id, df.plot.gene$gene_name)
  
  list(dds = dds, df.plot.gene = df.plot.gene)
}

data <- load_data()
dds_DGE <- data$dds
df.plot.gene <- data$df.plot.gene

ui <- fluidPage(
  titlePanel("Gene Expression Explorer"),
  sidebarLayout(
    sidebarPanel(
      textInput("gene_id", "Enter Gene ID or Name:", value = ""),
      actionButton("search", "Search")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Volcano Plot", plotOutput("volcano")),
        tabPanel("Box Plot", plotOutput("boxplot", height = "600px")),  # Make boxplot taller
        tabPanel("About", tags$div(id = "about", 
                                   tags$p("This Shiny app allows users to explore gene expression data from a study on autistic and neurotypical brains."),
                                   tags$p("Users can search for specific genes and view their expression counts across different samples."),
                                   tags$p("The app includes a volcano plot for visualizing differential expression and a box plot for viewing gene counts."),
                                   tags$p("This app is created as part of the coursework for the Omics module in the MSc Bioinformatics program at Birkbeck, University of London.")
        ))
      ),
      fluidRow(
        column(
          12, DT::dataTableOutput("table")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  # Render the table
  output$table <- DT::renderDataTable({
    df.plot.gene %>%
      mutate(across(where(is.numeric), ~ round(., 3)))
  }, options = list(pageLength = 5, selection = "single"))
  
  # Reactive function for gene selection
  df.click <- reactive({
    req(input$search)
    isolate({
      gene_id <- input$gene_id
      if (gene_id %in% df.plot.gene$ensembl_id | gene_id %in% df.plot.gene$gene_name) {
        df.plot.gene[df.plot.gene$ensembl_id == gene_id | df.plot.gene$gene_name == gene_id, ]
      } else {
        df.plot.gene[df.plot.gene$gene_name == "SRRM4", ]
      }
    })
  })
  
  get.df <- reactive({
    selected <- input$table_rows_selected
    if (length(selected)) {
      df.plot.gene[selected, ]
    } else {
      df.click()
    }
  })
  
  # Volcano plot
  output$volcano <- renderPlot({
    p <- ggplot(df.plot.gene, aes(log2FoldChange, -log10(padj), col=sig)) +
      geom_point(size=0.5) + 
      scale_color_manual(values = c("orange", "gray")) +
      theme_bw(base_size = 16)
    p$labels$colour <- "Significance"
    
    df.temp <- get.df()
    p <- p + geom_point(data = df.temp, size=1.5, col="magenta", shape = 23) +
      geom_text_repel(
        data = df.temp,
        aes(label = gene_label),
        col = "magenta"
      )
    p
  })
  
  # Boxplot
  output$boxplot <- renderPlot({
    df.temp <- get.df()
    gene_ensembl <- tail(df.temp, 1)$ensembl_id
    gene_title <- tail(df.temp, 1)$gene_label
    
    gene_data <- plotCounts(dds_DGE, gene_ensembl, "diagnosis", returnData = TRUE)
    
    q <- ggplot(gene_data, aes(diagnosis, count)) +
      geom_boxplot() +
      geom_jitter(aes(col = diagnosis), alpha = 0.4, width = 0.1) +
      labs(
        title = gene_title,
        x = "Diagnosis",
        y = "Normalized count"
      ) +
      guides(col = "none") +
      theme_bw(base_size = 16) +
      theme(plot.title = element_text(hjust = 0.5))
    q
  }, height = 600)  # Make boxplot taller
}

shinyApp(ui = ui, server = server)