library(shiny)
library(shinythemes)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(biomaRt)

# Define UI
ui <- fluidPage(
  theme = shinytheme("cerulean"),
  titlePanel("Gene Expression Explorer"),
  sidebarLayout(
    sidebarPanel(
      textInput("gene", "Enter Gene Name or Ensembl ID:", value = "SRRM4"),
      actionButton("plot", "Plot Gene Counts")
    ),
    mainPanel(
      plotOutput("genePlot")
    )
  )
)

# Define server logic
server <- function(input, output) {
  # Set correct working directory relative to the omics_coursework directory
  setwd("../input_data/p2")
  
  # Load data
  count_data <- read.table("GSE64018_countlevel_12asd_12ctl_edited.txt", header = TRUE, row.names = 1)
  metadata <- read.table("metadata.txt", header = TRUE, row.names = 1)
  metadata <- metadata %>% mutate(across(where(is.character), as.factor))
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~ diagnosis)
  dds$diagnosis <- relevel(dds$diagnosis, ref = "CTL")
  dds <- dds[rowSums(counts(dds)) > 10, ]
  dds <- DESeq(dds)
  
  observeEvent(input$plot, {
    gene_id <- input$gene
    res <- results(dds)
    res <- as.data.frame(res)
    res <- res %>% drop_na(padj)
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    g <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), filters = "ensembl_gene_id", mart = ensembl, values = rownames(res))
    res$ensembl_id <- rownames(res)
    res <- merge(res, g, by.x = "ensembl_id", by.y = "ensembl_gene_id")
    colnames(res)[10] <- "gene_name"
    
    if (gene_id %in% res$gene_name | gene_id %in% res$ensembl_id) {
      srrm4_id <- ifelse(gene_id %in% res$gene_name, res[res$gene_name == gene_id, "ensembl_id"], gene_id)
      srrm4_data <- plotCounts(dds, gene = srrm4_id, intgroup = "diagnosis", returnData = TRUE)
      output$genePlot <- renderPlot({
        ggplot(srrm4_data, aes(diagnosis, count)) +
          geom_boxplot() +
          geom_jitter(aes(col = diagnosis), alpha = 0.4, width = 0.1) +
          labs(x = "Diagnosis", y = "Normalized count") +
          guides(col = "none") +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5))
      })
    } else {
      output$genePlot <- renderPlot({
        plot.new()
        text(0.5, 0.5, "Gene not found", cex = 2, col = "red")
      })
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)
