#' @export
launch_rnaseq_data_explorer_gadget <- function(dataset_level_df = data.frame(), cluster_level_df = data.frame(), ENSEMBL_to_color_var_list = list()) {

  ui <- shinyUI(fluidPage(

    # Application title
    title = "Gene Expression Visualization in 3 Dimensions",

    fluidRow(
      column(6,
             h3("Expression Levels in Cells", align = 'center'),
             wellPanel(
               plotlyOutput("plot", width = '100%', height = "100%")
             ),
             selectInput(
               "select", label = h3("Cell level Descriptor"), choices = as.list(c('gene', names(dataset_level_df)[!grepl("^D[01-9]", names(dataset_level_df))][-1])), selected = 3
             )
      )
    ),
    fluidRow(
      h1("Gene", align = 'center'),
      DT::dataTableOutput('dt')#, align = 'center')
    )
  ))

  server <- shinyServer(function(input, output) {
    output$plot <- renderPlotly({
      if(length(input$dt_rows_selected) > 0 && input$select == 'gene')  {
        print(paste("Looking for", cluster_level_df[input$dt_rows_selected, 'ENSEMBL']))
        new_coloring_variable <- ENSEMBL_to_color_var_list[[cluster_level_df[input$dt_rows_selected, 'ENSEMBL']]]
      }
      else if(input$select != 'gene') {
        new_coloring_variable <- dataset_level_df[,input$select]
      }
      else if(length(input$dt_rows_selected) == 0 && input$select == 'gene')  {
        new_coloring_variable <- NULL
      }
      dims_inds <- grep("^D[123]", names(dataset_level_df))
      if(length(dims_inds) == 3)  {
        p<-plot_ly(x = dataset_level_df$D1, y = dataset_level_df$D2, z = dataset_level_df$D3, color = new_coloring_variable, text = paste(dataset_level_df$sample_name, cluster_level_df[input$dt_rows_selected, 'SYMBOL'], sep = ","), type = 'scatter3d', mode = 'markers')
      }
      else if(length(dims_inds) == 2) {
        p<-plot_ly(x = dataset_level_df$D1, y = dataset_level_df$D2, color = new_coloring_variable, text = paste(dataset_level_df$sample_name, cluster_level_df[input$dt_rows_selected, 'SYMBOL'], sep = ","), type = 'scatter', mode = 'markers')
      }
      else  {
        stop(paste('dataset_level_df needs do have either D1,D2 or D1,D2,D3 as columns'))
      }
      layout(p)
    })
    output$dt <- DT::renderDataTable(cluster_level_df, server =  TRUE, selection = 'single', filter = 'bottom', options = list(pageLength = 5, autoWidth = TRUE))
  })

  runGadget(ui, server)
}

#' @export
subsample_cell_ranger_data_gadget <- function(cell_by_gene_count_mat = matrix())  {
  n_genes <- apply(cell_by_gene_count_mat, 1, function(.x)  {
    sum(.x > 0)
  })
  n_reads <- apply(cell_by_gene_count_mat, 1, function(.x)  {
    sum(.x)
  })
  n_reads_and_genes_df <- data.frame(sample_name = rownames(cell_by_gene_count_mat), n_genes = n_genes, n_reads = n_reads, stringsAsFactors = FALSE)
  n_genes_plot <- ggplot2::ggplot(n_reads_and_genes_df, aes(n_genes)) + geom_histogram()
  n_reads_plot <- ggplot2::ggplot(n_reads_and_genes_df, aes(n_reads)) + geom_histogram()
  gridExtra::grid.arrange(n_genes_plot, n_reads_plot, ncol = 1)
}


