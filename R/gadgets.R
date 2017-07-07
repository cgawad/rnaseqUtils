

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


#' View expression of genes on 2D plots.
#'
#' @param ggplot_df a data.frame that must contain the following column names and associated types: sample_name, character; D1, numeric; D2, numeric
#' @param expr_mat a matrix with rownames corresponding to sample_name and with column names corresponding to ENSEMBL IDs. The matrix must be numeric
#' @param gene_selector_df a data.frame with an ENSEMBL column corresponding to the colnames of expr_mat and containing other informative columns such as GENENAME, GENE_TYPE, etc...
#' @return Returns null.
#' @export
explore_expression <- function(ggplot_df, expr_mat, gene_selector_df) {
  library(shiny)
  ui <- fluidPage(
    #gadgetTitleBar("Expression Display"),
    mainPanel(
      dataTableOutput('gene_selector'),
      plotOutput('vis')
    )
  )

  server <- function(input, output, session) {
    # Define reactive expressions, outputs, etc.

    output$gene_selector <- renderDataTable(gene_selector_df, selection = 'single')
    # When the Done button is clicked, return a value
    output$vis <- renderPlot({
      print(input$gene_selector)
      if(length(input$gene_selector_rows_selected) == 0)  {
        ggplot(ggplot_df, aes(x=D1, y=D2, shape = sample_replicate, color = sample_date)) + geom_point()
      }
      else  {
        final_df <- dplyr::left_join(ggplot_df,rnaseqUtils::get_df_from_named_char_vector(expr_mat[,gene_selector_df$ENSEMBL[input$gene_selector_rows_selected]], c('sample_name', 'expr')))
        print(head(final_df))
        ggplot(final_df, aes(x=D1, y=D2, color = expr)) + geom_point() + scale_color_gradient2(low='white', high = 'darkred')
      }
    })


  }
  runGadget(ui, server)
}

#' Explore a dataset with monocle
#'
#' This gadget uses monocle to plot a DDRTree and visualization of a branched
#' heatmap. The branch point and root of the tree can be changed interactively.
#'
#' @param expression_mat: This is a CellxGene matrix that will be used to run
#'   monocle.
#' @param max_cell_counts: Moncole can choke on large datasets, so this puts a
#'   cap on then number of cells from expression_mat that will be used to
#'   generate a 2D structure.
#' @param metadata_df: This is a data.frame that contains metatadata for the
#'   cells in expression_mats. It is used for coloring the monocle plots by
#'   additinal variables. It must at the very least have a sample_name column
#'   that matches the rownames of expression_mat.
#' @param species: this should be pretty obvious. Currently, 'mmusculus' and
#'   'hsapiens' are valid
#' @param normalizing_scale_factor: This is to un-normalize the data. This has
#'   to be done because of the different normalization procedure that monocle
#'   uses. This factor is just the median umi count that was used in the
#'   normalization of expression_mat. It can be found as an element of the list
#'   returned from get_normalized_expression_matrix.
#' @param cell_umi_counts: This is a named character vector of cell UMI counts,
#'   again required for un-normalizing the data in expression_mat. It can also
#'   be found as an element of the list returned from
#'   get_normalized_expression_matrix.
#'
#' @return This function launches an app that allows interactive exploration of
#'   a dataset. There is no object returned.
#'
#' @export
run_interactive_monocle_gadget <- function(expression_mat = matrix(), max_cell_counts = 1000, metadata_df = data.frame(), species = c('mmusculus', 'hsapiens'), normalizing_scale_factor = NULL, cell_umi_counts = NULL)  {
  library(shiny)
  library(monocle)
  if(is.null(normalizing_scale_factor))  {
    stop("Cannot 'un-normalize' (not a real word, I know) the matrix to the original count matrix without knowing the normalizing_scale_factor")
  }
  if(is.null(cell_umi_counts))  {
    stop("Cannot 'un-normalize' (not a real word, I know) the matrix to the original count matrix without knowing the original UMI counts")
  }
  if(length(unique(as.vector(foobar$top_1000_mat %% 1))) == 1)  {
    stop("expression_mat should be normalized on input! It is un-normalized during execution of this gadget")
  }
  if(is.null(names(cell_umi_counts)))  {
    stop("I pity the foo who doesn't use a named character vector for cell_umi_counts!")
  }

  expression_mat <- ((2^(expression_mat) -1)/normalizing_scale_factor) * cell_umi_counts[rownames(expression_mat)]

  if(species == 'mmusculus')  {
    f_data_df <- Biobase::AnnotatedDataFrame({dplyr::left_join(data.frame(GENEID = colnames(expression_mat), stringsAsFactors = FALSE), AnnotationDbi::select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79, keys = colnames(expression_mat), keytype = 'GENEID', columns = c('GENEID', 'GENENAME'))) %>% dplyr::mutate(GENENAME = ifelse(is.na(GENENAME), GENEID, GENENAME)) %>% dplyr::rename(ENSEMBL = GENEID, gene_short_name = GENENAME) %>% tibble::column_to_rownames('ENSEMBL')})
  }
  else if(species == 'hsapiens') {
    f_data_df <- Biobase::AnnotatedDataFrame({dplyr::left_join(data.frame(GENEID = colnames(expression_mat), stringsAsFactors = FALSE), AnnotationDbi::select(EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79, keys = colnames(expression_mat), keytype = 'GENEID', columns = c('GENEID', 'GENENAME'))) %>% dplyr::mutate(GENENAME = ifelse(is.na(GENENAME), GENEID, GENENAME)) %>% dplyr::rename(ENSEMBL = GENEID, gene_short_name = GENENAME) %>% tibble::column_to_rownames('ENSEMBL')})
  }
  else  {
    stop(paste(species, 'has not been handled yet'))
  }

  p_data_df <- Biobase::AnnotatedDataFrame(
    {data.frame(sample_name = rownames(expression_mat), stringsAsFactors = FALSE) %>% dplyr::left_join(., metadata_df) %>% tibble::column_to_rownames('sample_name')}# %>% dplyr::select(-sample_name)}
    )

  temp_cds <- monocle::newCellDataSet(as(t(expression_mat), 'sparseMatrix'), phenoData = p_data_df, featureData = f_data_df, expressionFamily = VGAM::negbinomial.size())

  temp_cds <- BiocGenerics::estimateSizeFactors(temp_cds)
  temp_cds <- BiocGenerics::estimateDispersions(temp_cds)

  if(nrow(expression_mat) > max_cell_counts)  {
    temp_cds <- monocle::reduceDimension(temp_cds[,sample(1:ncol(temp_cds), max_cell_counts, replace = FALSE)], max_components=2)
  }
  else  {
    temp_cds <- monocle::reduceDimension(temp_cds, max_components=2)
  }
  temp_cds <- monocle::orderCells(temp_cds)


  cols_to_select_list <- as.list(c('Pseudotime', 'State', colnames(metadata_df)))
  #cols_to_select = 1:length(cols_to_select_names)
  #names(cols_to_select) <- cols_to_select_names

  ui <- fluidPage(
    #gadgetTitleBar("Expression Display"),
    titlePanel(title = 'Monocle Explorer'),
    sidebarLayout(
      sidebarPanel(
        selectInput('color_by',
                    label = 'color_by',
                    choices = cols_to_select_list,
                    selected = 'Pseudotime'),
        selectInput('root_state',
                    label = 'root state',
                    choices = as.list(unique(pData(temp_cds)$State)),
                    selected = 1),
        numericInput('branch_point',
                    label = 'branch_point',
                    value = 1),
        selectInput("n_heatmap_entries",
                    label = "# heatmap entries",
                    choices = as.list(c(25,50,75,100,125,150,175,200)),
                    selected = 25)
      ),
      mainPanel(
        plotOutput('monocle_plot'),
        plotOutput('monocle_heatmap')
      )
    )
  )

  server <- function(input, output, session) {
    new_root <- reactive({
      temp_cds <<- monocle::orderCells(temp_cds, root_state = as.integer(input$root_state))
    })

    output$monocle_plot <- renderPlot({
      new_root()#temp_cds <- monocle::orderCells(temp_cds, root_state = as.integer(input$root_state))
      monocle::plot_cell_trajectory(temp_cds, color_by=input$color_by)
    })
    output$monocle_heatmap <- renderPlot({
      new_root()
      target_genes <- as.character({monocle::BEAM(temp_cds, branch_point = as.integer(input$branch_point), relative_expr = FALSE) %>% dplyr::arrange(qval)}$fd[1:as.integer(input$n_heatmap_entries)])
      target_ensembls <- rownames(fData(temp_cds))[fData(temp_cds)$gene_short_name %in% target_genes]
      #print(target_genes)
      monocle::plot_genes_branched_heatmap(temp_cds[target_ensembls,], branch_point = as.integer(input$branch_point), use_gene_short_name = TRUE, show_rownames = TRUE)
    })

  }

  runGadget(ui, server)
}
