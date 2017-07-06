#' Generate several plots and objects for dataset exploration
#'
#'
#' cell_dataframe: cell_dataframe consists of sample names of each cell to be
#' examined, as well as all other data that might be used for visualizing the
#' resulting cells. This data.frame must have a character column named
#' 'sample_name', which must match the row names of normalized_expression_matrix
#'
#' normalized_expression_matrix: this must be a CxG matrix, with C cells and G
#' genes, The C cells are the rownamews of the matrix and must match the
#' sample_name column from cell_dataframe. The column names must be ENSEMBL IDs.
#' The entries of the matrix must be the normalized since they will be used in tsne
#'
#' color_var_by: variable name to color cells by in overview ggplot images. This
#' must correspond to a column name in cell_data_dataframe.
#'

#' @export
examine_dataset <- function(cell_data_dataframe = data.frame(), normalized_expression_matrix = matrix(), .species = 'mmusculus', color_var = character())  {
  results_list = list()
  top_1000_dispersed_mat <- get_dispersion(normalized_expression_matrix)
  top_1000_dispersed_prcomp <- prcomp(top_1000_dispersed_mat, scale. = TRUE)
  top_1000_dispersed_rtsne_list <- lapply(c(15,30,45,60), function(.perplexity)  {Rtsne::Rtsne(top_1000_dispersed_prcomp$x[,1:50], perplexity = .perplexity, pca = FALSE)})
  results_list[['tsne_results']] <- top_1000_dispersed_rtsne_list
  tree_splits_list <- get_treecuts(top_1000_dispersed_mat)

  sample_name_cluster_id_df <- get_sample_name_to_cluster_id_df(tree_splits_list[[2]][[3]], normalized_expression_matrix)
  n50_enriched_genes_by_cluster_df <- get_enriched_genes_in_each_cluster(sample_name_cluster_id_df, expression_mat = normalized_expression_matrix, top_n = 50, species = .species)
  results_list[['n50_enriched_genes_by_cluster']] <- n50_enriched_genes_by_cluster_df

  clusters_scatterplot <- get_clustered_scatterplot(top_1000_dispersed_rtsne_list[[4]]$Y, sample_name_cluster_id_df)
  results_list[['clusters_scatterplot']] <- clusters_scatterplot


  prop_df <- get_prop_expressed_df(sample_to_cluster_df = sample_name_cluster_id_df, expression_mat = as.matrix(normalized_expression_matrix[,unique(n10_enriched_genes_by_cluster_df$ENSEMBL)]), species = .species)

  important_genes_heatmap_ggplot <- plot_ordered_heatmap(prop_df)
  results_list[['important_genes_heatmap']] <- important_genes_heatmap_ggplot

  temp_visualization_list <- lapply(top_1000_dispersed_rtsne_list, function(tsne_result)  {
    overview_df <- data.frame(tsne_result$Y, sample_name = rownames(normalized_expression_matrix), stringsAsFactors = FALSE) %>% dplyr::rename(D1 = X1, D2 = X2) %>% dplyr::left_join(., cell_data_dataframe)
    ggplot(overview_df, aes_string(x = 'D1', y = 'D2', color = color_var)) + geom_point()
  })
  results_list[['overview_plots']] <- temp_visualization_list

  return(results_list)
}
