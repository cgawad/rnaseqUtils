#' Generate several plots and objects for dataset exploration
#'
#' This function examines several facets of the dataset, inclusing filtering
#' cells based on metadata, dimension reduction, cluster analysis, and
#' visualization
#'
#' @param cell_dataframe: cell_dataframe consists of sample names of each cell
#'   to be examined, as well as all other data that might be used for
#'   visualizing the resulting cells. This data.frame must have a character
#'   column named 'sample_name', which must match the row names of
#'   normalized_expression_matrix
#'
#' @param normalized_expression_matrix: this must be a CxG matrix, with C cells
#'   and G genes, The C cells are the rownamews of the matrix and must match the
#'   sample_name column from cell_dataframe. The column names must be ENSEMBL
#'   IDs. The entries of the matrix must be the normalized since they will be
#'   used in tsne
#'
#' @param .species: currently either 'mmusculus' or 'hsapiens'
#'
#' @param color_var: variable name to color cells by in overview ggplot images.
#'   This must correspond to a column name in cell_data_dataframe.
#'
#' @return This function returns a list with the following elements: \itemize{
#'   \item \strong{top_1000_mat}: a matrix of the top 1000 overdispersed genes
#'   \item \strong{tsne_results}: a list of Rtsne outputs from the filtered
#'   matrix, each from a differnet perplexity (15,30,45,60) \item
#'   \strong{n50_enriched_genes_by_cluster}: a dataframe showing the top 50
#'   enriched genes in each cluster. Enrichment is based on proportion of cells
#'   expressing a gene \item \strong{updated_metadata_df}: the original
#'   cell_dataframe annotated with cluster_id \item
#'   \strong{clusters_scatterplot}: a scatterplot coloring the cells by
#'   cluster_id\item \strong{important_genes_heatmap}: a heatmap plot showing
#'   the enriched genes in each cluster \item \strong{overview_plots}: a list of
#'   plots corresponding to the tsne runs that are coloredc by @param color_by }

#' @export
examine_dataset <- function(cell_data_dataframe = data.frame(), normalized_expression_matrix = matrix(), .species = 'mmusculus', color_var = character())  {
  results_list = list()
  normalized_expression_matrix <- normalized_expression_matrix[cell_data_dataframe$sample_name,]
  top_1000_dispersed_mat <- get_dispersion(normalized_expression_matrix)
  results_list[['top_1000_mat']] <- top_1000_dispersed_mat
  top_1000_dispersed_prcomp <- prcomp(top_1000_dispersed_mat, scale. = TRUE)
  top_1000_dispersed_rtsne_list <- lapply(c(15,30,45,60), function(.perplexity)  {Rtsne::Rtsne(top_1000_dispersed_prcomp$x[,1:50], perplexity = .perplexity, pca = FALSE)})
  results_list[['tsne_results']] <- top_1000_dispersed_rtsne_list
  tree_splits_list <- get_treecuts(top_1000_dispersed_mat)

  sample_name_cluster_id_df <- get_sample_name_to_cluster_id_df(tree_splits_list[[2]][[3]], normalized_expression_matrix)
  n50_enriched_genes_by_cluster_df <- get_enriched_genes_in_each_cluster(sample_name_cluster_id_df, expression_mat = normalized_expression_matrix, top_n = 50, species = .species)
  results_list[['n50_enriched_genes_by_cluster']] <- n50_enriched_genes_by_cluster_df

  results_list[['updated_metadata_df']] <- dplyr::left_join(cell_data_dataframe, sample_name_cluster_id_df)
  clusters_scatterplot <- get_clustered_scatterplot(top_1000_dispersed_rtsne_list[[4]]$Y, sample_name_cluster_id_df)
  results_list[['clusters_scatterplot']] <- clusters_scatterplot


  prop_df <- get_prop_expressed_df(sample_to_cluster_df = sample_name_cluster_id_df, expression_mat = as.matrix(normalized_expression_matrix[,unique(n50_enriched_genes_by_cluster_df$ENSEMBL)]), species = .species)

  important_genes_heatmap_ggplot <- plot_ordered_heatmap(prop_df)
  results_list[['important_genes_heatmap']] <- important_genes_heatmap_ggplot

  temp_visualization_list <- lapply(top_1000_dispersed_rtsne_list, function(tsne_result)  {
    overview_df <- data.frame(tsne_result$Y, sample_name = rownames(normalized_expression_matrix), stringsAsFactors = FALSE) %>% dplyr::rename(D1 = X1, D2 = X2) %>% dplyr::left_join(., cell_data_dataframe)
    ggplot(overview_df, aes_string(x = 'D1', y = 'D2', color = color_var)) + geom_point()
  })
  results_list[['overview_plots']] <- temp_visualization_list

  return(results_list)
}

#' Process cellranger output
#'
#' Read cellranger output directory, filter the data, and return normalized
#' expression matrix
#'
#' @param run_to_path_df: A dataframe with two columns (1) run_name (character): The name
#' of the cellranger run, such as patient_id, sample_date, etc...
#' (2) cellranger_path (character): The path to the toplevel cellranger output
#'
#' @param .genome: a character, typically 'mm10', or 'hg19'.
#'
#' @param umi_limits: vector of length 2, with lower and upper bounds for umi_counts.
#'
#' @return This function reads a series of cellranger output directories. It merges the
#' results together into a large matrix anbd then filters out cells with UMI
#' counts outside of the range specified in umi_limits. Mitochondrial and
#' ribosomal protein-coding genes are then removed from the matrix, as are
#' ENSEMBL IDs with no expression across the remaining cells. Finally, the data
#' is scaled using a global scaling factor (total UMIs per cell) and
#' log2-transformed.
#'
#' @export
get_normalized_expression_matrix  <- function(run_to_path_df = data.frame(), .genome = c('mm10', 'hg19'), umi_limits = c(3e3,Inf))  {
  output_list = list()
  mat_list <- lapply(1:nrow(run_to_path_df), function(.ind)  {
    #sample_name <- sub("PAN_", "", sub('_[^_]+', '', basename(.filename)))
    sample_obj <- cellrangerRkit::load_cellranger_matrix(run_to_path_df$cellranger_path[.ind], genome = .genome)
    sample_mat <- as.matrix(exprs(sample_obj))
    print(class(sample_mat))
    colnames(sample_mat) <- paste(run_to_path_df$run_name[.ind], colnames(sample_mat), sep = "_")
    return(sample_mat)
  })
  print(sapply(mat_list, class))
  giant_mat <- t(do.call(cbind, mat_list))
  rm(mat_list)
  print(class(giant_mat))

  #Remove ribosomal and mt proteins
  if(.genome == 'mm10')  {
    giant_norpmt_mat <- remove_ribosomal_proteins_and_mitochondrial_genes_from_matrix(cell_by_ensembl_mat = giant_mat, species = 'mmusculus')
  }
  else if(.genome == 'hg19')  {
    giant_norpmt_mat <- remove_ribosomal_proteins_and_mitochondrial_genes_from_matrix(cell_by_ensembl_mat = giant_mat, species = 'hsapiens')
  }
  else  {
    stop(paste(.genome, 'not handled yet.'))
  }

  #Remove cells outside the umi count limits defined by umi_limits
  umi_counts <- rowSums(giant_norpmt_mat)
  filtered_mat <- giant_norpmt_mat[(umi_counts <= umi_limits[2]) & (umi_counts >= umi_limits[1]),]
  rm(umi_counts)

  #Remove non-expressed genes
  filtered_mat <- filtered_mat[,colSums(filtered_mat) > 0]
  filtered_mat <- log2((filtered_mat / rowSums(filtered_mat)) * median(rowSums(filtered_mat)) + 1)

  umi_counts <- rowSums(filtered_mat)
  names(umi_counts) <- rownames(filtered_mat)

  return(list(normalized_matrix = filtered_mat, umi_counts = umi_counts, normalization_scale_factor = median(umi_counts)))

}
