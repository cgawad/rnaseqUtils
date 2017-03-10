#' Generate Minimum spanning trees from expression matrices
#'
#' Computes the minimum spanning trees from a names list of expression matrices.
#'
#' @return create_mst_igraphs_for_samples retuens a list of lists. The outer list corresponds to each sample, whereas the inner list contains the following elements:
#' 1) the mst igraph object, and 2) the node by gene value of n values that were calculated from the node_aggregation_function, 3 a data.frame listing the cluster id of each cell in the sample
#'
#' @export
create_mst_igraphs_for_samples <- function(named_filtered_cell_by_ensembl_mat_list = list(), n_components, k, dist_measure, p, node_aggregation_function)  {
  mst_objects_list <- lapply(names(named_filtered_cell_by_ensembl_mat_list), function(.sample)  {
    tmp_mat <- named_filtered_cell_by_ensembl_mat_list[[.sample]]
    message(paste('Processing sample', .sample))
    #target_k <- ceiling()
    if(!is.null(n_components))  {
      tmp_prcomp <- prcomp(tmp_mat)
      tmp_ld_mat <- tmp_prcomp$x[,1:n_components]
      rownames(tmp_ld_mat) <- rownames(tmp_mat)
    }
    else  {
      tmp_ld_mat <- tmp_mat
    }
    tmp_kmeans_obj <- kmeans(tmp_ld_mat, k[[.sample]], iter.max = 40)
    tmp_kmeans_nodes_dist <- as.matrix(dist(tmp_kmeans_obj$centers, dist_measure, p))
    tmp_graph <- graph.adjacency(tmp_kmeans_nodes_dist, mode = 'undirected', weighted = TRUE)
    tmp_graph_mst <- mst(tmp_graph)
    cluster_df <- data.frame(sample_name = names(tmp_kmeans_obj$cluster), cluster = tmp_kmeans_obj$cluster, stringsAsFactors = FALSE)
    #create mean_values_of_nodes
    cluster_to_ensembl_summary_mat <- rnaseqUtils:::get_cluster_by_ensembl_mat_of_expr_summaries(tmp_mat, cluster_df, node_aggregation_function)
    return(list(tmp_graph_mst, cluster_to_ensembl_summary_mat, cluster_df))
  })
}

#' calculate cluster summaries
#'
#' calculate summary for each cluster defined in sample_to_cluster_df and output summary results.
#'
#' @return get_cluster_by_ensembl_mat_of_expr_summaries retuens a cluster by ensembl matrix consisting of the summmary value for each cluster and ensembl combination.
get_cluster_by_ensembl_mat_of_expr_summaries <- function(sample_by_ensembl_mat, sample_to_cluster_df = data.frame(), summary_function)  {
  cluster_id_to_matrix_rows_list <- lapply(1:length(unique(sample_to_cluster_df$cluster)), function(.id) {
    which(rownames(sample_by_ensembl_mat) %in% sample_to_cluster_df$sample_name[sample_to_cluster_df$cluster == .id])
  })
  cluster_by_ensembl_mat <- apply(sample_by_ensembl_mat, 2, function(.col)  {
    summarized_col <- sapply(cluster_id_to_matrix_rows_list, function(.cluster_row_ids)  {
      summary_function(.col[.cluster_row_ids])
    })
    return(summarized_col)
  })
  return(cluster_by_ensembl_mat)
}
