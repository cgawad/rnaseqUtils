# 'Load 10X matrix files
#'
#' Loads the matrix files from a 10X run
#'
#' @return  get_cellranger_matrices returns a list of count matrices corresponding to the matrix files of \code{filenames}. The names of the list are the basenames of \code{filenames}

#' @export
get_cellranger_matrices <- function(filenames = character(), sample_name_prefix = NULL, genome_build = 'hg19')  {
  matrix_list <- lapply(filenames, function(.filename)  {
    temp_gbm <- cellrangerRkit::load_cellranger_matrix(.filename, genome = genome_build)
    temp_gbm_mat <- as.matrix(exprs(temp_gbm))
    #gbm_gt1cell_mat <- gbm_mat[apply(as.matrix(gbm_mat), 1, function(.x)  {sum(.x > 1)}) > 1,]
    # gbm_gt1cell_df <- as.data.frame(as.matrix(gbm_gt1cell_mat), optional = TRUE, stringsAsFactors = FALSE)
    # gbm_gt1cell_df <- mutate(gbm_gt1cell_df, ENSEMBL = row.names(gbm_gt1cell_df))
    # gbm_gt1cell_df <-  reshape2::melt(gbm_gt1cell_df, id.vars = 'ENSEMBL', variable.name = 'sample_name', value.name = 'count')
    # gbm_gt1cell_df <- dplyr::rename(gbm_gt1cell_df, sample_name = variable, count = value)
    # gbm_gt1cell_df <- left_join(gbm_gt1cell_df, gene_and_umi_counts_per_cell_df) %>% mutate(log_tpht = log2((count/umi_count)*1e5 + 1))

    #get_matrix_from_df <- function(input_df, cast_formula = formula("sample_name ~ ENSEMBL"), value.var = 'log_tpht')  {
    #  temp_mat <- reshape2::acast(input_df, cast_formula, value.var = value.var)
    #  return(temp_mat)
    #}
    return(temp_gbm_mat)
  })
  names(matrix_list) <- basename(filenames)
  if(length(sample_name_prefix) > 0)  {
    for(.ind in 1:length(sample_name_prefix))  {
      colnames(matrix_list[[.ind]]) <- paste(sample_name_prefix[.ind], colnames(matrix_list[[.ind]]), sep = "_")
      #return(matrix_list)
    }
  }
  return(matrix_list)
}

# 'Load 10X matrix files
#'
#' Loads and merges the matrix files from a 10X run
#'
#' @return read_merge_and_return_sample_matrices returns a count matrices corresponding to the matrix files of \code{filenames}. The matrices of the original samples are merged with the rows corresponding to samples and the columns corresponding to ensembl ids.

#' @export
read_merge_and_return_sample_matrices <- function(filenames = character(), sample_name_prefix = NULL, genome_build = 'hg19')  {
  expr_mats <- get_cellranger_matrices(filenames, sample_name_prefix = sample_name_prefix, genome_build = genome_build)
  shared_rownames <-  unique(unlist(lapply(expr_mats, rownames)))
  raw_expr_mat <- t(do.call(cbind, lapply(expr_mats, function(.x)  {return(.x[shared_rownames,])})))
  rm(expr_mats, filenames, shared_rownames)
  return(raw_expr_mat)
}



#converts a named charcter vector into a 2-column data.frame with the column names specified by col_names
#' @export
get_df_from_named_char_vector <- function(char_vec = character(), col_names = character())  {
  temp_df <- as.data.frame(char_vec)
  names(temp_df) <- col_names[2]
  temp_df[col_names[1]] <- row.names(temp_df)
  temp_df <- temp_df[,2:1]
  row.names(temp_df) <- NULL
  return(temp_df)
}

#subsample the count matrix such that each cell (row) has threshold_level reads. Currently, cells with fewer reads are discarded, although in the future they will be upsampled based on similar cell read/UMI counts
#' @export
rarefy_count_matrix <- function(count_matrix = matrix(), threshold_level = numeric())  {
  counts_per_cell <- apply(count_matrix, 1, sum)
  count_sub_mat <- count_matrix[counts_per_cell >= threshold_level,]
  gene_by_sample_mat <- sapply(rownames(count_sub_mat), function(.sample_name)  {
    sample_prob <-  threshold_level / counts_per_cell[.sample_name]
    return(sapply(count_sub_mat[.sample_name,], function(.col)  {
      return(rbinom(rep(1, length(.col)), .col, sample_prob))
    }))
  })
  return(t(gene_by_sample_mat))
}

#This function takes a data.frame with numeric columns, computes the aic for clusterings with 1:25 clusters, identifies the cluster number such that the following clustering number has a higher AIC, and returns the cluster identity vector of each row of the input data.frame using that optimal cluster number.
#' @export
get_best_aic_cluster_assignment <- function(clustering_df = data.frame())  {
  aic_scores <- sapply(1:25, function(.x)  {
    mclust_obj <- Mclust(clustering_df, G =.x)
    return(2*mclust_obj$df - 2*mclust_obj$loglik)
  })
  n_clusters <- which(aic_scores[-length(aic_scores)] < aic_scores[-1])[1]
  n_clusters <- ifelse(length(n_clusters) > 0, n_clusters, length(aic_scores))
  return(Mclust(clustering_df, G = n_clusters)$classification)
}


#' Compute distance between samples
#'
#' Use cidr to compute inter-sample distances and return named distance matrix
#'
#' @return get_cidr_distance_matrix_from_raw_counts returns a named distance matrix showing the distance between samples. Both rows and columns are named.
#'
#' @export
get_cidr_distance_matrix_from_raw_counts <- function(cell_by_ensembl_mat = matrix()) {
  cidr_obj <- cidr::scDataConstructor(apply(cell_by_ensembl_mat, 1, function(.x)  {
    log2((.x/sum(.x)) * 1e5 + 1)
  }))
  cidr_obj <- cidr::determineDropoutCandidates(cidr_obj)
  cidr_obj <- cidr::wThreshold(cidr_obj)
  cidr_obj <- cidr::scDissim(cidr_obj)
  dist_mat <- cidr_obj@dissim
  dimnames(dist_mat) <- list(rownames(cell_by_ensembl_mat), rownames(cell_by_ensembl_mat))
  return(dist_mat)
}


#' Determine genes enriched in regions of a minimum spanning tree
#'
#' Calculate an enrichment score for each gene that indicates whether it is expressed in a particular non-random region of the mininimum-spanning tree
#'
#' This function computes an enrichment score for each gene in the mst using the followeing approach.
#' For each gene-mst pair, the sum of summary statistic differences between all pairs of nodes is determined, and is then divided by the range.
#' This value will be equal to or higher than one, with lower values indicating higher enrichment in the mst.
#' The basis for this is that the more randomly a gene is expressed, the higher the sum of the differences between neihboring nodes will be.
#' This will result in higher values. Values of one will occur when the full range of values occurs along a single edge in the graph,
#' bisecting the graph into high and low expressed regions for this particular gene.
#' Alternatively, a value of one can occur if a gene is gradually increased (or decreased) in one direction across its range, but never changing sign.
#'
#' The above value is then divided by the variance of the measurements in the graph.
#' This increases the enrichment score for genes that are enriched in large subgraphs relative to terminal nodes.
#'
#' @return get_df_of_cluster_enrichments_for_each_gene returnn a modified version  of cluster_level_df, which contains the enrichment scores for each gene from each mst,
#' or NA if the gene was not present in a particular mst
#'
#' @export
get_df_of_cluster_enrichments_for_each_gene <- function(cluster_level_df = data.frame(), mst_igraph_list = list()) {
  sample_level_cluster_enrichments_list <- lapply(names(mst_igraph_list), function(.prefix)  {
    mst_and_node_values_list <- mst_igraph_list[[.prefix]]
    current_mst_graph <- mst_and_node_values_list[[1]]
    enrichment_values <- apply(mst_and_node_values_list[[2]], 2, function(.col)  {
      min_max_dif_value <- abs(max(.col) - min(.col))
      edge_nodes_mat <- get.edges(current_mst_graph, E(current_mst_graph))
      #print(abs(.col[edge_nodes_mat[,1]] - .col[edge_nodes_mat[,2]]))
      total_edge_lengths <- sum(abs(.col[edge_nodes_mat[,1]] - .col[edge_nodes_mat[,2]]))
      log2(total_edge_lengths/(min_max_dif_value * var(.col)))
    })
    return(enrichment_values)
  })
  for(i in 1:length(mst_igraph_list))  {
    trans_df <- rnaseqUtils::get_df_from_named_char_vector(sample_level_cluster_enrichments_list[[i]], c('GENEID', paste0(names(mst_igraph_list)[i], "_clusenrich")))
    cluster_level_df <- dplyr::left_join(cluster_level_df, trans_df)
  }
  return(cluster_level_df)
}
