# 'Load 10X matrix files
#'
#' Loads the matrix files from a 10X run
#'
#' @return  get_cellranger_matrices returns a list of count matrices corresponding to the matrix files of \code{filenames}. The names of the list are the basenames of \code{filenames}

#' @export
get_cellranger_matrices <- function(filenames = character(), sample_name_prefix = NULL, genome_build = 'hg19')  {
  matrix_list <- lapply(filenames, function(.filename)  {
    temp_gbm <- cellrangerRkit::load_cellranger_matrix(.filename, genome = genome_build)
    temp_gbm_mat <- as.matrix(temp_gbm@mat)
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

#' Filter MT and Ribsomal protein genes
#'
#' Removes all mitochondrial and ribosomal protein genes in the count matrix and returns the reduces matrix.
#'
#' @return remove_ribosomal_proteins_and_mitochondrial_genes_from_matrix returns a reduced matrix with all mitochondrial and ribosomal protein genes removed

#' @export
remove_ribosomal_proteins_and_mitochondrial_genes_from_matrix <- function(cell_by_ensembl_mat = matrix(), species = 'hsapiens')  {
  library(biomaRt)
  if(!all(grepl("^ENS", colnames(cell_by_ensembl_mat))))  {
    stop("Column names should be ensembl ids and should start with 'ENS'")
  }
  if(species == 'hsapiens')  {
    ensembl <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
    useDataset("hsapiens_gene_ensembl", mart=ensembl)
    ensembl_gene_metadata_df <- biomaRt::select(ensembl, keys = colnames(cell_by_ensembl_mat), keytype = 'ensembl_gene_id', columns=c('ensembl_gene_id','chromosome_name', "hgnc_symbol", "description")) %>% dplyr::rename(ENSEMBL = ensembl_gene_id)
    ensembl_gene_metadata_MT_and_ribosomal_df <- dplyr::filter(ensembl_gene_metadata_df, grepl("^RP[SL][01-9]", ensembl_gene_metadata_df$hgnc_symbol) | (chromosome_name == 'MT'))
  }
  else {
    warning("This is not a valid option")
    stop()
  }
  ensembl_overlaps <- colnames(cell_by_ensembl_mat)[which(colnames(cell_by_ensembl_mat) %in% ensembl_gene_metadata_MT_and_ribosomal_df$ensembl_gene_id)]
  message(paste(length(ensembl_overlaps), 'of the', ncol(cell_by_ensembl_mat), 'ensembl ids are being removed'))
  return(cell_by_ensembl_mat[,!colnames(cell_by_ensembl_mat) %in% ensembl_overlaps])
}

#' @export
get_random_genes_from_matrix <- function(count_matrix = matrix(), distance_matrix = NULL, normalize_matrix = TRUE, range_of_counts = 2:50, n_sample_draws = 1000, tail_area = 0.2)  {
  #remove zero counts
  sample_names <- rownames(count_matrix)
  lt2_names <- colnames(count_matrix)[apply(count_matrix, 2, function(.col)  {sum(.col) <= 1})]
  gt1_mat <- count_matrix[,!colnames(count_matrix) %in% lt2_names]
  if(normalize_matrix)  {
    gt1_mat <- t(apply(gt1_mat, 1, function(.x) {log2(.x/sum(.x) * 1e6 + 1)}))
  }
  if(is.null(distance_matrix))  {
    cor_mat <- cor(t(gt1_mat))
    dist_mat <- 1 - cor_mat
  }
  else {
    dist_mat <- distance_matrix
  }

  cell_non_zero_gene_counts <- apply(gt1_mat, 2, function(.x)  {sum(.x >0)})
  valid_cell_non_zero_gene_counts <- cell_non_zero_gene_counts[(cell_non_zero_gene_counts %in% range_of_counts)]

  tested_genes_pval_and_dif_mat_list <- lapply(unique(valid_cell_non_zero_gene_counts), function(.count)  {

    mean_pw_cor <- sapply(1:n_sample_draws, function(.ind)  {
      sub_sample_names <- sample(sample_names, .count, replace = FALSE)
      sampled_names_mat <- dist_mat[sub_sample_names,sub_sample_names]
      get_mean_dist_from_matrix(sampled_names_mat)
    })
    mean_mean_pw_cor <- mean(mean_pw_cor)
    #print(mean_pw_cor)
    sd_mean_pw_cor <- sd(mean_pw_cor)

    genes_pval_and_dif_mat <- sapply(names(valid_cell_non_zero_gene_counts[valid_cell_non_zero_gene_counts == .count]), function(.gene_name)  {
      selected_sample_names <- rownames(gt1_mat)[which(gt1_mat[,.gene_name] > 0)]
      selected_mean_val<- get_mean_dist_from_matrix(dist_mat[selected_sample_names,selected_sample_names])
      return(c(mean(mean_pw_cor < selected_mean_val), log10(selected_mean_val/min(mean_pw_cor))))
    })
    #print(head(genes_pval_and_dif_mat))
    return(genes_pval_and_dif_mat)
  })
  return(t(do.call(cbind, tested_genes_pval_and_dif_mat_list)))
  #return(unique(c(unlist(names(genes_to_remove_list)), lt2_names)))
}

#` Calculate mean value from matrix
#'
#' The function calculates and returns the mean of the upper-triangular portion of a matrix
#'
#' @return get_mean_dist_from_matrix returns the upper-triangular mean of temp_mat, exluding the zeros from the diagonal or the lower-triangular portion of the matrix
#'
get_mean_dist_from_matrix <- function(temp_mat)  {
  temp_mat[!upper.tri(temp_mat)] <- 0
  #print(temp_mat)
  mean(temp_mat[temp_mat != 0])
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
