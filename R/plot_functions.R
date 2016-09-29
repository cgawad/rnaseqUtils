plot_igraph_for_gene <- function(expr_object, mst_graph, kmeans_obj, graph_layout_matrix = NULL, ensembl_id = character(), color_scale_divisor = NULL, color_vec = colorRampPalette(c('red', 'white', 'blue'))(100), color_function = mean, size_function = length)  {
  print(ensembl_id)
  list_of_expr_values_by_cluster_index <- lapply(1:length(unique(kmeans_obj$cluster)),  function(.cluster_num) {
    cluster_level_expr_values <- expr_object[which(rownames(expr_object) %in% names(kmeans_obj$cluster[which(kmeans_obj$cluster == .cluster_num)])), ensembl_id]
    return(cluster_level_expr_values)
  })

  color_variable <- sapply(list_of_expr_values_by_cluster_index, color_function)
  print(color_variable)
  color_values <-rep('grey', length(color_variable))
  if(is.null(color_scale_divisor))  {
    color_scale_divisor <- max(color_variable)
  }
  scaled_vals <- ceiling(color_variable/color_scale_divisor *  100)
  print(scaled_vals)
  color_values[which(color_variable != 0 & !is.na(color_variable))] <- color_vec[scaled_vals[which(scaled_vals != 0 & !is.na(scaled_vals))]]

  print(color_values)

  size_variable <- sapply(list_of_expr_values_by_cluster_index, size_function)

  if(is.null(graph_layout_matrix))  {
    graph_layout_matrix <<- layout.auto(mst_graph)
  }
  plot(mst_graph, vertex.size = size_variable/10, vertex.color = color_values, layout = graph_layout_matrix, vertex.label=NA)
  debug_df<- data.frame(vid = as.character(V(mst_graph)), measure_var = color_variable, color = color_values, size = size_variable, stringsAsFactors = FALSE)
  debug_df$cluster_members <- sapply(1:length(unique(kmeans_obj$cluster)), function(.x)  {paste(names(kmeans_obj$cluster[kmeans_obj$cluster == .x]), collapse=",")})
  return(debug_df)
}

split_samples_and_generate_plot <- function(ensembl_id = character(), sample_list = list(), kmeans_obj_list = list(), mst_object_list = list(), layout_list = list(), sample_names = character(), color_function = mean)  {

  gene_symbol <- paste(select(org.Hs.eg.db, keys = ensembl_id, keytype = 'ENSEMBL', columns = 'SYMBOL')$SYMBOL, collapse = ",")
  max_measure_value <- NULL
  min_non_zero_measure_var <- 1e20

  par(mfrow = c(1, length(sample_list)), oma = c(0,0,2,0))
  for(i in 1:length(sample_list))  {
    debug_fl_igraph_df <- plot_igraph_for_gene(sample_list[[i]], kmeans_obj = kmeans_obj_list[[i]], mst_graph = mst_object_list[[i]], graph_layout_matrix = layout_list[[i]], ensembl_id = ensembl_id)
    title(sample_names[[i]])
    max_measure_value <- max(max_measure_value, max(debug_fl_igraph_df$measure_var))
    min_non_zero_measure_var <- min(min_non_zero_measure_var, min(debug_fl_igraph_df$measure_var[(debug_fl_igraph_df$measure_var != 0) & !is.na(debug_fl_igraph_df$measure_var)]))
  }
  for(i in 1:length(sample_list))  {
    debug_fl_igraph_df <- plot_igraph_for_gene(sample_list[[i]], kmeans_obj = kmeans_obj_list[[i]], mst_graph = mst_object_list[[i]], graph_layout_matrix = layout_list[[i]], ensembl_id = ensembl_id, color_scale_divisor = max_measure_value)
    title(sample_names[[i]])
  }
  mtext(paste0("Mean log expression of ", gene_symbol, " across clusters") , outer = TRUE, cex = 1.5)
}

