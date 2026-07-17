graph_cluster_plot_EBseq <- function (g, expr_list, meta_list, p_time = NULL) {
  require(igraph)
  # comm_obj: community from igraph 
  # g: graph object from igraph 
  # comm_n: number of communities (set to NULL to not cut)
  # co_id: either the results of get_co_clusters or a vector with co_clusters as value and genes as names
  # expr_list: list of expression objects  
  
  v_genes <- V(g)$genes
  v_names <-  V(g)$name
  names(v_genes) <- v_names
  
  flat_v_names <- unlist(lapply(1:length(v_genes), 
                                function(n) rep(names(v_genes[n]), length(v_genes[[n]]))))
  
  flat_v_names <- unlist(flat_v_names)

  gene2cluster <- data.frame(gene = unlist(v_genes), cluster_id =  flat_v_names)

  #print(head(gene2cluster))
  #print(table(gene2cluster$cluster_id))
  p <- plot_cluster_line_EBseq(expr_list, 
                         meta_list, 
                         gene2cluster,
                         p_time) 
  return(p)
} # end community cut   