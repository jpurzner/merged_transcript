recursive_edgeeater <- function(g, end_num, memb_cutoff = 0) {
  # iterates through the 
  
  v_num <- length(V(g))
  memb_byV <- numeric(length = v_num-end_num)
  gene_num_byV <- numeric(length = v_num-end_num)
  count = 1 
  while (v_num >= end_num) { 
    g_rank <- graph.strength(g, mode ="in") + V(g)$weight
    vdel <- which.min(g_rank)
    print(vdel)
    g <- edgeeat(g, vdel, memb_cutoff = 0) 
    memb_byV[count] <- sum(V(g)$weight)
    gene_num_byV[count] <- length(unlist(V(g)$genes))
    # test internal consistency of membership values
    #if (test_graph_membership(g)) {
    #  warning('vertex weight not consistent with gene membership')  
    #}
    v_num <- length(V(g))
    count = count + 1; 
  }
  metric_df <- data.frame(memb_byV = memb_byV, gene_num_byV = gene_num_byV)
  return(list(graph = g, metrif_df = metric_df))
}