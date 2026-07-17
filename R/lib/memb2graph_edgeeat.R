memb2graph_edgeeat <- function(memb_prod, memb_prod_rank, steps = 5) {
  require(plyr)
  require(igraph)
  require(ggplot2)
  
  
  # trying to find out why NA values showing up in graph
  #print(memb_prod)
  
  # generate a df of membership by rank for each gene 
  memb_prod_byrank <- lapply(1:(dim(memb_prod_rank)[1]), function(n) memb_prod[n,match(1:steps,memb_prod_rank[n,])])
  memb_prod_byrank <- data.frame(matrix(unlist(memb_prod_byrank), nrow=dim(memb_prod_rank)[1], byrow=T))
  row.names(memb_prod_byrank) <- row.names(memb_prod)
  colnames(memb_prod_byrank) <- seq(1:steps)
  byrank_molten <- melt(memb_prod_byrank)
  plot(ggplot(byrank_molten, aes(x = variable, y = value)) + geom_violin())
 
  
  # generate a df of coclust by rank for each gene
  coclust_index <- lapply(1:(dim(memb_prod_rank)[1]), function(n) colnames(memb_prod_rank)[match(1:steps,memb_prod_rank[n,])])
  coclust_index <- data.frame(matrix(unlist(coclust_index), nrow=dim(memb_prod_rank)[1], byrow=T))
  colnames(coclust_index) <- colnames(memb_prod_byrank)
  row.names(coclust_index) <- row.names(memb_prod_byrank)
  #print(head(coclust_index))
  
  #Trouble shooting in case the row.names give error
  #print(length(row.names(memb_prod_byrank)))
  #print(length(c((1+((2-1)*dim(memb_prod_byrank)[1])):(dim(memb_prod_byrank)[1]*2))))
  
  # generate a flat df with node1 node2 and membership transfer for each gene 
  
  # in this approach we are going represent each edge as a contribution of membership to another node 
  # to do this each gene will contribute the value of the membership for rank 2,3,4,5 
  # these will then be summed as before to determine the total membership contributin of the edge 
  
  edge_df <- lapply(c(2:steps), 
                    function(n) data.frame(
                      gene = row.names(memb_prod_byrank), 
                      from = coclust_index[,n-1],
                      to = coclust_index[,n],
                      membership = memb_prod_byrank[,n],  
                      row.names = c((1+((n-1)*dim(memb_prod_byrank)[1])):(dim(memb_prod_byrank)[1]*n))))
  
  edge_df <- do.call("rbind", edge_df)
  
  # collapse by sum along all unige N1 and N2 
  edge_df_collapse <- ddply(edge_df, c("from", "to"), summarise, 
                            weight = sum(membership),
                            genes = paste(gene, collapse =  ","),
                            gene_membership = paste(membership, collapse =  ","))
  
  # set vertices as the gene number * membership 
  vertex_df <-  data.frame(
                      gene = row.names(memb_prod_byrank), 
                      name = coclust_index[,1],
                      membership = memb_prod_byrank[,1], 
                      row.names = c(1:(dim(memb_prod_byrank)[1])))
  
  
  
  vertex_df_collapse <- ddply(vertex_df, "name", summarise, 
                              weight = sum(membership),
                              genes = paste(gene, collapse =  ","),
                              gene_membership = paste(membership, collapse =  ","))
  #return(vertex_df_collapse)
  
  vertex_genes <- strsplit(vertex_df_collapse[,3], ",")
  vertex_gene_membership <- strsplit(vertex_df_collapse[,4], ",")
  vertex_gene_membership <- lapply(vertex_gene_membership, as.numeric)
  
  vertex_df_collapse$genes <- NULL
  vertex_df_collapse$gene_membership <- NULL
  
  #print(head(vertex_df_collapse))
    
  # prune all edges that are not connected represented in the rank1 only vertexes. 
  edge_df_collapse <- subset(edge_df_collapse, (edge_df_collapse$from %in% vertex_df_collapse$name) & 
                               (edge_df_collapse$to %in% vertex_df_collapse$name))
  
  edge_genes <- strsplit(edge_df_collapse[,4], ",")
  edge_gene_membership <- strsplit(edge_df_collapse[,5], ",")
  edge_gene_membership <- lapply(edge_gene_membership, as.numeric)
  edge_df_collapse$genes <- NULL
  edge_df_collapse$gene_membership <- NULL
  
  g <- graph.data.frame(edge_df_collapse, directed=TRUE, vertices=vertex_df_collapse)
  g <- set.vertex.attribute(g, "genes", value = vertex_genes)
  g <- set.vertex.attribute(g, "gene_membership", value = vertex_gene_membership)
  g <- set.edge.attribute(g, "genes", value = edge_genes)
  g <- set.edge.attribute(g, "gene_membership", value = edge_gene_membership)
  
  
  return(g)
  
}