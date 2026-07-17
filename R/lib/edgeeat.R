edgeeat <- function(g, vdel, memb_cutoff = 0) {
  
  # removes a vertex and adds the membership contribution from the lost 
  # vertex to the remaining vertex, these  
  
  # alarm if any vertex is na
  if (any(is.na(V(g)$weight))) {
    warning("Vertex has <NA> value")
    
  }
  
  #print(vdel)
  # find adjacent vertices
  nVtx <- neighbors(g, vdel, mode = "out") 
  #print(nVtx)
  nVtx <- nVtx[! is.na(nVtx$name)]
  dEdg <- incident(g, vdel, mode = "out")
  #print(dEdg)
  dEdg <- dEdg[! is.na(dEdg$name)]
  
  # handle na values in the vertices, not sure where they are being generated 
  
  # need to handle 2 cases 
  # 1 if the edge carries membership for genes that have a higher rank elsewhere
  # then the membership should not be transfered.  
  start_all_genes <- unlist(V(g)$genes)
  #print(length(all_genes))
  Vdel_genes <- V(g)[[vdel]]$genes
  #print(length(Vdel_genes))
  all_genes <- setdiff(start_all_genes, Vdel_genes)
  #print(length(all_genes))
  
  # troubleshooting: check that there are no gene duplications in the edges 
  # from the vertex to be deleted, there shouldn't be since the input graph is 
  # based upon directional hops
  #----------------------------------------------------------------------------------
  if(any(duplicated(unlist(dEdg$genes)))) stop("genes in edges from deleted vertex not unique")
    
  edge_genes <- get.edge.attribute(g, "genes", index=dEdg) 
  edge_num <- length(edge_genes)
  #handle case where no remaining edges
  if (length(edge_genes) > 0 ) {
    
    next_rank_index <- lapply(c(1:edge_num), function(n) !(edge_genes[[n]] %in% all_genes))
    edge_gene_memb <- get.edge.attribute(g, "gene_membership", index=dEdg) 
    edge_gene_memb_next <- lapply(c(1:edge_num), function(n) edge_gene_memb[[n]][next_rank_index[[n]]])
    edge_gene_name_next <- lapply(c(1:edge_num), function(n) edge_genes[[n]][next_rank_index[[n]]])
    #print(edge_gene_name)
    #print(edge_gene_name_next)
  
    # at this point filter low membership genes 
    if (memb_cutoff > 0) {
      memb_cut_index <- lapply(c(1:edge_num), function(n) (edge_gene_memb[[n]] > memb_cutoff) )
      edge_gene_memb_next <- lapply(c(1:edge_num), function(n) edge_gene_memb[[n]][memb_cut_index[[n]]])
      edge_gene_name_next <- lapply(c(1:edge_num), function(n) edge_genes[[n]][memb_cut_index[[n]]])
    } # end if for memb_cutoff > 0 
    
    # sum the membership scores
    edge_gene_memb_next_sum <- lapply(c(1:length(edge_gene_memb_next)), function(n) sum(edge_gene_memb_next[[n]]))
  
    # add the weight of the edges to the new vertex 
    #print(nVtx$name)
    #print(unlist(edge_gene_memb_next_sum))
    V(g)[nVtx$name]$weight <- (nVtx$weight + unlist(edge_gene_memb_next_sum))
    
    # add genes to the new vertex 
    V(g)[nVtx$name]$genes <- lapply(c(1:length(edge_genes)), function(n) c(nVtx[[n]]$genes, edge_gene_name_next[[n]]))
    
    # ----------------------------------------------------------------------------------------------
    # 2 destroying the vertex will break the hierarchy to other vertices 
    # so an edge needs to be created spanning the other vertex
    # to find these compare the in / out edges for the deleted vertex 
    # ----------------------------------------------------------------------------------------------
    # overall strategy is to flatten all the lists and then use the order to merge the in / out data
    # take the edge genes that are not top rank after deletion 
    
    reassign_index <- lapply(c(1:edge_num), function(n) (edge_genes[[n]] %in% all_genes))
    reassign_out <- lapply(c(1:edge_num), function(n) edge_genes[[n]][reassign_index[[n]]])
    names(reassign_out) <- ends(g,dEdg)[,2]
    # takes the list and makes a flat vector with all the names (there should be a base function for this)
    flat_vtx_out <- unlist(lapply(1:length(reassign_out), 
                           function(n) rep(names(reassign_out[n]), length(reassign_out[[n]]))))
    
    flat_genes_out <- unlist(reassign_out)
    # membership score is the lower out membership 
    reassign_out_memb <- lapply(c(1:edge_num), function(n) edge_gene_memb[[n]][reassign_index[[n]]])
    flat_memb_out <- unlist(reassign_out_memb)
    # create a df so we can merge across gene names  
    reassign_out_df <- data.frame(gene = flat_genes_out, to = flat_vtx_out, membership = flat_memb_out) 
    
    
    #------------------------------------------------------------------------------------------  
    # same as above except with the in edges
    dEdgIn <- incident(g, vdel, mode = "in")
    # stop if no genes for in edges 
    if (length(dEdgIn) > 0) {
      reassign_in <- dEdgIn$genes
      names(reassign_in) <-  ends(g,dEdgIn)[,1]
      flat_vtx_in <- unlist(lapply(1:length(reassign_in), 
                                  function(n) rep(names(reassign_in[n]), length(reassign_in[[n]]))))
      flat_genes_in <- unlist(reassign_in)
      reassign_in_df <- data.frame(gene = flat_genes_in, from = flat_vtx_in) 
      # gene based table of membership scores and to/from data  
      # handle the data as in memb2graph_edgeeat()
      reassign_df <- merge(x = reassign_in_df, y = reassign_out_df, by.x = 1, by.y = 1, all = F)
    
      reassign_df_collapse <- ddply(reassign_df, c("from", "to"), summarise, 
                              weight = sum(membership),
                              genes = paste(gene, collapse =  ","),
                              gene_membership = paste(membership, collapse =  ","))
      
      # stop if no overlap between in / out edges 
      if (dim(reassign_df_collapse)[1] > 0) {    
      
        reassign_edge_genes <- strsplit(as.character(reassign_df_collapse[,4]), ",")
        reassign_edge_gene_membership <- strsplit(as.character(reassign_df_collapse[,5]), ",")
        reassign_edge_gene_membership <- lapply(reassign_edge_gene_membership, as.numeric)
        reassign_df_collapse$genes <- NULL
        reassign_df_collapse$gene_membership <- NULL
    
        # reassign all the attributes to the new edges
        # get current weight of existing edges (0 value means no edge)
        current_weight <- g[from = as.character(reassign_df_collapse$from), 
          to = as.character(reassign_df_collapse$to), 
          attr="weight"]  
    
        # modify the weight, which creates the non existing edges 
        g[from = as.character(reassign_df_collapse$from), 
          to = as.character(reassign_df_collapse$to), 
          attr="weight"] <- current_weight +  reassign_df_collapse$weight 
        # add the genes to the edge 
        current_genes <- g[from = as.character(reassign_df_collapse$from), 
          to = as.character(reassign_df_collapse$to), 
          attr="genes"]
    
        current_reassign_genes <- mapply(c, current_genes, reassign_edge_genes, SIMPLIFY=FALSE)
        current_reassign_genes <- lapply(1:length(current_reassign_genes), 
                                     function (n) current_reassign_genes[[n]][!is.na(current_reassign_genes[[n]])])
    
        g[from = as.character(reassign_df_collapse$from), 
          to = as.character(reassign_df_collapse$to), 
          attr="genes"] <- current_reassign_genes
    
        # add the gene_membership to each edge 
        current_gene_membership <- g[from = as.character(reassign_df_collapse$from), 
                       to = as.character(reassign_df_collapse$to), 
                       attr="gene_membership"]
    
        current_reassign_gene_membership <- mapply(c, current_gene_membership, reassign_edge_gene_membership, SIMPLIFY=FALSE)
        current_reassign_gene_membership <- lapply(1:length(current_reassign_gene_membership), 
                                     function (n) current_reassign_gene_membership[[n]][!is.na(current_reassign_gene_membership[[n]])])
    
        g[from = as.character(reassign_df_collapse$from), 
          to = as.character(reassign_df_collapse$to), 
          attr="gene_membership"] <- current_reassign_gene_membership
      } # end if no overlap between in / out genes 
    } # end if for no genes in in edges  
  } # end if for no genes in out edges 
  
  # delete vertex from graph 
  g <- delete_vertices(g, vdel)
  
  
  # troubleshooting, the total number of genes should not change
  # -------------------------------------------------------------------
  #if (memb_cutoff == 0) {
  #  end_all_genes <- unlist(V(g)$genes)
  #  if (length(end_all_genes) != length(start_all_genes)) {
  #    stop("genes are lost with no memb_cutoff: indicates graph problem")
  #  }
  #}

  
  # sanity check, do our edge membership values make sense
  #--------------------------------------------------------------------
  #next_rank_index_not <- lapply(c(1:length(edge_genes)), function(n) edge_genes[[n]] %in% all_genes)
  #edge_gene_memb_next_not <- lapply(c(1:length(edge_genes)), function(n) edge_gene_memb[[n]][next_rank_index_not[[n]]])
  #edge_gene_memb_next_sum_not <- lapply(c(1:length(edge_gene_memb_next_not)), function(n) sum(edge_gene_memb_next_not[[n]]))
  
  # generate a df containing the edge membership to each connected vertex 
  #edge_df <- data.frame(neighbor_vtx = nVtx$name, 
  #                      edge_memb = dEdg$weight, 
  #                      edge_new_memb = unlist(edge_gene_memb_next_sum), 
  #                      edge_new_memb_not = unlist(edge_gene_memb_next_sum_not))
  #
  
  
  return(g)
  
  
}