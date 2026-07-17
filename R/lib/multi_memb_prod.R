multi_memb_prod <- function(all_memb_matched) {
  require(pbapply)
  # product of number of clusters for each dataset 
  co_cluster_number <- sapply(1:length(all_memb_matched), function (d) dim(all_memb_matched[[d]])[2])
  # generate the index matrix 
  co_cluster_index <- expand.grid(lapply(co_cluster_number, function(x) seq(1:x)))
  # create hyphenated names 
  idmaker = function(vec){
    return(paste(vec, collapse="-"))
  }
  co_id_label <- apply(as.matrix(co_cluster_index), 1, idmaker)
  
  memb_prod <- function(i) {
    ival <- sapply(1:(dim(co_cluster_index)[2]), 
                function (ds) all_memb_matched[[ds]][,co_cluster_index[i,ds]])
    return(apply(ival, 1, prod))
  }
  
  memb_prod <- pbsapply(1:dim(co_cluster_index)[1], function(i) memb_prod(i))
  colnames(memb_prod) <- co_id_label
  row.names(memb_prod) <- row.names(all_memb_matched[[1]])
  return(memb_prod)
  
} # end of multi_overlap