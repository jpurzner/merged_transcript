time_expr_align <- function(expr_align_list, return_shift = FALSE) {
  
  # takes the interpolated dataset and aligns the dataset 
  # return_shift TRUE will provide a data frame of values to subtract both datasets by 
  
  # only works for 2 datasets search 2DS FIX to find problems 
  
  gene_num <-  dim(expr_align_list[[1]])[1]
  # checks the first gene and determines the overlap 
  ov <- complete.cases(sapply(c(1:length(expr_align_list)), function(n) expr_align_list[[n]][1,]))
  
  # keep only the overlaping segments 
  ds_ov <- lapply(expr_align_list, function(x) x[,ov])
  
  # 2DS FIX 
  shift_by_gene <- apply(ds_ov[[1]] - ds_ov[[2]], 1, median)
  

  
  # 2DS FIX 
  expr_align_list[[2]] <- sweep(expr_align_list[[2]], 1, shift_by_gene, "+")
  
  # now that time series aligned find the common minimum and maximum 
  
  all_min <- apply((sapply(expr_align_list, function(i) apply(i, 1, function(j) min(j, na.rm = TRUE)))), 1, min)
  all_max <- apply((sapply(expr_align_list, function(i) apply(i, 1, function(j) max(j, na.rm = TRUE)))), 1, max)
  
  if (return_shift) {
    shift <- data.frame(row.names = names(shift_by_gene), x = rep(0, length(shift_by_gene)), y = shift_by_gene)
    zero_shift <- ((all_max - all_min) / 2) + all_min 
    shift <- shift - zero_shift
    return(shift)
  }
    
  # shift the data centered upon zero 
  expr_align_list <-  lapply(expr_align_list, function(i) sweep(i, 1,((all_max - all_min) / 2) + all_min , "-") )
  
  

  # sorry for this one liner, here is the breakdown 
  # get all ds for a gene in a data frame  
  # sapply(1:length(expr_align_list), function(k) expr_align_list[[k]])
  
  # get the row mean from the above sapply 
  # apply(    ,1, function (i) mean(i, na.rm = TRUE )
  
  # go through each gene 
  # sapply(c(1:gene_num), function(j)     )
  # again sorry!
  
  return(expr_align_list)
}