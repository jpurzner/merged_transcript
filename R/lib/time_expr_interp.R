time_expr_interp <- function(expr_list, time_ref, meta = NULL ) {
  
  # expr_list is a list of dataframes that contain genes (row) data sets (columns)
  #
  #  time_ref in the following format 
  #    scott hatten
  #   E15  TRUE  FALSE
  #   E16 FALSE  FALSE
  #   E17 FALSE  FALSE 
  #   E18 FALSE  FALSE
  #   E19 FALSE  FALSE
  #   E20 FALSE  FALSE
  #
  # generated like this 
  # time_ref  <- data.frame( time = c(str_c("E", c(15:20)), str_c("P", c(0:56)))) 
  # time_ref$scott <- time_ref$time %in% c("E15", "P1", "P7", "P14")
  # time_ref$hatten <- time_ref$time %in% c("P0","P7","P12","P18","P21","P56")
  #
  # meta: used to indicates replicates that should be pooled to make the time course 
  #
  # uses approx function to interpolate the data
  # 
  
  # get sizes 
  gene_num <-  dim(expr_list[[1]])[1]
  ds_num <- length(expr_list)
  interp_size <- dim(time_ref)[1]
  
  # not needed since using apply 
  # make empty list of data frames 
  #empty_df <- matrix(0L, nrow = gene_num, ncol = interp_size)
  #empty_df <- as.data.frame(empty_df)
  #row.names(empty_df) <- row.names(expr_list[[1]])
  #colnames(empty_df) <- row.names(time_ref)
  #gene_ap <- replicate(ds_num, empty_df, simplify=FALSE)
  
  gene_ap <- list() 
  
  
  for (ds in 1:ds_num) {
    temp <- apply(expr_list[[ds]], 1,  function(x) approx(x = which(time_ref[,ds]), y =x, xout = c(1:interp_size))$y)
    #temp <- sapply(c(1:gene_num),  function(x) approx(x = which(time_ref[,ds]), y = expr_list[[ds]][x,], xout = c(1:interp_size)))
    temp <- t(temp)
    colnames(temp) <- row.names(time_ref)
    gene_ap[[ds]] <- temp  
  }

  return(gene_ap)  
}