expr_half_max <- function(expr_align_list) {

# calculates the half max time point between the 
  
  
find_y0 <- function (j) {
  j <- j[!is.na(j)]
  if (max(j) < 0) { 
    return(NA)
  } else { 
    # discard points that are not between the max and min 
    start_j <- min(which.min(j), which.max(j)) 
    end_j <- max(which.min(j), which.max(j)) 
    j_trim <- j[c(start_j:end_j)]
    mid_exp <- (max(j) + min(j))/2
    max_i <- tryCatch(approx(x = j_trim, y = c(1:length(j_trim)), xout = mid_exp)$y, error=function(e) NA)
    max_i <- max_i + start_j
    return(max_i)
  }
}

if (is.list(expr_align_list)) {
  ds_time <- lapply(expr_align_list, function (x) which(!is.na(x[1,])))  
  all_time <- length(expr_align_list[[1]][1,])
  half_max <- sapply(c(1:length(expr_align_list)), 
                     function(i) apply(expr_align_list[[i]], 1, 
                                       function(j) min(find_y0(j), max(ds_time[[i]]), na.rm = TRUE) ))  
} else { 
  ds_time <- which(!is.na(expr_align_list[1,]))
  all_time <- length(expr_align_list[1,])
  half_max <- apply(expr_align_list, 1,  function(j) min(find_y0(j), max(ds_time), na.rm = TRUE) ) 
} # end if 
return(half_max)
}