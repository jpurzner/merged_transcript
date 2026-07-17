expr_half_max_min <- function(expr_align_list) {

# calculates the time of largest increase and largest decrease 
  
  
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
    largest_dir <- ifelse( which.min(j_trim) > which.max(j_trim), -1, 1)
    largest_diff <- max(j_trim) - min(j_trim)
    max_i <- tryCatch(approx(x = j_trim, y = c(1:length(j_trim)), xout = mid_exp)$y, error=function(e) NA)
    max_i <- max_i + start_j

    # the the 2 excluded segments for change
    leftover_i <- list( c(1:(start_j)), 
                        c(end_j:length(j)) 
                      )
    
    leftover_diff <- c( max(j[leftover_i[[1]]]) - min(j[leftover_i[[1]]]),
                        max(j[leftover_i[[2]]]) - min(j[leftover_i[[2]]]) 
                       )
    
    leftover_dir <- c( ifelse( which.min(j[leftover_i[[1]]]) > which.max(j[leftover_i[[1]]]), -1, 1), 
                       ifelse( which.min(j[leftover_i[[2]]]) > which.max(j[leftover_i[[2]]]), -1, 1)
                     )

    # handle case when no different direction is available
    if (all(leftover_dir ==  largest_dir) |  max(c(leftover_diff[leftover_dir != largest_dir], 0), na.rm = TRUE) == 0 ) {
      leftover_max_diff <- 0 
      leftover_early_late <- 0
      leftover_max_dir <- largest_dir * -1
      leftover_j <- 0 
      leftover_mid_exp <- 0  
      leftover_max_i <- 0 
      
      
    } else {
      # set incorrect direction to NA
      leftover_diff[leftover_dir == largest_dir]  <- NA;
      leftover_max_diff <- max(leftover_diff, na.rm = TRUE)
      leftover_early_late <- which.max(leftover_diff)
      leftover_max_dir <- leftover_dir[leftover_early_late]
      
#      print(leftover_dir)
#      print(leftover_diff)
#      print(leftover_i)
#      print(leftover_dir != largest_dir)
#      print(leftover_early_late) 
#      print(leftover_i[[leftover_early_late]])
#      print(j[leftover_i[[leftover_early_late]]])
      
      leftover_j <- j[leftover_i[[leftover_early_late]]]
      leftover_mid_exp <- (max(leftover_j) + min(leftover_j))/2         
      
  
      if (length(leftover_j > 2))  {   
        
        leftover_max_i <- tryCatch(approx(x = leftover_j, 
                                          y = c(1:length(leftover_j)), 
                                          xout = leftover_mid_exp)$y
                                   , error=function(e) 0)
        
        leftover_max_i <- leftover_max_i + min(leftover_i[[leftover_early_late]]) - 1 
        
      } else {
        leftover_max_i <- 0 
      }    
      
            
    }

    

    
    if (largest_dir == 1) {  
          res <-  c(max_i,  largest_diff, largest_dir, leftover_max_i, leftover_max_diff, leftover_max_dir)
    } else {
          res <- c(leftover_max_i, leftover_max_diff, leftover_max_dir, max_i,  largest_diff, largest_dir)
    }
          
    names(res) <- c("up_t50", "up_log2fc","up_dir", "down_t50", "down_log2fc","down_dir" )
#    print(res)
    
    return(res)
    
#    return(max_i)
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
  half_max <- t(apply(expr_align_list, 1,  function(j) find_y0(j)))
} # end if 
return(half_max)
}