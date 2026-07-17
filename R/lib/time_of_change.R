time_of_change <- function(df, thres = 0.2) {
  # returns a data.frame with inforamtion about how timecourse data is increasing or decreasing
  # threshold the diff1 here, set anything near 0 to 0 
  
  
  diff1 <- t(diff(t(as.matrix(df)), lag = 1))
  
  d1sum <- rowSums(diff1)
  
  # take index and get the max/min values 
  diff1_i_max <- apply(diff1, 1, which.max)
  diff1_i_min <- apply(diff1, 1, which.min)
  diff1_max <- apply(diff1, 1, max)
  diff1_min <- apply(diff1, 1, min)
  
  
  d1 <- as.matrix(diff1)
  #print(head(d1))
  d1_f <- cut(d1 , c(-Inf, -1*thres, thres, Inf), labels = c(-1, 0 ,1 ))   
  #print(head(d1_f))
  d1_f <- matrix(d1_f, ncol = dim(d1)[2], byrow = FALSE)
  #print(head(d1_f))
  d1_rle <- apply(d1_f, 1, rle)
  
  d1sum <- rowSums(df)
  # want to know what the longest run is 
  
  mod_rle <-  function(g_rle, minimum = FALSE) {
    # crete a mask array with only the values for the longest non-zero run 
    
    # use a temporary variable to remove unwanted 0 or 1, -1 
    tmp_length <- g_rle$lengths
    tmp_length[g_rle$values == "0"] <- 0  
    if (minimum) {
      tmp_length[g_rle$values == "1"] <- 0 
      mi <- which.max(tmp_length)
      g_rle$values[g_rle$values == "1"] <- "garbage"
    } else {
      tmp_length[g_rle$values == "-1"] <- 0
      mi <- which.max(tmp_length)
      g_rle$values[g_rle$values == "-1"] <- "garbage"
    }
    g_rle$values[!(c(1:length(g_rle$values)) %in%  mi)] <- "garbage"
    return(g_rle) 
  }
  
  
  #print(head(d1_rle))

  
  #print(length(d1_rle))
  d1_rle_max  <-  lapply(d1_rle, mod_rle)
  #print(length(d1_rle_max))
  #print(head(d1_rle_max))
  d1_rle_max <- lapply(d1_rle_max, inverse.rle)
  #print(length(d1_rle_max))
  #print(head(d1_rle_max))
  d1_rle_max<- matrix(unlist(d1_rle_max), ncol = dim(diff1)[2], byrow = TRUE)
  #d1_rle_max <- matrix(unlist(lapply(d1_rle_max, inverse.rle)), ncol = 3, byrow = FALSE)
  #print(dim(d1_rle_max))
  d1_rle_max[d1_rle_max %in% c("0", "garbage")] <- 0
  #print(dim(d1_rle_max))
  d1_rle_max[d1_rle_max == "1"] <- 1
  #print(dim(d1_rle_max))
  d1_rle_max <- apply(d1_rle_max, 2, function(x) as.numeric(as.character(x))) 
  
  d1_span_max <- rowSums(diff1 * d1_rle_max )
  
  
  #print(dim(d1_rle_max))
  #print(head(d1_rle_max))
  
  d1_rle_min  <-  lapply(d1_rle, mod_rle, minimum = TRUE)
  #print(head(d1_rle_min))
  d1_rle_min <- matrix(unlist(lapply(d1_rle_min, inverse.rle)), ncol = dim(diff1)[2], byrow = TRUE)
  d1_rle_min[d1_rle_min %in% c("0", "garbage")] <- 0
  d1_rle_min[!(d1_rle_min %in% c("0", "garbage"))] <- 1
  d1_rle_min <- apply(d1_rle_min, 2, function(x) as.numeric(as.character(x))) 
  
  d1_span_min <- rowSums(diff1 * d1_rle_min )
  
  #print(head(d1_rle_min))
  
  max_vals <- apply(df, 1, max)
  min_vals <- apply(df, 1, min)

  
  tof <- data.frame(max.val = max_vals, 
                    min.val = apply(df, 1, min), 
                    max.i = apply(df, 1, which.max), 
                    min.i = apply(df, 1, which.min), 
                    change.total = max_vals - min_vals,
                    inc.start = apply(d1_rle_max, 1, function(x) min(which(x == 1), 100)), 
                    inc.stop = apply(d1_rle_max, 1, function(x) max(which(x == 1), -100)),
                    inc.span = rowSums(d1_rle_max),
                    inc.span.total = d1_span_max,
                    inc.max.step = diff1_max,
                    inc.maxi = diff1_i_max,
                    dec.start = apply(d1_rle_min, 1, function(x) min(which(x == 1), 100)), 
                    dec.stop = apply(d1_rle_min, 1, function(x) max(which(x == 1), -100)), 
                    dec.span = rowSums(d1_rle_min),
                    dec.span.total = d1_span_min,
                    dec.min.step = diff1_min,
                    dec.mini = diff1_i_min
            
  )
  
  tof[abs(tof) == 100] <- 0
  return(tof)
  
}