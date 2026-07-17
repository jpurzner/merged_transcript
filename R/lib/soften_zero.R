soften_zero <- function(counts, scale_max = 1.5, breaks = 100) {
  
  histogram <- hist(counts, plot = TRUE, breaks = breaks)
  max_value <- max(histogram$counts[2:length(histogram$counts)])
  zero_indices <- which(counts == 0)
  num_zeros <- sum((counts == 0))
  indices_to_remove <- sample(zero_indices,num_zeros - (max_value * scale_max))
  counts[indices_to_remove] <- NA
  return(counts)
  
}