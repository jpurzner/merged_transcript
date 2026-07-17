wins_norm_histones <- function(histone_rlog, median_subtract = FALSE, low = 500, high = 20) {

  

  
library(reshape2)
library(pheatmap)
  
  
# from http://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode



estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

col_metrics <- function(x) {
  range_df <- data.frame(mins = apply(x, 2, min), 
                         maxes = apply(x, 2, max),
                         median = apply(x, 2, median),
                         mode = apply(x, 2, estimate_mode))
  return(range_df)
}                        

# winsorize the data 
library(DescTools)

winsorize_by_points <- function(df, low_wins, high_wins) {
  # Handle edge case: if low_wins <= 0, don't clip lower end
  if (low_wins > 0) {
    low_cut <- apply(df, 2, function(x) sort(x)[low_wins])
  } else {
    low_cut <- apply(df, 2, min, na.rm = TRUE)
  }

  # Handle edge case: if high_wins <= 0, don't clip upper end
  if (high_wins > 0) {
    high_cut <- apply(df, 2, function(x) sort(x, decreasing = TRUE)[high_wins])
  } else {
    high_cut <- apply(df, 2, max, na.rm = TRUE)
  }

  # Use pmax/pmin for manual clipping (avoids DescTools API changes)
  wins_df <- sapply(1:ncol(df), function(i) {
    x <- df[,i]
    pmin(pmax(x, low_cut[i]), high_cut[i])
  })
  colnames(wins_df) <- colnames(df)
  row.names(wins_df) <- row.names(df)

  return(wins_df)

}

if (median_subtract) {
  # median subtract the data to better align the transformed values 
  histone_rlog_medsubt <- sweep(histone_rlog, 2, apply(histone_rlog, 2, median), `-`)
} else {
  histone_rlog_medsubt <- histone_rlog
}

histone_rlog_medsubt_wins <- winsorize_by_points(histone_rlog, low, high)
normalize <- function(x) {(x-min(x))/(max(x)-min(x))}
histone_rlog_medsubt_wins_norm <- apply(histone_rlog_medsubt_wins, 2, normalize)

return(as.data.frame(histone_rlog_medsubt_wins_norm))

} # end of function 
