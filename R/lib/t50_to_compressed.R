t50_to_compressed <- function (t50, time_ref) {
  require(zoo)
  
  t_ref <- time_ref
  t_ref$i <- 1:nrow(time_ref)
  t_ref$short_i <- cumsum(time_ref[,1] | time_ref[,2])
  t_ref$short_i[!(time_ref[,1] | time_ref[,2])] <- NA
  t_ref$short_i <- na.approx(t_ref$short_i)
  af <- approxfun(t_ref$i, t_ref$short_i)
  return(af(t50))
} 