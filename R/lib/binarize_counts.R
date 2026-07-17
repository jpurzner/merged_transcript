binarize_counts <- function (counts, log_data = FALSE, exclude_extreme = FALSE, return_cutoff = FALSE, return_mean = FALSE, return_cutoff_mean = FALSE,  mu= c(0.3, 0.65), label_name = NULL) {
  # take counts and seperate signal from background using mixtools with a mixture of normal distributions
  # bulk of code found at http://stats.stackexchange.com/questions/57993/how-to-explain-how-i-divided-a-bimodal-distribution-based-on-kernel-density-esti 
  # modified to use next gen sequencing count data
  
  
  
	require(mixtools)

	# log transform the count data 
	if (log_data) {
    x <- log10(counts + 1)
	} else { 
	  x <- counts
	}
    
  if (exclude_extreme) {
    x <- counts[counts > 0.1 & counts < 1]
  }
    
    
  
  model <- normalmixEM(x= x, k=2, mu = mu, lambda =c(.6, .4),sigma = c(0.3, 0.2), arbvar = FALSE ,  maxit =10000, maxrestarts = 200, epsilon = 1e-20)
	#model <- normalmixEM(x= x, k=2, mu = NULL, lambda =NULL,sigma = NULL, maxit =10000, maxrestarts = 10, epsilon = 1e-16)
	index.lower <- which.min(model$mu)  # Index of component with lower mean

	find.cutoff <- function(proba=0.5, i=index.lower) {
    	## Cutoff such that Pr[drawn from bad component] == proba
    	f <- function(x) {
        	proba - (model$lambda[i]*dnorm(x, model$mu[i], model$sigma[i]) /
                     (model$lambda[1]*dnorm(x, model$mu[1], model$sigma[1]) + model$lambda[2]*dnorm(x, model$mu[2], model$sigma[2])))
        	}
        # find minimum of f and set lower limit of uniroot
        low <- x[which.min(f(x))]
        up <- max(x)
        
        cutoff <- try(uniroot(f=f, lower=low, upper=up)$root)
        if (class(cutoff) == "try-error") {
          print("warning: failed to find root")
          cutoff <- (model$mu[1] + model$mu[2])/ 2
        }
        
        return(cutoff)  # Careful with division by zero if changing lower and upper
	} # end of find cutoffs
 
	cutoffs <- c(find.cutoff(proba=0.5, i = index.lower), find.cutoff(proba=0.75, i = index.lower))  # Around c(1.8, 1.5)
	
# plot using ggplot2 
x <- as.data.frame(x)
	
norm_funct <- function(x, lambda, mu, sigma) {
	return(lambda*dnorm(x, mu, sigma))
}
x$nd1 <- norm_funct(x$x, model$lambda[1], model$mu[1], model$sigma[1])
x$nd2 <- norm_funct(x$x, model$lambda[2], model$mu[2], model$sigma[2])


  m <- ggplot(x) + 
    geom_histogram(aes(x = x, y = ..density..), binwidth = 0.01, color = "black", fill = "grey", alpha = 0.1) +
    geom_line(aes(x = x, y = nd1), color = "dodgerblue1", size = 1, alpha = 0.8) +
    geom_line(aes(x = x, y = nd2), color = "seagreen4", size = 1, alpha = 0.8) +
    geom_vline(xintercept = cutoffs[1], linetype = "dotdash", alpha = 0.5) +
    theme(axis.title.y = element_text(size = 16, angle = 90),
          axis.text.y = element_text(size = 16),
          axis.text.x = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    xlab("log10 + 1 counts") + 
    ylab("density")
  
  # Add title if label_name is not NULL
  if (!is.null(label_name)) {
    m <- m + ggtitle(label_name)
  }
  
print(m)

# determine cutoff for un-transformed counts
if (log_data) {
  cutoffs <- (10^cutoffs) - 1
}


# binarize the data based upon this cutoff
results <- as.data.frame(counts)
results$bin <- 0 
results[results$counts >= cutoffs[2], 2] <- 1

if (return_cutoff)  {
  return(cutoffs)
} else if (return_mean == TRUE) {
  
  index.upper <- which.max(model$mu) 
  return(model$mu[index.upper])
  
} else if (return_cutoff_mean == TRUE) {  
  
  index.upper <- which.max(model$mu) 
  index.lower <- which.min(model$mu) 
  
  # Log-likelihood
  loglik <- model$loglik
  
  # Calculate AIC/BIC manually
  n <- length(data)
  k <- length(model$lambda) + 2*length(model$mu) # Number of parameters
  aic <- -2*model$loglik + 2*k
  bic <- -2*model$loglik + log(n)*k
  
  all_vals <- c(cutoffs, model$mu[index.upper], model$mu[index.lower], loglik, aic, bic) 
  names(all_vals) <- c("cutoff_50", "cutoff_75", "mean_upper", "mean_lower", "loglik", "aic", "bic")
  return(all_vals)
  
} else {
  return(results)  
}
  

}# end of binarize_counts