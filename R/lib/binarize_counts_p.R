binarize_counts_p <- function(counts, log_data = FALSE, exclude_extreme = FALSE, return_all = FALSE,
                              # New parameters for mixture model control
                              alpha = NULL, beta = NULL, lambda = NULL,
                              maxit = 10000, maxrestarts = 20, epsilon = 1e-8,
                              auto_init = TRUE, verbose = FALSE) {
  # Fit a 2-component gamma mixture model to separate signal from background
  # Modified to auto-estimate starting parameters and be more robust

  require(mixtools)
  require(ggplot2)

  # log transform the count data if requested
  if (log_data) {
    x <- log10(counts + 1)
  } else {
    x <- counts
  }

  if (exclude_extreme) {
    x <- x[x > 0 & x < 1]
  }

  # Remove NAs and negative values (gamma requires positive)
  x_clean <- x[!is.na(x) & x >= 0]

  # Handle zeros: gamma distribution requires positive values
  # Add small offset to zeros
  min_nonzero <- min(x_clean[x_clean > 0], na.rm = TRUE)
  x_for_fit <- ifelse(x_clean == 0, min_nonzero / 10, x_clean)

  # Auto-estimate starting parameters if not provided
  if (auto_init || is.null(alpha) || is.null(beta) || is.null(lambda)) {

    # Estimate from data quantiles
    q25 <- quantile(x_for_fit, 0.25, na.rm = TRUE)
    q50 <- quantile(x_for_fit, 0.50, na.rm = TRUE)
    q75 <- quantile(x_for_fit, 0.75, na.rm = TRUE)
    data_mean <- mean(x_for_fit, na.rm = TRUE)
    data_var <- var(x_for_fit, na.rm = TRUE)

    # Proportion of low values (potential "background" component)
    prop_low <- mean(x_for_fit < q25, na.rm = TRUE)

    # For gamma: mean = alpha*beta, var = alpha*beta^2
    # So: beta = var/mean, alpha = mean/beta = mean^2/var

    # Component 1: lower/background (concentrated near zero or low values)
    mean1 <- q25
    var1 <- (q25)^2 / 2  # assume tighter distribution
    if (var1 == 0) var1 <- 0.01
    beta1 <- max(var1 / mean1, 0.001)
    alpha1 <- max(mean1 / beta1, 0.1)

    # Component 2: higher/signal
    mean2 <- q75
    var2 <- data_var
    if (var2 == 0) var2 <- 0.01
    beta2 <- max(var2 / mean2, 0.001)
    alpha2 <- max(mean2 / beta2, 0.1)

    # Set defaults if not provided
    if (is.null(alpha)) alpha <- c(alpha1, alpha2)
    if (is.null(beta)) beta <- c(beta1, beta2)
    if (is.null(lambda)) lambda <- c(0.5, 0.5)

    if (verbose) {
      cat("Auto-estimated starting parameters:\n")
      cat("  alpha:", alpha, "\n")
      cat("  beta:", beta, "\n")
      cat("  lambda:", lambda, "\n")
    }
  }

  # Try fitting with multiple random restarts
  best_model <- NULL
  best_loglik <- -Inf

  for (restart in 1:maxrestarts) {
    # Add some noise to starting parameters for restarts > 1
    if (restart > 1) {
      noise <- runif(2, 0.5, 1.5)
      alpha_try <- alpha * noise
      beta_try <- beta * noise
      lambda_try <- c(runif(1, 0.3, 0.7), NA)
      lambda_try[2] <- 1 - lambda_try[1]
    } else {
      alpha_try <- alpha
      beta_try <- beta
      lambda_try <- lambda
    }

    model <- tryCatch({
      suppressWarnings(
        gammamixEM(x = x_for_fit, k = 2,
                   alpha = alpha_try,
                   beta = beta_try,
                   lambda = lambda_try,
                   maxit = maxit,
                   maxrestarts = 3,
                   epsilon = epsilon,
                   verb = FALSE)
      )
    }, error = function(e) {
      if (verbose) cat("Restart", restart, "failed:", e$message, "\n")
      return(NULL)
    })

    if (!is.null(model) && !is.null(model$loglik)) {
      if (model$loglik > best_loglik) {
        best_loglik <- model$loglik
        best_model <- model
        if (verbose) cat("Restart", restart, "- loglik:", model$loglik, "\n")
      }
    }
  }

  if (is.null(best_model)) {
    stop("Failed to fit gamma mixture model after ", maxrestarts, " attempts. Try adjusting starting parameters.")
  }

  model <- best_model

  # Identify which component has lower mean
  mean1 <- model$gamma.pars[1,1] * model$gamma.pars[2,1]
  mean2 <- model$gamma.pars[1,2] * model$gamma.pars[2,2]
  index.lower <- which.min(c(mean1, mean2))

  if (verbose) {
    cat("\nFitted model:\n")
    cat("  Component 1 - alpha:", model$gamma.pars[1,1], "beta:", model$gamma.pars[2,1],
        "lambda:", model$lambda[1], "mean:", mean1, "\n")
    cat("  Component 2 - alpha:", model$gamma.pars[1,2], "beta:", model$gamma.pars[2,2],
        "lambda:", model$lambda[2], "mean:", mean2, "\n")
    cat("  Lower component:", index.lower, "\n")
  }

  # Find cutoff function
  find.cutoff <- function(proba = 0.5, i = index.lower) {
    f <- function(x) {
      proba - (model$lambda[i] * dgamma(x, model$gamma.pars[1,i], 1/model$gamma.pars[2,i]) /
                 (model$lambda[1] * dgamma(x, model$gamma.pars[1,1], 1/model$gamma.pars[2,1]) +
                    model$lambda[2] * dgamma(x, model$gamma.pars[1,2], 1/model$gamma.pars[2,2])))
    }

    # Find where f changes sign
    x_grid <- seq(min(x_for_fit), max(x_for_fit), length.out = 1000)
    f_vals <- sapply(x_grid, f)

    # Find crossing point
    sign_changes <- which(diff(sign(f_vals)) != 0)
    if (length(sign_changes) > 0) {
      low <- x_grid[sign_changes[1]]
      up <- x_grid[min(sign_changes[1] + 1, length(x_grid))]
      cutoff <- tryCatch(
        uniroot(f = f, lower = low, upper = up)$root,
        error = function(e) (mean1 + mean2) / 2
      )
    } else {
      cutoff <- (mean1 + mean2) / 2
    }

    return(cutoff)
  }

  cutoffs <- c(find.cutoff(proba = 0.5, i = index.lower),
               find.cutoff(proba = 0.75, i = index.lower))

  if (verbose) cat("Cutoffs (p=0.5, p=0.75):", cutoffs, "\n")

  # Plot using ggplot2
  plot_df <- data.frame(x = x_for_fit)

  gamma_funct <- function(x, lambda, alpha, beta) {
    return(lambda * dgamma(x, shape = alpha, rate = 1/beta))
  }

  plot_df$nd1 <- gamma_funct(plot_df$x, model$lambda[1], model$gamma.pars[1,1], model$gamma.pars[2,1])
  plot_df$nd2 <- gamma_funct(plot_df$x, model$lambda[2], model$gamma.pars[1,2], model$gamma.pars[2,2])

  # Determine appropriate binwidth
  data_range <- diff(range(plot_df$x, na.rm = TRUE))
  binwidth <- data_range / 50

  m <- ggplot(plot_df) +
    geom_histogram(aes(x = x, y = after_stat(density)), binwidth = binwidth,
                   color = "black", fill = "grey", alpha = 0.3) +
    geom_line(aes(x = x, y = nd1), color = "dodgerblue1", linewidth = 1, alpha = 0.8) +
    geom_line(aes(x = x, y = nd2), color = "seagreen4", linewidth = 1, alpha = 0.8) +
    geom_vline(xintercept = cutoffs[1], linetype = "dotdash", alpha = 0.5) +
    geom_vline(xintercept = cutoffs[2], linetype = "dashed", alpha = 0.5, color = "red") +
    theme(axis.title.y = element_text(size = 16, angle = 90),
          axis.text.y = element_text(size = 16),
          axis.text.x = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    xlab("Expression value") +
    ylab("Density") +
    ggtitle(paste0("Gamma mixture model (loglik: ", round(model$loglik, 2), ")"))

  print(m)

  # Transform cutoffs back if log transform was used
  if (log_data) {
    cutoffs <- (10^cutoffs) - 1
  }

  # Binarize the data based upon cutoff
  results <- data.frame(counts = counts)
  results$bin <- as.integer(counts >= cutoffs[2])

  if (return_all) {
    return(list(results = results, model = model, cutoffs = cutoffs, plot = m))
  } else {
    return(results)
  }

} # end of binarize_counts_p
