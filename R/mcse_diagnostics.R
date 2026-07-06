# =============================================================================
# MCSE diagnostics
# =============================================================================
# This file contains helper functions for MCSE-based convergence monitoring.
# The functions are mainly used to check the stability of posterior inclusion
# probabilities for variable indicators gamma and group indicators eta.
# =============================================================================

mcse_max_halfwidth <- function(gamma_mat, alpha = 0.05) {
  # --------------------------------------------------------------------------
  # Purpose
  # --------------------------------------------------------------------------
  # Compute the maximum Monte Carlo standard error (MCSE) half-width for
  # posterior inclusion probabilities using the consistent batch means method.
  #
  # This function is mainly used to monitor the convergence of indicator chains,
  # such as variable indicators gamma or group indicators eta.
  #
  # Inputs:
  #   gamma_mat : d x S matrix of posterior samples of indicators.
  #               Each row is one indicator chain, and each column is one MCMC
  #               iteration.
  #   alpha     : Significance level for the confidence interval.
  #               Default is 0.05.
  #
  # Outputs:
  #   max_hw    : Maximum MCSE half-width across all indicator chains.
  #   mcse      : MCSE estimate for each indicator chain.
  #   a         : Number of batches.
  #   b         : Batch size.
  # --------------------------------------------------------------------------
  
  S <- ncol(gamma_mat)
  d <- nrow(gamma_mat)
  
  if (S < 50 || d == 0) {
    return(list(
      max_hw = Inf,
      mcse   = rep(Inf, d),
      a      = NA_integer_,
      b      = NA_integer_
    ))
  }
  
  # Batch size and number of batches
  b <- floor(S^(1/3))
  a <- floor(S / b)
  
  if (a <= 1) {
    return(list(
      max_hw = Inf,
      mcse   = rep(Inf, d),
      a      = a,
      b      = b
    ))
  }
  
  # Use the first a * b samples so that all batches have equal length.
  S_use <- a * b
  gamma_use <- gamma_mat[, 1:S_use, drop = FALSE]
  
  # Compute MCSE for each indicator chain.
  mcse_vec <- apply(gamma_use, 1, function(x) {
    
    chain_mean <- mean(x)
    
    batch_means <- colMeans(matrix(
      x,
      nrow = b,
      ncol = a
    ))
    
    batch_var <- sum((batch_means - chain_mean)^2) / (a - 1)
    
    # Estimated asymptotic variance
    sigma2_hat <- b * batch_var
    
    # MCSE of the posterior mean estimate
    sqrt(sigma2_hat / S_use)
  })
  
  # Half-width of the confidence interval
  tval <- qt(1 - alpha / 2, df = a - 1)
  half_width <- tval * mcse_vec
  
  return(list(
    max_hw = max(half_width),
    mcse   = mcse_vec,
    a      = a,
    b      = b
  ))    
}
