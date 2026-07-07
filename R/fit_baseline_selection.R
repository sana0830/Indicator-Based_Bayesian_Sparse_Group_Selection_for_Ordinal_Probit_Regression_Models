# =============================================================================
# Baseline Bayesian selection methods for ordinal probit regression models
# =============================================================================
# This file implements baseline indicator-based Bayesian selection methods for
# ordinal probit regression models.
#
# Main function:
#   fit_baseline_selection()
#
# If group_list = list(), ordinary variable selection is performed.
# If group_list is provided, group selection is performed.
#
# In the group selection method, once a group is selected, all variables in the
# selected group are included.
# =============================================================================


fit_baseline_selection <- function(
  fitting_data,
  iter = 3000,
  a = Inf,
  b = 10,
  theta = 0.5,
  seed_num = 1,
  group_list = list(),
  use_mcse = TRUE,
  check_every = 1000,
  target_hw = 0.06,
  burn_proposed = 1000,
  extra_keep = 2000,
  print_every = 1000,
  verbose = TRUE
) {
  # --------------------------------------------------------------------------
  # Fit baseline Bayesian variable selection or group selection for ordinal
  # probit regression.
  #
  # Data format:
  #   fitting_data is an n x (p + 1) data frame or matrix.
  #   The first p columns are predictors, and the last column is the ordinal
  #   response coded as 1, ..., K.
  #
  # Model setting:
  #   group_list = list()      : ordinary variable selection.
  #   group_list is provided   : group selection.
  #
  # Output:
  #   A list containing posterior samples of gamma, beta, tau, eta, and
  #   information about burn-in and stopping iteration.
  # --------------------------------------------------------------------------
  
  if (!"msm" %in% loadedNamespaces()) {
    library(msm)
  }
  if (!"MASS" %in% loadedNamespaces()) {
    library(MASS)
  }
  
  if (!is.numeric(fitting_data[[ncol(fitting_data)]])) {
    stop("The last column of fitting_data must be the ordinal response coded as 1, ..., K.")
  }
  
  Y_data <- as.numeric(fitting_data[[ncol(fitting_data)]])
  
  if (!all(sort(unique(Y_data)) == seq_len(max(Y_data)))) {
    stop("The ordinal response must be coded as consecutive integers: 1, ..., K.")
  }
  
  if (length(b) != 1 || b <= 0) {
    stop("b must be a positive scalar.")
  }
  
  if (length(theta) != 1) {
    stop("theta must be a scalar.")
  }
  theta <- clamp_prob(theta)
  
  set.seed(seed_num)
  t0 <- Sys.time()
  
  final_burn     <- NA_integer_
  stop_iter      <- iter
  mcse_iteration <- NA_integer_
  last_iter      <- iter
  
  # Initialize data dimensions and model parameters
  n_size <- nrow(fitting_data)
  p_size <- ncol(fitting_data) - 1
  X_data <- as.matrix(fitting_data[, 1:p_size])
  K      <- max(Y_data)
  
  # Complete group_list when group selection is used.
  if (length(group_list) > 0) {
    
    existing_numbers <- unlist(group_list)
    
    if (any(existing_numbers < 1) || any(existing_numbers > p_size)) {
      stop("group_list contains invalid predictor indices.")
    }
    
    if (anyDuplicated(existing_numbers)) {
      stop("Each predictor can only appear in one group.")
    }
    
    missing_numbers <- setdiff(seq_len(p_size), existing_numbers)
    
    if (length(missing_numbers) > 0) {
      group_list <- c(group_list, as.list(missing_numbers))
    }
    
    g_size <- length(group_list)
    eta <- rep(0L, g_size)
    
  } else {
    
    g_size <- 0L
    eta <- NULL
  }
  
  beta  <- rep(0, p_size)
  gamma <- rep(0, p_size)
  
  if (a == Inf) {
    tau <- sort(runif(K - 1, 0, 1))
  } else {
    tau <- sort(rnorm(K - 1, mean = 0, sd = a))
  }
  
  # Allocate storage for posterior samples
  beta_record  <- array(0, dim = c(p_size, iter))
  gamma_record <- array(0, dim = c(p_size, iter))
  tau_record   <- array(0, dim = c(K - 1, iter))
  eta_record   <- if (g_size > 0) array(0, dim = c(g_size, iter)) else NULL
  
  # Initialize latent responses Y_star
  Y_star_data <- rep(0, n_size)
  
  for (i in seq_len(n_size)) {
    
    if (Y_data[i] == 1) {
      Y_star_data[i] <- -abs(rnorm(1)) + tau[1]
      
    } else if (Y_data[i] == K) {
      Y_star_data[i] <- abs(rnorm(1)) + tau[K - 1]
      
    } else {
      lower <- tau[Y_data[i] - 1]
      upper <- tau[Y_data[i]]
      Y_star_data[i] <- runif(1, lower, upper)
    }
  }
  
  # Current linear predictor
  xbeta <- as.numeric(X_data %*% beta)
  
  # Main MCMC loop
  for (counter in seq_len(iter)) {
    
    # ------------------------------------------------------------------------
    # Update beta and gamma
    # ------------------------------------------------------------------------
    
    if (g_size == 0) {
      
      # Ordinary variable selection
      for (j in seq_len(p_size)) {
        
        xj <- X_data[, j]
        
        # Partial residual excluding variable j
        R_j <- Y_star_data - xbeta + xj * beta[j]
        
        xTx <- drop(crossprod(xj))
        
        sigma2_tilde_j <- b^2 / (1 + b^2 * xTx)
        sigma_tilde_j  <- sqrt(sigma2_tilde_j)
        
        mu_tilde_j <- b^2 * sum(R_j * xj) / (1 + b^2 * xTx)
        
        # Posterior inclusion probability of gamma_j
        logG <- 
          0.5 * log(sigma2_tilde_j) -
          log(b) +
          0.5 * mu_tilde_j^2 / sigma2_tilde_j
        
        p_j <- plogis(log(theta) + logG - log1p(-theta))
        
        old_beta_j <- beta[j]
        
        if (runif(1) < p_j) {
          beta[j]  <- rnorm(1, mean = mu_tilde_j, sd = sigma_tilde_j)
          gamma[j] <- 1L
        } else {
          beta[j]  <- 0
          gamma[j] <- 0L
        }
        
        # Incremental update of the linear predictor
        xbeta <- xbeta + xj * (beta[j] - old_beta_j)
      }
      
    } else {
      
      # Group selection
      for (g in seq_len(g_size)) {
        
        group_g <- group_list[[g]]
        len_g   <- length(group_g)
        
        Xg <- X_data[, group_g, drop = FALSE]
        old_beta_g <- beta[group_g]
        
        # Partial residual excluding group g
        R_g <- Y_star_data - xbeta + as.numeric(Xg %*% old_beta_g)
        
        # Posterior distribution under the active group model
        Sig_inv <- crossprod(Xg) + diag(1 / b^2, nrow = len_g)
        bvec    <- crossprod(Xg, R_g)
        
        Rchol <- chol(Sig_inv)
        tmp   <- forwardsolve(t(Rchol), bvec)
        mu_g  <- backsolve(Rchol, tmp)
        
        logdet_Sig_inv <- 2 * sum(log(diag(Rchol)))
        quad <- drop(crossprod(bvec, mu_g))
        
        # Posterior inclusion probability of eta_g
        logG <- 
          -len_g * log(b) -
          0.5 * logdet_Sig_inv +
          0.5 * quad
        
        p_g <- plogis(log(theta) + logG - log1p(-theta))
        
        if (runif(1) < p_g) {
          
          eta[g] <- 1L
          gamma[group_g] <- 1L
          
          if (len_g == 1L) {
            beta[group_g] <- rnorm(
              1,
              mean = as.numeric(mu_g),
              sd   = sqrt(1 / as.numeric(Sig_inv))
            )
          } else {
            Sig <- chol2inv(Rchol)
            beta[group_g] <- drop(MASS::mvrnorm(
              1,
              mu = as.numeric(mu_g),
              Sigma = Sig
            ))
          }
          
        } else {
          
          eta[g] <- 0L
          gamma[group_g] <- 0L
          beta[group_g]  <- 0
        }
        
        # Incremental update of the linear predictor
        xbeta <- xbeta + as.numeric(Xg %*% (beta[group_g] - old_beta_g))
      }
    }
    
    # ------------------------------------------------------------------------
    # Update cutpoints tau
    # ------------------------------------------------------------------------
    
    if (a == Inf) {
      
      for (k in seq_len(K - 1)) {
        
        lower <- if (any(Y_data == k)) {
          max(Y_star_data[Y_data == k])
        } else {
          -Inf
        }
        
        upper <- if (any(Y_data == k + 1)) {
          min(Y_star_data[Y_data == k + 1])
        } else {
          Inf
        }
        
        if (lower < upper) {
          tau[k] <- runif(1, lower, upper)
        }
      }
      
    } else {
      
      tau_star <- rep(0, K - 1)
      
      for (k in seq_len(K - 1)) {
        
        if (k == 1) {
          lower <- -1e20
          upper <- tau[k + 1]
          
        } else if (k == (K - 1)) {
          lower <- tau_star[k - 1]
          upper <- 1e20
          
        } else {
          lower <- tau_star[k - 1]
          upper <- tau[k + 1]
        }
        
        tau_star[k] <- msm::rtnorm(
          1,
          mean  = tau[k],
          sd    = a,
          lower = lower,
          upper = upper
        )
      }
      
      logr <- -(sum(tau_star^2) - sum(tau^2)) / (2 * a^2)
      
      for (i in seq_len(n_size)) {
        
        mu <- xbeta[i]
        
        if (Y_data[i] == 1) {
          
          p_star <- pnorm(tau_star[1], mean = mu, sd = 1)
          p_now  <- pnorm(tau[1],      mean = mu, sd = 1)
          
        } else if (Y_data[i] == K) {
          
          p_star <- 1 - pnorm(tau_star[K - 1], mean = mu, sd = 1)
          p_now  <- 1 - pnorm(tau[K - 1],      mean = mu, sd = 1)
          
        } else {
          
          lo_star <- tau_star[Y_data[i] - 1]
          up_star <- tau_star[Y_data[i]]
          lo_now  <- tau[Y_data[i] - 1]
          up_now  <- tau[Y_data[i]]
          
          p_star <- pnorm(up_star, mean = mu, sd = 1) -
            pnorm(lo_star, mean = mu, sd = 1)
          
          p_now <- pnorm(up_now, mean = mu, sd = 1) -
            pnorm(lo_now, mean = mu, sd = 1)
        }
        
        p_star <- max(p_star, 1e-300)
        p_now  <- max(p_now,  1e-300)
        
        logr <- logr + log(p_star) - log(p_now)
      }
      
      if (log(runif(1)) < min(0, logr)) {
        tau <- tau_star
      }
    }
    
    # ------------------------------------------------------------------------
    # Update latent responses Y_star
    # ------------------------------------------------------------------------
    
    for (i in seq_len(n_size)) {
      
      mean_Y_star <- xbeta[i]
      
      if (Y_data[i] == 1) {
        
        Y_star_data[i] <- msm::rtnorm(
          1,
          mean  = mean_Y_star,
          sd    = 1,
          lower = -Inf,
          upper = tau[1]
        )
        
      } else if (Y_data[i] == K) {
        
        Y_star_data[i] <- msm::rtnorm(
          1,
          mean  = mean_Y_star,
          sd    = 1,
          lower = tau[K - 1],
          upper = Inf
        )
        
      } else {
        
        lower <- tau[Y_data[i] - 1]
        upper <- tau[Y_data[i]]
        
        Y_star_data[i] <- msm::rtnorm(
          1,
          mean  = mean_Y_star,
          sd    = 1,
          lower = lower,
          upper = upper
        )
      }
    }
    
    # Store current MCMC samples
    beta_record[, counter]  <- beta
    gamma_record[, counter] <- gamma
    tau_record[, counter]   <- tau
    
    if (!is.null(eta_record)) {
      eta_record[, counter] <- eta
    }
    
    # Print progress
    if (counter %% print_every == 0 && verbose) {
      cat(sprintf(
        "iter %d / %d | tau = [%s] | elapsed %.1f mins\n",
        counter, iter,
        paste(sprintf("%.3f", tau), collapse = ", "),
        as.numeric(difftime(Sys.time(), t0, units = "mins"))
      ))
    }
    
    # MCSE stopping rule
    if (use_mcse) {
      
      if (counter >= burn_proposed + 50 && counter %% check_every == 0) {
        
        idx <- (burn_proposed + 1):counter
        
        if (!is.null(eta_record)) {
          
          eta_mat <- eta_record[, idx, drop = FALSE]
          res_eta <- mcse_max_halfwidth(eta_mat, alpha = 0.05)
          max_hw_now <- res_eta$max_hw
          
          if (counter %% print_every == 0 && verbose) {
            cat(sprintf(
              "iter = %d, max MCSE half-width of eta = %.6f\n",
              counter, max_hw_now
            ))
          }
          
        } else {
          
          gamma_mat <- gamma_record[, idx, drop = FALSE]
          res_g <- mcse_max_halfwidth(gamma_mat, alpha = 0.05)
          max_hw_now <- res_g$max_hw
          
          if (counter %% print_every == 0 && verbose) {
            cat(sprintf(
              "iter = %d, max MCSE half-width of gamma = %.6f\n",
              counter, max_hw_now
            ))
          }
        }
        
        if (is.na(final_burn) && max_hw_now < target_hw) {
          
          final_burn     <- counter
          stop_iter      <- min(iter, counter + extra_keep)
          mcse_iteration <- counter
          
          if (verbose) {
            cat(sprintf("MCSE criterion reached at iter = %d\n", mcse_iteration))
            cat(sprintf("maximum MCSE half-width = %.6f\n", max_hw_now))
            cat(sprintf("fixed burn-in = %d; stop_iter = %d\n",
                        final_burn, stop_iter))
          }
        }
      }
      
      if (!is.na(final_burn) && counter >= stop_iter) {
        last_iter <- counter
        break
      }
    }
  }
  
  if (is.na(final_burn)) {
    final_burn <- burn_proposed
    last_iter  <- iter
    
    if (use_mcse && verbose) {
      cat("MCSE criterion was not reached; keeping samples after burn_proposed.\n")
    }
  }
  
  kept_idx <- (final_burn + 1):last_iter
  
  result <- list(
    gamma_record   = gamma_record[, kept_idx, drop = FALSE],
    beta_record    = beta_record[, kept_idx, drop = FALSE],
    tau_record     = tau_record[, kept_idx, drop = FALSE],
    eta_record     = if (!is.null(eta_record)) eta_record[, kept_idx, drop = FALSE] else NULL,
    burn_in        = final_burn,
    last_iter      = last_iter,
    mcse_iteration = mcse_iteration
  )
  return(result)
}
