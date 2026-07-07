# =============================================================================
# Bayesian sparse group selection for ordinal probit regression models
# =============================================================================
# This file implements the proposed indicator-based Bayesian sparse group
# selection method for ordinal probit regression models.
#
# The main function is:
#   fit_ordinal_sparse_group()
#
# Supporting functions include:
#   compute_group_logZ()       : computes group-level log Bayes factors.
#   cgs_transit()              : updates within-group variable indicators.
#   update_groups_one_sweep()  : performs one group-wise Gibbs update.
# =============================================================================

compute_group_logZ <- function(
  Xg, 
  Res, 
  tau_g, 
  prob, 
  sigma = 1){
  # --------------------------------------------------------------------------
  # Purpose
  # --------------------------------------------------------------------------
  # Compute the log Bayes factor Z_g for one group by enumerating all possible
  # variable subsets within the group.
  #
  # The returned logZ is used to sample the group indicator eta_g.
  #
  # Inputs:
  #   Xg     : n x p_g design matrix for group g.
  #   Res    : Partial residual excluding group g.
  #   tau_g  : Slab standard deviations for variables in group g.
  #   prob   : Prior exclusion probabilities for variables in group g.
  #   sigma  : Error standard deviation. In ordinal probit, usually fixed at 1.
  #
  # Outputs:
  #   logzj  : Log Bayes factor log(Z_g).
  # --------------------------------------------------------------------------
  
  Xg <- as.matrix(Xg)
  Res <- as.numeric(Res)
  
  n <- nrow(Xg)
  p <- ncol(Xg)
  
  if (p == 0) {
    return(list(logzj = 0))
  }
  
  if (length(Res) != n) stop("length of Res must be equal to nrow(Xg).")
  if (sigma <= 0) stop("sigma must be positive.")
  
  tau_g <- rep(tau_g, length.out = p)
  prob  <- rep(prob,  length.out = p)
  
  if (any(tau_g <= 0)) stop("tau_g must be positive.")
  
  # Keep prob in (0, 1) to avoid log(0).
  eps <- 1e-12
  prob <- pmin(pmax(prob, eps), 1 - eps)
  
  idx_all <- seq_len(p)
  log_terms <- numeric(0)
  
  
  for (m in 1:p) {
    for (S in combn(idx_all, m, simplify = FALSE)) {
      S <- unlist(S)
      
      # Log prior weight:
      # prod_{j in S} (1 - prob_j) * prod_{j notin S} prob_j
      logw <- sum(log1p(-prob[S])) + sum(log(prob[setdiff(idx_all, S)]))
      Xk   <- Xg[, S, drop = FALSE]
      
      # Q = X_S' X_S / sigma^2 + T_S^{-1}
      Tinv <- diag(1 / (tau_g[S]^2), nrow = length(S))
      Q <- crossprod(Xk) / sigma^2 + Tinv
      
      # b = X_S' Res / sigma^2
      b_vec <- crossprod(Xk, Res) / sigma^2
      
      # Use Cholesky decomposition for numerical stability.
      R <- chol(Q)                                  
      ytmp  <- forwardsolve(t(R),  b_vec)                 
      Qi_b  <- backsolve(R, ytmp) 
      
      logdetQ <- 2 * sum(log(diag(R)))
      
      # log term for subset S:
      # logw + 1/2 b'Q^{-1}b - 1/2 log|Q| - 1/2 log|T_S|  
      log_term <- 0.5 * drop(t(b_vec) %*% Qi_b) - 0.5 * logdetQ - 0.5 * sum(log(tau_g[S]^2))
      log_terms <- c(log_terms, logw + log_term)
    }
  }
  
  # Empty subset term:
  log_base <- sum(log(prob))
  
  # logZ = log(exp(log_base) + sum_S exp(log_terms_S))
  logZ <- log_sum_exp(c(log_base, log_terms))
  
  return(list(logzj = as.numeric(logZ)))
}

cgs_transit <- function(
  basis, 
  Y_data, 
  beta_init = NULL, 
  r_init = NULL,
  slab_sd = 1,
  rho = 0.5,
  iteration = 200){
  # --------------------------------------------------------------------------
  # Purpose
  # --------------------------------------------------------------------------
  # Update variable indicators and regression coefficients within one active
  # group by a component-wise Gibbs sampler.
  #
  # This function treats Y_data as the partial residual for the selected group.
  # The error variance is fixed at 1, as in the ordinal probit latent model.
  #
  # Inputs:
  #   basis      : n x p_g design matrix for one group.
  #   Y_data     : n-vector partial residual for the selected group.
  #   beta_init  : Initial coefficient values for variables in the group.
  #   r_init     : Initial variable indicators for variables in the group.
  #   slab_sd    : Slab standard deviation for active coefficients.
  #   rho        : Prior exclusion probability for variable indicators.
  #   iteration  : Number of component-wise Gibbs iterations.
  #
  # Outputs:
  #   beta_final : Updated coefficient vector.
  #   r_final    : Updated variable indicator vector.
  # --------------------------------------------------------------------------

  # Basic input checks
  basis  <- as.matrix(basis)
  Y_data <- as.numeric(Y_data)
  
  n <- nrow(basis)
  p <- ncol(basis)
  
  if (length(Y_data) != n) {
    stop("length of Y_data must be equal to nrow(basis).")
  }
  
  # Remove zero-norm columns and fill them back after updating.
  keep_idx <- which(colSums(basis^2) > 0)
  
  if (length(keep_idx) < p) {
    warning("there is at least one regressor whose norm is equal to zero.")
  }
    
  X <- as.matrix(basis[, keep_idx, drop = FALSE])
  p2 <- ncol(X)
  
  slab_sd  <- rep(slab_sd,  length.out = p)[keep_idx]
  rho  <- rep(rho,  length.out = p)[keep_idx]
  
  # Initialization
  beta <- if (is.null(beta_init)) numeric(p2) else as.numeric(beta_init[keep_idx])
  r    <- if (is.null(r_init))    numeric(p2) else as.numeric(r_init[keep_idx])
  
  # Error standard deviation is fixed at 1 in the ordinal probit latent model.
  Sigma_error <- 1
  
  # Squared norm of each predictor column.
  basis_norm <- colSums(X^2)
  
  # Component-wise Gibbs sampler
  for (counter in seq_len(iteration)) {

      for (i in seq_len(p2)) {
        
        # Partial residual:
        # Y_data - X beta + x_i beta_i = Y_data - X_{-i} beta_{-i}
        Res <- as.numeric(Y_data - X %*% beta + X[, i] * beta[i])
       
        denom      <- basis_norm[i] * slab_sd[i]^2 + Sigma_error^2
        
        sigma2_star  <- (Sigma_error^2 * slab_sd[i]^2) / denom
        sigma_star  <- sqrt(sigma2_star)
        mu_star <- (sum(Res * X[, i]) * slab_sd[i]^2) / denom
        
        # z_inv is the inverse of the Bayes factor G_i.
        z_inv <- sqrt(slab_sd[i]^2 / sigma2_star) * exp(-mu_star^2 / (2 * sigma2_star))
       
        # Posterior inclusion probability:
        p_act      <- (1 - rho[i]) / (rho[i] * z_inv + (1 - rho[i]))  
        
        if (runif(1) > p_act) {
          beta[i] <- 0
          r[i]    <- 0
        } else {
          beta[i] <- mu_star + sigma_star * rnorm(1)
          r[i]    <- 1
        }
     }
   }
    
  # Fill updated values back to the original p-dimensional vectors.
  beta_final <- numeric(p) 
  r_final <- numeric(p)
  
  beta_final[keep_idx] <- beta
  r_final[keep_idx]    <- r
  
  return(list(
    beta_final = beta_final,
    r_final    = r_final
  )) 
}

update_groups_one_sweep <- function(
  X, 
  y_aug, 
  beta, 
  r,
  slab_sd, 
  rho, 
  theta, 
  group_list,
  iter_inn) {
  # --------------------------------------------------------------------------
  # Purpose
  # --------------------------------------------------------------------------
  # Perform one group-wise Gibbs sweep for sparse group selection.
  #
  # For each group, this function first samples the group indicator eta_g.
  # If eta_g = 0, all variables in the group are excluded.
  # If eta_g = 1, gamma_gj and beta_gj within the group are updated by the
  # component-wise Gibbs sampler.
  #
  # Inputs:
  #   X          : n x p full design matrix.
  #   y_aug      : Current latent response vector Y_star.
  #   beta       : Current regression coefficient vector.
  #   r          : Current variable indicator vector.
  #   slab_sd    : Slab standard deviations for all variables.
  #   rho        : Prior exclusion probabilities for variable indicators.
  #   theta      : Prior exclusion probabilities for group indicators.
  #   group_list : List of predictor indices defining groups.
  #   iter_inn   : Number of within-group Gibbs updates.
  #
  # Outputs:
  #   beta       : Updated regression coefficient vector.
  #   r          : Updated variable indicator vector.
  #   eta        : Updated group indicator vector.
  # --------------------------------------------------------------------------
  X <- as.matrix(X)
  y_aug <- as.numeric(y_aug)

  n <- nrow(X)
  p <- ncol(X)
  G <- length(group_list)
  
  # Basic input checks
  if (length(y_aug) != n) {
    stop("length of y_aug must be equal to nrow(X).")
  }
  
  if (length(beta) != p || length(r) != p) {
    stop("length of beta and r must be equal to ncol(X).")
  }
  
  if (G == 0) {
    stop("group_list must contain at least one group.")
  }
  
  all_index <- unlist(group_list)
  
  if (any(all_index < 1) || any(all_index > p)) {
    stop("group_list contains invalid predictor indices.")
  }
  
  
  slab_sd <- rep(slab_sd, length.out = p)
  rho     <- rep(rho,     length.out = p)
  theta   <- rep(theta,   length.out = G)
  eta <- integer(G)
  
  for (g in seq_len(G)) {
    index <- group_list[[g]]
    
    if (length(index) == 0) next
    
    # Partial residual excluding group g:
    Xg  <- X[, index, drop = FALSE]
    Res <- as.numeric(y_aug - X %*% beta + Xg %*% beta[index])
    
    # Compute log Bayes factor log(Z_g)
    out <- compute_group_logZ(
      Xg = Xg, 
      Res = Res, 
      tau_g = slab_sd[index],
      prob = rho[index], 
      sigma = 1
      )
    
    logZ <- out$logzj
    
    # Posterior inclusion probability of group g:
    logit_p <- log(1 - theta[g]) - log(theta[g]) + logZ
    p_act <- plogis(logit_p)
    
    if (is.na(p_act) || !is.finite(p_act)) p_act <- 0
    
    if (runif(1) >= p_act) {
      
      # Inactive group: set all variables in this group to zero.
      eta[g] <- 0
      beta[index] <- 0
      r[index]    <- 0
      
    } else {
      
      # Active group: update variables within the group.
      res_inner <- cgs_transit(
        basis      = Xg,
        Y_data     = Res,
        beta_init  = beta[index],
        r_init     = r[index],
        slab_sd   = slab_sd[index],
        rho        = rho[index],
        iteration  = iter_inn,
      )
      eta[g] <- 1
      beta[index] <- res_inner$beta_final
      r[index]    <- res_inner$r_final
    }
  }
  return(list(
    beta = beta,
    r    = r,
    eta  = eta
  ))
}

fit_ordinal_sparse_group <- function(
  fitting_data, 
  iter = 3000, 
  iter_inn = 2000,
  a = Inf, 
  b = 10, 
  theta = 0.5, 
  rho = 0.5, 
  seed_num = 1, 
  group_list = list(), 
  use_mcse = TRUE, 
  check_every = 1000, 
  target_hw = 0.06, 
  burn_proposed = 1000, 
  extra_keep = 2000,
  print_every = 1000,
  verbose = TRUE){
  # --------------------------------------------------------------------------
  # Purpose
  # --------------------------------------------------------------------------
  # Fit an ordinal probit regression model with indicator-based Bayesian
  # variable selection or sparse group selection.
  #
  # Data format:
  #   fitting_data is an n x (p + 1) data frame or matrix.
  #   The first p columns are predictors, and the last column is the ordinal
  #   response coded as 1, ..., K.
  #
  # If group_list = list(), ordinary variable selection is used.
  # If group_list is provided, sparse group selection is used.
  #
  # The MCMC sampler updates beta, gamma, tau, latent responses Y_star, and
  # group indicators eta when group_list is provided.
  #
  # --------------------------------------------------------------------------
  # Inputs
  # --------------------------------------------------------------------------
  # fitting_data   : Data matrix. First p columns are predictors; last column is Y.
  # iter           : Maximum number of MCMC iterations.
  # iter_inn       : Number of within-group Gibbs updates for active groups.
  # a              : Cutpoint prior/proposal sd. If Inf, use uniform update.
  # b              : Slab sd for active regression coefficients.
  # theta          : Prior exclusion probability for selection indicators.
  # rho            : Prior exclusion probability for variables within active groups.
  # seed_num       : Random seed.
  # group_list     : List of predictor indices defining groups.
  # use_mcse       : Whether to use MCSE stopping rule.
  # check_every    : Frequency for checking MCSE.
  # target_hw      : Target maximum MCSE half-width.
  # burn_proposed  : Initial burn-in before checking MCSE.
  # extra_keep     : Extra samples kept after MCSE criterion is reached.
  # print_every    : Frequency for printing progress.
  # verbose        : Whether to print progress.
  #
  # --------------------------------------------------------------------------
  # Outputs
  # --------------------------------------------------------------------------
  # A list with:
  #   gamma_record : p x S posterior samples of variable indicators.
  #   beta_record  : p x S posterior samples of regression coefficients.
  #   tau_record   : (K - 1) x S posterior samples of cutpoints.
  #   eta_record   : G x S posterior samples of group indicators, if grouped.
  #   burn_in      : Final burn-in iteration.
  #   last_iter    : Last MCMC iteration used.
  #
  # Here S is the number of retained posterior samples after burn-in.
  # --------------------------------------------------------------------------
  
  if (!"msm" %in% loadedNamespaces()){
    library(msm)
  }
  if (!"MASS" %in% loadedNamespaces()){
    library(MASS)
  }
  if (!is.numeric(fitting_data[[ncol(fitting_data)]])) {
    stop("The last column of `fitting_data` must be the ordinal response coded as 1, ..., K.")
  }
  Y_data <- as.numeric(fitting_data[[ncol(fitting_data)]])
  if (!all(sort(unique(Y_data)) == seq_len(max(Y_data)))) {
    stop("The ordinal response must be coded as consecutive integers: 1, ..., K.")
  }

  set.seed(seed_num)
  t0 <- Sys.time() 
  final_burn    <- NA_integer_
  stop_iter     <- iter   
  
  # Initialize data dimensions and model parameters
  n_size = nrow(fitting_data)          
  p_size = ncol(fitting_data)-1         
  X_data = as.matrix(fitting_data[, 1:p_size])
  K = max(Y_data)                     
  beta = rep(0, p_size)
  gamma = rep(0, p_size)
  
  if (a == Inf){
    tau = sort(runif(K-1, 0, 1))
  }else{
    tau = sort(rnorm(K-1, mean = 0, sd = a))
  }
  
  # Allocate storage for posterior samples
  beta_record = array(0, dim = c(p_size, iter))
  gamma_record = array(0, dim = c(p_size, iter))
  tau_record = array(0, dim = c(K-1, iter))
  
  if (length(group_list) > 0) {
    eta_record <- array(0, dim = c(length(group_list), iter))
  }else{
    eta_record <- NULL 
  }
  
  # Initial values of latent variables
  Y_star_data = rep(0, n_size)
  
  for (i in 1:n_size){
    if (Y_data[i] == 1){
      Y_star_data[i] = -abs(rnorm(1)) + tau[1]
    }
    else if (Y_data[i] == K){
      Y_star_data[i] = abs(rnorm(1)) + tau[K-1]
    }
    else{
      lower = tau[Y_data[i]-1]
      upper = tau[Y_data[i]]
      Y_star_data[i] = runif(1, lower, upper)
    }
  }
  temp_start = 1
  
  # Main MCMC loop
  for (counter in temp_start:iter) {
    # Draw gamma, beta
    if (length(group_list) == 0){
      for (j in 1:p_size){
        R_j = Y_star_data - as.numeric(X_data %*% beta) + X_data[, j] * beta[j]
        mu_tilde_j = b^2 * sum(R_j * X_data[, j]) / (1 + b^2 * sum(X_data[, j] * X_data[, j]))
        sigma_2_tilde_j = b^2 / (1 + b^2 * sum(X_data[, j] * X_data[, j]))
        sigma_tilde_j = sqrt(sigma_2_tilde_j)
        G_j = sigma_tilde_j / b * exp(1/2 * mu_tilde_j^2 / sigma_2_tilde_j)
        # Posterior inclusion probability of variable j
        if (G_j > 10^20){
          p_j = 1
        }
        else{
          p_j = (1 - theta) * G_j / ((1 - theta) * G_j + theta)
        }
        # Sample gamma_j and beta_j
        unif = runif(1)
        if (p_j < unif){
          beta[j] = 0
          gamma[j] = 0
        }
        else{
          beta[j] = rnorm(1, mu_tilde_j, sigma_tilde_j)
          gamma[j] = 1
        }
      }
    } else {
      # Draw eta, gamma, and beta
      res_gw <- update_groups_one_sweep(
        X        = X_data,
        y_aug    = Y_star_data,
        beta     = beta,
        r        = gamma,         
        slab_sd  = rep(b,  length.out = p_size), 
        rho      = rep(rho, length.out = p_size), 
        theta    = rep(theta, length(group_list)),
        group_list = group_list,
        sigma_w  = rep(1, n_size),
        Sigma_error = 1,
        iter_inn  = iter_inn
      )
      beta  <- res_gw$beta
      gamma <- res_gw$r
      eta   <- res_gw$eta
      eta_record[, counter] <- eta
    }
    # Draw tau
    if (a == Inf){
      for (k in 1:(K-1)){
        lower <- if (any(Y_data == k))     max(Y_star_data[Y_data == k])     else -Inf
        upper <- if (any(Y_data == k + 1)) min(Y_star_data[Y_data == k + 1]) else  Inf
        if (lower < upper) {
          tau[k] <- runif(1, lower, upper)
        } else {
          tau[k] <- tau[k] 
        }
      }
    }
    else{
      tau_star = rep(0, (K-1))
      for (k in 1:(K-1)){
        if (k == 1){
          lower = -10^20
          upper = tau[k+1]
        }
        else if (k == (K-1)){
          lower = tau_star[k-1]
          upper = 10^20
        }
        else{
          lower = tau_star[k-1]
          upper = tau[k+1]
        }
        tau_star[k] = rtnorm(1, mean=tau[k], sd=a, lower=lower, upper=upper)
      }
      logr <- - (sum(tau_star^2) - sum(tau^2)) / (2 * a^2)   
      # Accept or reject the proposed cutpoints
      for (i in 1:n_size) {
        mu <- sum(X_data[i,] * beta)
        
        if (Y_data[i] == 1) {
          p_star <- pnorm(tau_star[1], mean = mu, sd = 1)
          p_now  <- pnorm(tau[1],      mean = mu, sd = 1)
          
        } else if (Y_data[i] == K) {
          p_star <- 1 - pnorm(tau_star[K-1], mean = mu, sd = 1)
          p_now  <- 1 - pnorm(tau[K-1],      mean = mu, sd = 1)
          
        } else {
          lo_s <- tau_star[Y_data[i]-1]; up_s <- tau_star[Y_data[i]]
          lo   <- tau[Y_data[i]-1];      up   <- tau[Y_data[i]]
          
          p_star <- pnorm(up_s, mean = mu, sd = 1) - pnorm(lo_s, mean = mu, sd = 1)
          p_now  <- pnorm(up,   mean = mu, sd = 1) - pnorm(lo,   mean = mu, sd = 1)
        }
        # avoid log(0)
        p_star <- max(p_star, 1e-300)
        p_now  <- max(p_now,  1e-300)
        
        logr <- logr + (log(p_star) - log(p_now))
      }
      
      if (log(runif(1)) < min(0, logr)) {
        tau <- tau_star
      }
      
      
    }
    # Draw Y_star
    for (i in 1:n_size){
      mean_Y_star = sum(X_data[i,] * beta)
      if (Y_data[i] == 1){
        Y_star_data[i] = rtnorm(1, mean=mean_Y_star, sd=1, lower=-Inf, upper=tau[1])
      }
      else if (Y_data[i] == K){
        Y_star_data[i] = rtnorm(1, mean=mean_Y_star, sd=1, lower=tau[K-1], upper=Inf)
      }
      else{
        lower = tau[Y_data[i]-1]
        upper = tau[Y_data[i]]
        Y_star_data[i] = rtnorm(1, mean=mean_Y_star, sd=1, lower=lower, upper=upper)
      }
    }
    
    # Store current MCMC samples
    beta_record[, counter] = beta
    gamma_record[, counter] = gamma
    tau_record[, counter] = tau
    
    # print diagnostics
    if (counter %% print_every == 0 && verbose == TRUE){
      
      cat(sprintf("iter %d / %d | tau = [%s] | elapsed %.1f mins\n",
                  counter,iter,
                  paste(sprintf("%.3f", tau), collapse = ", "),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    }
  
    # MCSE Convergence Check
    if(use_mcse){
      if (counter >= burn_proposed + 50 && counter %% check_every == 0) {
        idx <- (burn_proposed + 1):counter   
        gamma_mat <- gamma_record[, idx, drop = FALSE]
        res_g   <- mcse_max_halfwidth(gamma_mat, alpha = 0.05)
        if (length(group_list) > 0){
          eta_mat   <- eta_record[,   idx, drop = FALSE]
          res_eta <- mcse_max_halfwidth(eta_mat,   alpha = 0.05)
          max_hw_now <- max(res_g$max_hw, res_eta$max_hw)
        }else{
          max_hw_now <- res_g$max_hw
        }
        # First time arriving the standard. Determine final burn-in & set stop position
        if (is.na(final_burn) && max_hw_now < target_hw) {
          stop_iter  <- min(iter, counter + extra_keep)
          final_burn <- counter
          print(sprintf("maximum MCSE of eta and gamma is %.6f", max_hw_now))
          print(sprintf("MCSE criterion reached at iter = %d",  final_burn))
        }
      }
      if (!is.na(final_burn) && counter >= stop_iter) {
        last_iter <- counter
        break
      }
    }
  }
  # Keep post-burn-in samples
  if (is.na(final_burn)) {
    final_burn <- burn_proposed   
    last_iter  <- iter
  }
  kept_idx <- (final_burn + 1):last_iter
  result <- list(
    gamma_record = gamma_record[, kept_idx, drop = FALSE],
    beta_record  = beta_record[,  kept_idx, drop = FALSE],
    tau_record   = tau_record[,   kept_idx, drop = FALSE],
    burn_in      = final_burn,
    last_iter    = last_iter
  )
  if (!is.null(eta_record)) {
    result$eta_record <- eta_record[, kept_idx, drop = FALSE]
  }
  return(result)
}



