# =============================================================================
# Accelerated Bayesian sparse group selection for ordinal probit models
# =============================================================================
# This file implements the accelerated version of the proposed Bayesian sparse
# group selection method.
#
# The main function is:
#   fit_ordinal_sparse_group_fast()
#
# Supporting functions include:
#   compute_group_logBF_allin()    : computes all-in group-level log Bayes factors.
#   update_active_group_reparam()  : updates variables within an active group.
# =============================================================================


compute_group_logBF_allin <- function(Xg, Res_g, slab_sd_g) {
  # --------------------------------------------------------------------------
  # Compute the all-in group-level log Bayes factor for updating eta_g.
  #
  # This version treats all variables in the group as active when computing the
  # group-level Bayes factor.
  # --------------------------------------------------------------------------
  
  Xg    <- as.matrix(Xg)
  Res_g <- as.numeric(Res_g)
  
  p_g <- ncol(Xg)
  
  if (p_g == 0) {
    return(list(logBF = 0))
  }
  
  slab_sd_g <- pmax(slab_sd_g, 1e-12)
  slab_var  <- slab_sd_g^2
  
  A    <- crossprod(Xg) + diag(1 / slab_var, nrow = p_g)
  bvec <- crossprod(Xg, Res_g)
  
  Rchol <- chol(A)
  tmp   <- forwardsolve(t(Rchol), bvec)
  Ainvb <- backsolve(Rchol, tmp)
  
  logdetA <- 2 * sum(log(diag(Rchol)))
  quad    <- drop(crossprod(bvec, Ainvb))
  
  logBF <- 0.5 * quad - 0.5 * logdetA - sum(log(slab_sd_g))
  
  return(list(logBF = as.numeric(logBF)))
}


update_active_group_reparam <- function(
  X, 
  y_aug, 
  beta, 
  gamma, 
  tilde_beta,
  xbeta, 
  g_idx, 
  slab_sd, 
  rho
) {
  # --------------------------------------------------------------------------
  # Update variable indicators and coefficients within one active group.
  #
  # The current linear predictor xbeta is updated incrementally to avoid
  # repeatedly computing X %*% beta.
  # --------------------------------------------------------------------------
  
  for (j in g_idx) {
    
    xj <- X[, j]
    
    # Partial residual excluding variable j
    R_j <- y_aug - xbeta + xj * beta[j]
    
    # Conditional distribution of tilde_beta_j given gamma_j = 1
    xTx <- drop(crossprod(xj))
    
    sigma2_tilde_j <- 1 / (xTx + 1 / slab_sd[j]^2)
    mu_tilde_j     <- sigma2_tilde_j * drop(crossprod(xj, R_j))
    
    # Posterior inclusion probability of gamma_j
    log_zj <- 
      0.5 * log(sigma2_tilde_j) -
      log(slab_sd[j]) +
      0.5 * mu_tilde_j^2 / sigma2_tilde_j
    
    logit_p <- log1p(-rho[j]) + log_zj - log(rho[j])
    p_j <- plogis(logit_p)
    
    old_beta_j <- beta[j]
    
    if (runif(1) < p_j) {
      gamma[j]      <- 1
      tilde_beta[j] <- rnorm(1, mean = mu_tilde_j, sd = sqrt(sigma2_tilde_j))
      beta[j]       <- tilde_beta[j]
    } else {
      gamma[j]      <- 0
      tilde_beta[j] <- 0
      beta[j]       <- 0
    }
    
    # Incremental update of the linear predictor
    xbeta <- xbeta + xj * (beta[j] - old_beta_j)
  }
  
  return(list(
    beta       = beta,
    gamma      = gamma,
    tilde_beta = tilde_beta,
    xbeta      = xbeta
  ))
}


fit_ordinal_sparse_group_fast <- function(
  fitting_data,
  iter = 3000,
  a = 1,
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
  verbose = TRUE) {
  # --------------------------------------------------------------------------
  # Fit the accelerated Bayesian sparse group selection model for ordinal
  # probit regression.
  #
  # This version updates the group indicator eta_g using an all-in group-level
  # Bayes factor and updates variables within active groups by an incremental
  # component-wise Gibbs update.
  #
  # Data format:
  #   fitting_data is an n x (p + 1) data frame or matrix.
  #   The first p columns are predictors, and the last column is the ordinal
  #   response coded as 1, ..., K.
  #
  # Main inputs:
  #   iter        : Maximum number of MCMC iterations.
  #   a           : Cutpoint proposal standard deviation.
  #   b           : Slab standard deviation for active coefficients.
  #   theta       : Prior exclusion probability for group indicators.
  #   rho         : Prior exclusion probability for variables within active groups.
  #   group_list  : List of predictor indices defining groups.
  #   use_mcse    : Whether to use the MCSE stopping rule.
  #
  # Output:
  #   A list containing posterior samples of gamma, beta, tilde_beta, tau, eta,
  #   and information about burn-in and stopping iteration.
  # --------------------------------------------------------------------------  
  
  if (!is.numeric(fitting_data[[ncol(fitting_data)]])) {
    stop("Y 必須是最後一欄且為數值類別 1..K")
  }
  if (any(!fitting_data[[ncol(fitting_data)]] %in%
          1:max(fitting_data[[ncol(fitting_data)]]))) {
    stop("Y 必須是 1..K")
  }
  if (length(group_list) == 0) {
  stop("group_list must be provided for sparse group selection.")
  }
  t0 <- Sys.time()
  final_burn     <- NA_integer_
  stop_iter      <- iter
  mcse_iteration <- NA_integer_
  
  set.seed(seed_num)
  
  if (!"msm" %in% loadedNamespaces()) library(msm)
  
  # -------------------------------------------------------
  # basic objects
  # -------------------------------------------------------
  n_size <- nrow(fitting_data)
  p_size <- ncol(fitting_data) - 1
  
  X_data <- as.matrix(fitting_data[, 1:p_size])
  Y_data <- as.numeric(fitting_data[, p_size + 1])
  K      <- max(Y_data)
  
  # ensure each variable belongs to some group
  if (length(group_list) != 0) {
    existing_numbers <- unlist(group_list)
    missing_numbers  <- setdiff(1:p_size, existing_numbers)
    for (num in missing_numbers) {
      group_list <- append(group_list, list(c(num)))
    }
  }
  g_size <- length(group_list)
  
  # prior hyperparameters
  slab_sd <- if (length(b) == 1) rep(b, p_size) else b
  if (length(slab_sd) != p_size) stop("b 必須是長度 1 或 p 的向量")
  
  rho_vec <- if (length(rho) == 1) rep(rho, p_size) else rho
  if (length(rho_vec) != p_size) stop("rho 必須是長度 1 或 p 的向量")
  rho_vec <- clamp_prob(rho_vec)
  
  theta_vec <- if (length(theta) == 1) rep(theta, max(1, g_size)) else theta
  if (g_size > 0 && length(theta_vec) != g_size) {
    stop("grouped case 下，theta 必須是長度 1 或 group 數的向量")
  }
  if (g_size > 0) theta_vec <- clamp_prob(theta_vec)
  
  # -------------------------------------------------------
  # initialize states
  # -------------------------------------------------------
  beta       <- rep(0, p_size)
  tilde_beta <- rep(0, p_size)
  gamma      <- rep(0, p_size)
  eta        <- if (g_size > 0) rep(0, g_size) else NULL
  
  if (is.infinite(a)) {
    tau <- sort(runif(K - 1, 0, 1))
  } else {
    tau <- sort(rnorm(K - 1, mean = 0, sd = a))
  }
  
  # -------------------------------------------------------
  # storage
  # -------------------------------------------------------
  beta_record       <- array(0, dim = c(p_size, iter))
  tildebeta_record  <- array(0, dim = c(p_size, iter))
  gamma_record      <- array(0, dim = c(p_size, iter))
  tau_record        <- array(0, dim = c(K - 1, iter))
  eta_record        <- if (g_size > 0) array(0, dim = c(g_size, iter)) else NULL
  
  # -------------------------------------------------------
  # initialize latent Y*
  # -------------------------------------------------------
  Y_star_data <- rep(0, n_size)
  for (i in 1:n_size) {
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
  
  # current linear predictor
  xbeta <- as.numeric(X_data %*% beta)
  
  # -------------------------------------------------------
  # main MCMC
  # -------------------------------------------------------
  for (counter in 1:iter) {
    
    # =====================================================
    # update beta / gamma / eta
    # =====================================================
      for (g in seq_len(g_size)) {
        g_idx <- group_list[[g]]
        Xg    <- X_data[, g_idx, drop = FALSE]
        
        # residual excluding group g
        group_contrib <- as.numeric(Xg %*% beta[g_idx])
        R_g <- Y_star_data - xbeta + group_contrib
        
        # --------------------------------------------------
        # Version A: partially-collapsed eta_g update (all-in)
        # --------------------------------------------------
        bf_out <- compute_group_logBF_allin(
          Xg        = Xg,
          Res_g     = R_g,
          slab_sd_g = slab_sd[g_idx]
        )
        
        # P(eta_g = 1 | ...) = (1-theta)*BF / ((1-theta)*BF + theta)
        logit_eta <- log1p(-theta_vec[g]) + bf_out$logBF - log(theta_vec[g])
        p_eta     <- plogis(logit_eta)
        
        old_beta_g <- beta[g_idx]
        
        if (runif(1) < p_eta) {
          eta[g] <- 1
          
          upd <- update_active_group_reparam(
            X          = X_data,
            y_aug      = Y_star_data,
            beta       = beta,
            gamma      = gamma,
            tilde_beta = tilde_beta,
            xbeta      = xbeta,
            g_idx      = g_idx,
            slab_sd    = slab_sd,
            rho        = rho_vec
          )
          
          beta       <- upd$beta
          gamma      <- upd$gamma
          tilde_beta <- upd$tilde_beta
          xbeta      <- upd$xbeta
          
        } else {
          eta[g] <- 0
          
          # switch off the entire group
          beta[g_idx]       <- 0
          tilde_beta[g_idx] <- 0
          gamma[g_idx]      <- 0
          
          # update xbeta by removing old group contribution
          xbeta <- xbeta - as.numeric(Xg %*% old_beta_g)
        }
      }
    
    
    # =====================================================
    # update cutpoints tau (Chang MH / uniform version)
    # =====================================================
    if (is.infinite(a)) {
      
      for (k in 1:(K - 1)) {
        lower <- if (any(Y_data == k))     max(Y_star_data[Y_data == k])     else -Inf
        upper <- if (any(Y_data == k + 1)) min(Y_star_data[Y_data == k + 1]) else  Inf
        
        if (lower < upper) {
          tau[k] <- runif(1, lower, upper)
        } else {
          tau[k] <- tau[k]
        }
      }
      
    } else {
      tau_star <- rep(0, K - 1)
      
      for (k in 1:(K - 1)) {
        if (k == 1) {
          lower <- -1e20
          upper <- tau[k + 1]
        } else if (k == (K - 1)) {
          lower <- tau_star[k - 1]
          upper <-  1e20
        } else {
          lower <- tau_star[k - 1]
          upper <- tau[k + 1]
        }
        tau_star[k] <- msm::rtnorm(1, mean = tau[k], sd = a,
                                   lower = lower, upper = upper)
      }
      
      logr <- -(sum(tau_star^2) - sum(tau^2)) / (2 * a^2)
      
      for (i in 1:n_size) {
        mu <- xbeta[i]
        
        if (Y_data[i] == 1) {
          p_star <- pnorm(tau_star[1], mean = mu, sd = 1)
          p_now  <- pnorm(tau[1],      mean = mu, sd = 1)
          
        } else if (Y_data[i] == K) {
          p_star <- 1 - pnorm(tau_star[K - 1], mean = mu, sd = 1)
          p_now  <- 1 - pnorm(tau[K - 1],      mean = mu, sd = 1)
          
        } else {
          lo_s <- tau_star[Y_data[i] - 1]
          up_s <- tau_star[Y_data[i]]
          lo   <- tau[Y_data[i] - 1]
          up   <- tau[Y_data[i]]
          
          p_star <- pnorm(up_s, mean = mu, sd = 1) - pnorm(lo_s, mean = mu, sd = 1)
          p_now  <- pnorm(up,   mean = mu, sd = 1) - pnorm(lo,   mean = mu, sd = 1)
        }
        
        p_star <- max(p_star, 1e-300)
        p_now  <- max(p_now,  1e-300)
        
        logr <- logr + log(p_star) - log(p_now)
      }
      
      if (log(runif(1)) < min(0, logr)) {
        tau <- tau_star
      }
    }
    
    # =====================================================
    # update latent Y*
    # =====================================================
    for (i in 1:n_size) {
      mean_Y_star <- xbeta[i]
      
      if (Y_data[i] == 1) {
        Y_star_data[i] <- msm::rtnorm(1, mean = mean_Y_star, sd = 1,
                                      lower = -Inf, upper = tau[1])
      } else if (Y_data[i] == K) {
        Y_star_data[i] <- msm::rtnorm(1, mean = mean_Y_star, sd = 1,
                                      lower = tau[K - 1], upper = Inf)
      } else {
        lower <- tau[Y_data[i] - 1]
        upper <- tau[Y_data[i]]
        Y_star_data[i] <- msm::rtnorm(1, mean = mean_Y_star, sd = 1,
                                      lower = lower, upper = upper)
      }
    }
    
    # =====================================================
    # record
    # =====================================================
    beta_record[, counter]      <- beta
    tildebeta_record[, counter] <- tilde_beta
    gamma_record[, counter]     <- gamma
    tau_record[, counter]       <- tau
    if (g_size > 0) eta_record[, counter] <- eta
    
    # print
    if (counter %% print_every == 0 && verbose == TRUE) {
      cat(sprintf(
        "iter %d / %d | tau = [%s] | elapsed %.1f mins\n",
        counter, iter,
        paste(sprintf("%.3f", tau), collapse = ", "),
        as.numeric(difftime(Sys.time(), t0, units = "mins"))
      ))
    }
  
    # =====================================================
    # MCSE stopping
    # =====================================================
    if (use_mcse) {
      if (counter >= burn_proposed + 50 && counter %% check_every == 0) {
        idx <- (burn_proposed + 1):counter
        
        gamma_mat <- gamma_record[, idx, drop = FALSE]
        res_g     <- mcse_max_halfwidth(gamma_mat, alpha = 0.05)
        
        if (g_size > 0) {
          eta_mat <- eta_record[, idx, drop = FALSE]
          res_eta <- mcse_max_halfwidth(eta_mat, alpha = 0.05)
          max_hw_now <- max(res_g$max_hw, res_eta$max_hw)
        } else {
          max_hw_now <- res_g$max_hw
        }
        
        if (is.na(final_burn) && max_hw_now < target_hw) {
          final_burn     <- counter
          stop_iter      <- min(iter, counter + extra_keep)
          mcse_iteration <- counter
          if(verbose){
          print(sprintf("maximum MCSE of eta and gamma is %.6f", max_hw_now))
          print(sprintf("MCSE criterion reached at iter = %d", mcse_iteration))
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
  }
  
  kept_idx <- (final_burn + 1):last_iter
  
  result <- list(
    gamma_record     = gamma_record[, kept_idx, drop = FALSE],
    beta_record      = beta_record[, kept_idx, drop = FALSE],
    tildebeta_record = tildebeta_record[, kept_idx, drop = FALSE],
    tau_record       = tau_record[, kept_idx, drop = FALSE],
    eta_record       = if (g_size > 0) eta_record[, kept_idx, drop = FALSE] else NULL,
    burn_in          = final_burn,
    last_iter        = last_iter,
    mcse_iteration   = mcse_iteration
  )
  return(result)
}




