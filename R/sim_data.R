#' Simulate RCT and OS Data
#'
#' Simulate randomized clinical trial and observational data
#' with shared properties.
#'
#' @param n_r number of observations in RCT data
#' @param n_o number of observations in OS data
#' @param p total number of covariates to create
#' @param support_fraction percent of covariates associated with observed betas
#' @param covariate_effect error in the outcome model. In the covariate mismatch
#' setting this should be a named vector with values for \sQuote{Z} and \sQuote{V}.
#' @param delta_fraction percent of betas to shift for RCT
#' @param treatment_prop percent of OS to receive treatment
#' @param shift_x set to TRUE to introduce covariate shift
#' @param shift_y set to TRUE to introduce outcome shift
#' @param model function to model outcome; defaults to a linear model. Can
#' also be quadratic, sine, or user specified.
#' @param percent_modifier_to_cover when greater than zero, use missing
#' effect modifier setting; incompatible with covariate mismatch setting
#' @param frac_U fraction of RCT data not found in OS data
#' @param frac_V fraction of OS data not found in RCT data
#'
#' @return
#' \item{simulated data, components and settings}{
#' \describe{
#' \item{X_RCT}{RCT covariates}
#' \item{X_OS}{OS covariates}
#' \item{components}{components created in covariate mismatch settting}
#' \item{a_r}{treatment assignment for RCT data}
#' \item{y_r}{0/1 outcome for RCT data}
#' \item{a_o}{treatment assignment for OS data}
#' \item{y_o}{0/1 outcome for OS data}
#' \item{beta_RCT}{betas for RCT data}
#' \item{beta_OS}{betas for OS data}
#' \item{idx}{indices used to track U/Z/V in covariate mismatch settting}
#' \item{exp_settings}{preserved function argument values}
#' }}
#'
#' @examples
#' x1 <- simulate_rct_and_os_data(
#'   n_r = 300, n_o = 1000, p = 100,
#' )
#' # simulate outcome with quadratic function
#' x2 <- simulate_rct_and_os_data(
#'   n_r = 300, n_o = 1000, model = 'quadratic'
#' )
#' # simulate missing effect modifiers
#' x3 <- simulate_rct_and_os_data(
#'   n_r = 300, n_o = 1000, percent_modifier_to_cover = 0.2
#' )
#' # simulate with covariate mismatch
#' x4 <- simulate_rct_and_os_data(
#'   n_r = 300, n_o = 1000, support_fraction = 1/20,
#'   covariate_effect = c(Z = 2/3, V = 1),
#'   frac_U = 0.3, frac_V = 0.3
#' )
#'
#' @export

simulate_rct_and_os_data <- function(
  n_r = 1000,
  n_o = 10000,
  p = 100,
  support_fraction = 1/10,
  covariate_effect = 2/3,
  delta_fraction = 2/100,
  treatment_prop = 1/3,
  shift_x = TRUE,
  shift_y = TRUE,
  model = c('linear','quadratic','sinusoidal'),
  percent_modifier_to_cover = NULL,
  frac_U = NULL,
  frac_V = NULL
){
  # save arguments
  exp_settings <- as.list(environment())
  # w.r.t. covariate mismatch
  # U: only in RCT
  # V: only in OS
  # Z: shared variables
  frac_U <- frac_U %||% 0
  frac_V <- frac_V %||% 0
  cv_mismatch <- frac_U > 0 || frac_V > 0
  if(cv_mismatch) {
    if (frac_U < 0 || frac_V < 0 || frac_U + frac_V > 1) {
      stop("fractions must be ≥0 and sum to ≤1")
    }
    cov_effect_Z <- unname(covariate_effect['Z'])
    cov_effect_V <- unname(covariate_effect['V'])
    if(is.na(cov_effect_Z) || is.na(cov_effect_V)) {
      stop("covariate_effect must have Z and V set")
    }
    p_u <- floor(frac_U * p)
    p_v <- floor(frac_V * p)
    p_z <- p - p_u - p_v
    idx_u <- seq.int(1, length.out = p_u)
    idx_z <- seq.int(p_u + 1, length.out = p_z)
    idx_v <- seq.int(p_u + p_z + 1, length.out = p_v)
    rct_size <- p_u + p_z
    rct_idx <- c(idx_u, idx_z)
    # turn off missing effect modifiers?
    percent_modifier_to_cover <- NULL
  } else {
    rct_size <- p
    cov_effect_Z <- covariate_effect
  }

  if(is.function(model)) {
    mf <- model
  } else {
    mf <- switch(match.arg(model),
      linear = linearFun,
      quadratic = quadraticFun,
      sinusoidal = sineFun
    )
  }

  # Parameters
  if(cv_mismatch) {
    rho <- 0.4
    Sigma <- matrix(rho, p, p)
    diag(Sigma) <- 1

    sup_Z <- seq_len(ceiling(support_fraction * rct_size))
    # pu+pv OR pz+pv?
    sup_V <- seq_len(ceiling(support_fraction * (p_u+p_v)))
  } else {
    rho <- 0.75
    Sigma <- rho ^ abs(outer(1:p, 1:p, FUN = "-"))
    main_support_size <- seq_len(support_fraction * p)
  }

  mu_o <- numeric(p)
  mu_r <- mu_o
  if(shift_x) {
    ## shift mean of multivariate distribution
    ## sample from unif(-0.5, -0.25) |U| unif(0.25, 0.5)
    if(cv_mismatch) {
      mu_r[idx_u] <- perturb(mu_r[idx_u], idx_u, 0.5)
    } else {
      ## select 10% of covariates to shift mean from zero/mu_o
      mu_r <- perturb(mu_r, sample(p, p/10), 0.5)
    }
  }

  ## zero mean multivariate normal distribution with dense covariance matrix
  X_full_o <- MASS::mvrnorm(n_o, mu_o, Sigma)
  X_full_r <- MASS::mvrnorm(n_r, mu_r, Sigma)

  ## random treatment assignment
  ## at least ~1/3 of samples receive treatment
  base_logit <- -log(1 / treatment_prop - 1)
  if(cv_mismatch) {
    ## carve out observed pieces
    Z_o <- X_full_o[, idx_z, drop = FALSE]
    V_o <- X_full_o[, idx_v, drop = FALSE]

    Z_r <- X_full_r[, idx_z, drop = FALSE]
    U_r <- X_full_r[, idx_u, drop = FALSE]

    p_z10 <- min(p_z, 10)
    lin_term_o <- base_logit + Z_o[, seq_len(p_z10), drop = FALSE] %*% runif(p_z10, .25, .5)
  } else {
    ## logistic model uses 10 random covariates
    #let the first 10 features determine the assignments, the weight of the rest is zero
    lin_term_o <- base_logit + X_full_o %*% c(runif(10, min=.25, max=.5), numeric(max(p-10,0)))[seq_len(p)]
  }
  # linear term=0 => prob=0.5
  e_o <- 1 / (1 + exp(-lin_term_o))
  a_o <- rbinom(n_o, 1, e_o)
  # RCT: 1:1 randomization
  a_r <- rbinom(n_r, 1, 0.5)

  ## ---- coefficient generation for OS arm ------------------------------
  beta0_o <- beta1_o <- numeric(p)
  if(cv_mismatch) {
    ##    Z: fixed effect size (cov_effect_Z)
    beta0_o[rct_idx] <- perturb(beta0_o[rct_idx], sup_Z, cov_effect_Z)
    beta1_o[rct_idx] <- perturb(beta0_o[rct_idx], sup_Z, cov_effect_Z)

    ##    V: user-controlled effect size (cov_effect_V)
    beta0_o[idx_v] <- perturb(beta0_o[idx_v], sup_V, cov_effect_V)
    beta1_o[idx_v] <- perturb(beta0_o[idx_v], sup_V, cov_effect_V)
  } else {
    ## simulate outcome y^o(a)
    ## select 10% (support_fraction) of coefficients associated with beta_observed, arm 0/1
    ## sample from unif(-0.67, -0.33) |U| unif(0.33, 0.67)
    ## noise term is rnorm(0, 1/3) <-- `outcome_model` uses covariate_effect * 0.5 => 1/3
    beta0_o <- perturb(beta0_o, main_support_size, cov_effect_Z)
    ## sample from unif(-0.67, -0.33) |U| unif(0.33, 0.67)
    beta1_o <- perturb(beta0_o, main_support_size, cov_effect_Z)
  }

  ## outcome shift
  beta0_r <- beta0_o
  beta1_r <- beta1_o
  ## shift mean of outcome models y^r(a)
  ## select 2% (delta_fraction) of beta^r_a and perturb
  ## sample from unif(-1, -0.5) |U| unif(0.5, 1)
  if(shift_y){
    #d: size of support of delta
    delta_support_size <- rct_size * delta_fraction
    if(cv_mismatch) {
      beta0_r[rct_idx] <- perturb(beta0_r[rct_idx], sample(rct_size, delta_support_size), 1)
      beta1_r[rct_idx] <- perturb(beta1_r[rct_idx], sample(rct_size, delta_support_size), 1)
    } else {
      beta0_r <- perturb(beta0_r, sample(rct_size, delta_support_size), 1)
      beta1_r <- perturb(beta1_r, sample(rct_size, delta_support_size), 1)
    }
  }

  Y_o <- outcome_model(X_full_o, beta0_o, beta1_o, cov_effect_Z, a_o, mf)
  Y_r <- outcome_model(X_full_r, beta0_r, beta1_r, cov_effect_Z, a_r, mf)

  ## missing effect modifiers
  if(!is.null(percent_modifier_to_cover) && percent_modifier_to_cover > 0 && percent_modifier_to_cover < 1) {
    modifiers_index <- which(beta0_o | beta1_o | beta0_r | beta1_r)
    number_modifier_to_cover <- floor(percent_modifier_to_cover * length(modifiers_index))
    if(number_modifier_to_cover > 0L) {
      index_unobserved_modifiers <- sample(modifiers_index, number_modifier_to_cover)
      x_o_selected <- setdiff(seq_len(ncol(X_full_o)), index_unobserved_modifiers)
      x_r_selected <- setdiff(seq_len(ncol(X_full_r)), index_unobserved_modifiers)
      X_full_o <- X_full_o[,x_o_selected, drop=FALSE]
      X_full_r <- X_full_r[,x_r_selected, drop=FALSE]
    }
  }
  if(cv_mismatch) {
    ## ---------- impute V̂ for the RCT ----------
    ## Fit V ~ Z on OS (one OLS per column; closed-form matrix works too)
    B_hat <- MASS::ginv(t(Z_o) %*% Z_o) %*% (t(Z_o) %*% V_o)   # p_z × p_v
    Vhat_r <- Z_r %*% B_hat                                    # n_r × p_v
    B2_hat <- MASS::ginv(t(Z_r) %*% Z_r) %*% (t(Z_r) %*% U_r)  # p_z × p_v
    Uhat_o <- Z_o %*% B2_hat                                   # n_r × p_v

    components <- list(
      U = U_r,
      V = V_o,
      Z_r = Z_r,
      Z_o = Z_o,
      Vhat = Vhat_r,
      Uhat = Uhat_o,
      U0 = matrix(0, n_o, p_u)
    )
    rct_full <- c('U', 'Z_r', 'Vhat')
    os_full <- c('U0', 'Z_o', 'V')
    X_full_r <- do.call(cbind, components[rct_full])
    X_full_o <- do.call(cbind, components[os_full])
    idx <- list(U = idx_u, Z = idx_z, V = idx_v)
  } else {
    components = NULL
    idx = NULL
  }

  list(
    X_RCT = X_full_r,
    X_OS = X_full_o,
    components = components,
    ## treatments & outcomes
    a_r = a_r, y_r = Y_r,
    a_o = a_o, y_o = Y_o,
    ## bookkeeping
    beta_RCT = list(b0 = beta0_r, b1 = beta1_r),
    beta_OS  = list(b0 = beta0_o, b1 = beta1_o),
    idx = idx,
    exp_settings = exp_settings
  )
}
