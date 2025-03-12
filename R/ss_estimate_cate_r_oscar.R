#' estimate_cate
#'
#' Estimate conditional average treatment effect with R-OSCAR method under sample-split.
#'
#' @keywords internal
#' @param X NEEDINFO
#' @param A NEEDINFO
#' @param Y NEEDINFO
#' @param propensity_vec NEEDINFO
#' @param K NEEDINFO
#' @param outcome_model_obs_coefs NEEDINFO
#' @param m_of_x_mat NEEDINFO
#' @param normalize_weights NEEDINFO

ss_estimate_cate_r_oscar <-
function(X, A, Y, propensity_vec,
                                   K = 5, outcome_model_obs_coefs,
                                   m_of_x_mat = NULL, normalize_weights = TRUE) {
  s_cates_coeffs <- matrix(NA, nrow = ncol(X)+1, ncol = K) #CATEs will be the columns
  # To have a balance of treatment in each split use AR for creating folds.
  index <- caret::createFolds(A,
                       k = 5,
                       list = TRUE,
                       returnTrain = FALSE)
  for (i in seq_len(K)) {
    idx_cur <- index[[i]]
    calibrated_outcome_model_rct <- calibrate_outcome_model(
      X = X[-idx_cur,,drop = FALSE],
      A = A[-idx_cur],
      Y = Y[-idx_cur],
      mu_x_a_coefs = outcome_model_obs_coefs
    )
    cal_mu_x_a_coefs <- calibrated_outcome_model_rct$mu_x_a_coefs
    init_obs_cate_coefs <- cal_mu_x_a_coefs[,2] - cal_mu_x_a_coefs[,1]
    
    X_cur <- X[idx_cur,,drop = FALSE]
    cate_r_oscar <- estimate_cate_r_oscar(
      X = X_cur,
      A = A[idx_cur],
      Y = Y[idx_cur],
      propensity_vec = propensity_vec[idx_cur],
      init_cate_coefs = init_obs_cate_coefs,
      normalize_weights = normalize_weights,
      m_of_x_mat = calibrated_outcome_model_rct$outcome_function(X_cur)
    )
    s_cates_coeffs[, i] <- as.numeric(cate_r_oscar$cate_coefficients)
  }
  averaged_cates_coeffs <- as.matrix(rowMeans(s_cates_coeffs))
  ## function to predict CATE for new data
  return_cate <- function(x)
  {
    unname(drop(cbind(1, x) %*% averaged_cates_coeffs))
  }
  list(
    cate_coefficients = averaged_cates_coeffs,
    cate_function = return_cate,
    fitted_cate = return_cate(X)
  )
}
