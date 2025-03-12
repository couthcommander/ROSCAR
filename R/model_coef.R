#' Compute Model Coefficients
#'
#' CATE estimates
#'
#' @param rct NEEDINFO
#' @param os NEEDINFO
#' @param randomForest NEEDINFO
#' @param sampleSplit NEEDINFO
#'
#' @return
#' \item{estimates by method}{
#' \describe{
#' \item{cate_coefficients}{cate estimate; not available in
#' naive/RACER random forest approach}
#' \item{cate_function}{function to predict CATE for new data}
#' \item{fitted_cate}{fitted CATE values}
#' \item{delta_coef_mat}{difference between coefficient estimates;
#' only available in OSCAR approach}
#' \item{cate_delta_coefficients}{difference between the two CATEs;
#' only available in R-OSCAR approach}
#' \item{cate_model}{CATE model; only available in naive random
#' forest approach}
#' }}
#'
#' @examples
#' data(os, rct)
#' modNoOs <- model_coef(rct)
#' print(modNoOs, max = 40)
#' modOs <- model_coef(rct, os)
#' print(modOs, max = 40)
#'
#' @export

model_coef <-
function(rct, os = NULL, randomForest = FALSE, sampleSplit = FALSE) {
  rctX <- rct$X
  rctY <- rct$Y
  rctA <- rct$A
  propensity_vec <- rep(mean(rctA), length(rctY))

  # 1) fit outcome regression model for each treatment arm in RCT data
  outcome_model_rct <- fit_outcome_model_by_arm(X = rctX, Y = rctY, A = rctA)
  if(randomForest) {
    omRF <- fit_outcome_model_by_arm_rf(X = rctX, Y = rctY, A = rctA)
    muXrf <- omRF$mu_x
  }

  if(!is.null(os)) {
    # 2) fit outcome regression model for each treatment arm in obs data
    outcome_model_obs <- fit_outcome_model_by_arm(X = os$X, Y = os$Y, A = os$A)
    # 3) calibrate the above regression functions to the RCT data
    calibrated_outcome_model_rct <- calibrate_outcome_model(X = rctX,
                                                            A = rctA,
                                                            Y = rctY,
                                                            mu_x_a_coefs = outcome_model_obs$mu_x_a_coefs)
    init_obs_cate_coefs <- calibrated_outcome_model_rct$mu_x_a_coefs[,2] -
      calibrated_outcome_model_rct$mu_x_a_coefs[,1]
  }

  # 4) estimate CATE
  naive_args <- list(
    X = rctX,
    A = rctA,
    Y = rctY,
    propensity_vec = propensity_vec,
    normalize_weights = TRUE,
    m_of_x_mat = NULL
  )
  racer_args <- utils::modifyList(naive_args, list(m_of_x_mat = outcome_model_rct$mu_x))

  if(!is.null(os)) {
    r_oscar_args <- utils::modifyList(naive_args, list(
      init_cate_coefs = init_obs_cate_coefs,
      m_of_x_mat = calibrated_outcome_model_rct$mu_x
    ))
    oscar_args <- utils::modifyList(naive_args, list(mu_x_a_coefs = outcome_model_obs$mu_x_a_coefs))
    if(sampleSplit) {
      ss_r_oscar_args <- utils::modifyList(naive_args, list(
        K = 5,
        outcome_model_obs_coefs = outcome_model_obs$mu_x_a_coefs,
        m_of_x_mat = calibrated_outcome_model_rct$mu_x
      ))
    }
  }
  
  estimators <- list(
    naive = naive_args,
    racer = racer_args
  )
  if(randomForest) {
    estimators[['naive_RF']] <- c(naive_args, useRF = TRUE)
    racer_rf <- c(racer_args, useRF = TRUE)
    racer_rf[['m_of_x_mat']] <- muXrf
    estimators[['racer_RF']] <- racer_rf
  }
  if(!is.null(os)) {
    estimators[['oscar']] <- oscar_args
    estimators[['r_oscar']] <- r_oscar_args
    if(sampleSplit) {
      estimators[['ss_r_oscar']] <- ss_r_oscar_args
    }
  }

  o <- lapply(estimators, function(i) {
    do.call(estimate_cate, i)
  })
  class(o) <- 'cate_model'
  o
}
