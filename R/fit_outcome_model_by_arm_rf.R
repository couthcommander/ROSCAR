#' Fit Outcome Model with Random Forest
#'
#' Fit outcome regression model for each treatment arm.
#'
#' @keywords internal
#' @param X Data set covariates
#' @param A Vector with treatment assignment
#' @param Y Vector with outcome variable

fit_outcome_model_by_arm_rf <-
function(X, A, Y)
{
  ## Convert X to data frame with column names
  X <- as.data.frame(X)
  names(X) <- paste0("Var", seq_len(ncol(X)))

  trt_idx <- tapply(seq_along(A), A, I)
  ## list to store trained models -- one for each treatment
  rf_models <- lapply(trt_idx, function(i) {
    ## train random forest model
    train_rf(Y[i], X[i,,drop = FALSE])
  })

  ## function to predict mu_x for new data
  return_regression_mat <- function(x)
  {
    vapply(rf_models, stats::predict, numeric(nrow(x)), newdata = x)
  }
  
  ## predict mu_x for the training data
  mu_x_mat <- return_regression_mat(X)
  
  list(mu_x = mu_x_mat, 
       outcome_models = rf_models, 
       outcome_function = return_regression_mat)
}
