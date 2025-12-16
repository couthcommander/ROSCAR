#' Fit Outcome Model
#'
#' Fit outcome regression model for each treatment arm.
#'
#' @keywords internal
#' @param X Data set covariates
#' @param A Vector with treatment assignment
#' @param Y Vector with outcome variable

fit_outcome_model_by_arm <-
function(X, A, Y)
{
  trt_idx <- tapply(seq_along(A), A, I)
  ## outcome regression coefficient matrix -- one column for each treatment
  regression_coef_mat <- vapply(trt_idx, function(i) {
    lasso_outcome_mod <- glmnet::cv.glmnet(x = X[i,,drop = FALSE], y = Y[i])
    unname(drop(as.vector(stats::predict(lasso_outcome_mod, type = "coef", s = "lambda.min"))))
  }, numeric(ncol(X)+1))

  return_regression_mat <- function(x)
  {
    mu_x_mat <- cbind(1, x) %*% regression_coef_mat
    colnames(mu_x_mat) <- colnames(regression_coef_mat)
    mu_x_mat
  }
  
  list(mu_x = return_regression_mat(X),
       mu_x_a_coefs = regression_coef_mat, 
       outcome_function = return_regression_mat)
}
