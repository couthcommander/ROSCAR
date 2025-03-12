#' Calibrate Outcome Model
#'
#' Calibrate the outcome regression functions to (RCT) data.
#'
#' @keywords internal
#' @param X NEEDINFO
#' @param A NEEDINFO
#' @param Y NEEDINFO
#' @param mu_x_a_coefs NEEDINFO

calibrate_outcome_model <-
function(X, A, Y, mu_x_a_coefs)
{
  trt_idx <- tapply(seq_along(A), A, I)
  ## delta coefficient matrix
  delta_coef_mat <- vapply(seq_along(trt_idx), function(l) {
    i <- trt_idx[[l]]
    ## get preliminary prediction of outcome function based on prior estimate
    pred_cur <- drop(cbind(1, X[i,,drop = FALSE]) %*% mu_x_a_coefs[,l])
    lasso_outcome_mod <- glmnet::cv.glmnet(x = X[i,,drop = FALSE], y = Y[i], offset = pred_cur)
    unname(drop(as.vector(stats::predict(lasso_outcome_mod, type = "coef", s = "lambda.min"))))
  }, numeric(ncol(X)+1))
  colnames(delta_coef_mat) <- names(trt_idx)
  delta_coefs_current <- delta_coef_mat[,ncol(delta_coef_mat)]

  ## overall outcome regression coefficients
  cal_regression_coef_mat <- mu_x_a_coefs + delta_coef_mat

  return_regression_mat <- function(x)
  {
    mu_x_mat <- cbind(1, x) %*% cal_regression_coef_mat
    colnames(mu_x_mat) <- colnames(cal_regression_coef_mat)
    mu_x_mat
  }

  list(mu_x = return_regression_mat(X),
       mu_x_a_coefs = cal_regression_coef_mat,
       delta_coefs = delta_coefs_current,
       outcome_function = return_regression_mat)
}
