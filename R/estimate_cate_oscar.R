#' estimate_cate
#'
#' Estimate conditional average treatment effect with R-OSCAR method under sample-split.
#'
#' @keywords internal
#' @param X NEEDINFO
#' @param A NEEDINFO
#' @param Y NEEDINFO
#' @param propensity_vec NEEDINFO
#' @param mu_x_a_coefs NEEDINFO

estimate_cate_oscar <-
function(X, A, Y, propensity_vec, 
                              mu_x_a_coefs) ## regression function coefficients for each trt level
{
  
  stopifnot("A must be either 0 or 1" = all(sort(unique(A)) == c(0,1)))
  
  propens_ratio <- ifelse(A == 1, (1-propensity_vec)/propensity_vec, propensity_vec/(1-propensity_vec))

  propens_weights <- ifelse(A == 1, 1/propensity_vec, 1/(1-propensity_vec))
  
  normalize_weights <- FALSE
  
  if (normalize_weights)
  {
    mean_propens <- mean(propens_weights)
    propens_weights <- propens_weights / mean_propens
  } else
  {
    mean_propens <- 1
  }

  if (is.null(mu_x_a_coefs))
  {
    m_of_x_avg <- numeric(NROW(X))
    mu_x_a_coefs <- 0
  } else
  {
    ## estimated main effects from obs data
    m_of_x_mat <- cbind(1, X) %*% mu_x_a_coefs
    ## compute shift function
    m_of_x_avg <- (propens_ratio + 1) *  (A * m_of_x_mat[,2] + 
                                            (1 - A) * m_of_x_mat[,1])
  }
  y_tilde <- (Y * propens_weights - m_of_x_avg) * mean_propens
  
  trt_idx <- tapply(seq_along(A), A, I)
  delta_coef_mat <- vapply(trt_idx, function(i) {
    lasso_outcome_mod <- glmnet::cv.glmnet(x = (1 + propens_ratio[i]) * cbind(1, X[i,,drop = FALSE]), ## need to force an intercept b/c need to mult by factor
                                   y = y_tilde[i],
                                   intercept = FALSE,
                                   penalty.factor = c(0, rep(1, NCOL(X))) ## dont penalize main effect of trtment
    )
    unname(drop(as.vector(stats::predict(lasso_outcome_mod, type = "coef", s = "lambda.min"))))[-1]
  }, numeric(ncol(X)+1))
  
  # estimate cate
  overall_regression_coefs <- mu_x_a_coefs + delta_coef_mat
  
  ## estimated coefficients for the linear estimate of the CATE
  cate_coefficients <- matrix(overall_regression_coefs[,2] - overall_regression_coefs[,1], ncol = 1)
  
  ## function to predict CATE for new data
  return_cate <- function(x)
  {
    unname(drop(cbind(1, x) %*% cate_coefficients))
  }
  
  list(cate_coefficients = cate_coefficients,
       cate_function = return_cate,
       fitted_cate = return_cate(X),
       delta_coef_mat = delta_coef_mat)
}
