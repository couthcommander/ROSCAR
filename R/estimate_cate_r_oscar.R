#' estimate_cate
#'
#' Estimate conditional average treatment effect with R-OSCAR method.
#'
#' @keywords internal
#' @param X NEEDINFO
#' @param A NEEDINFO
#' @param Y NEEDINFO
#' @param propensity_vec NEEDINFO
#' @param init_cate_coefs NEEDINFO
#' @param m_of_x_mat NEEDINFO
#' @param normalize_weights NEEDINFO

estimate_cate_r_oscar <-
function(X, A, Y, propensity_vec, 
                                   init_cate_coefs,
                                   m_of_x_mat = NULL,
                                   normalize_weights = TRUE)
{
  stopifnot("A must be either 0 or 1" = all(sort(unique(A)) == c(0,1)))

  propens_weights <- ifelse(A == 1, 1/propensity_vec, 1/(1-propensity_vec))

  A_1_m1 <- (2 * A - 1)

  if (normalize_weights)
  {
    mean_propens <- mean(propens_weights)
    propens_weights <- propens_weights / mean_propens
  } else
  {
    mean_propens <- 1
  }

  if (is.null(m_of_x_mat))
  {
    m_of_x_avg <- numeric(NROW(X))
  } else
  {
    m_of_x_avg <- rowSums(m_of_x_mat * cbind(propensity_vec, 1 - propensity_vec))
  }

  ## initial/preliminary cate estimate for the training data
  init_cate_est <- drop(cbind(1, X) %*% init_cate_coefs)

  ## construct working response...
  ## here we subtract the preliminary estimate of the CATE and fit the residual
  y_tilde <- (A_1_m1 * (Y - m_of_x_avg) * propens_weights) * mean_propens - init_cate_est

  ## fit model for the DIFFERENCE between the preliminary CATE estimate
  ## and the final CATE estimate. The final CATE estimate is then the prelim est + the difference (delta)
  cate_model <- glmnet::cv.glmnet(x = X, y = y_tilde)

  ## coefficients for the difference between the two CATES
  cate_delta_coefficients <- stats::predict(cate_model, type = "coef", s = "lambda.min")

  ## final estimate of the CATE...
  ## estimated coefficients for the linear estimate of the CATE
  cate_coefficients <- as.matrix(cate_delta_coefficients + init_cate_coefs)

  ## fitted CATE values for training data
  fitted_cate <- stats::predict(cate_model, newx = X, type = "response", s = "lambda.min") + init_cate_est
  fitted_cate <- unname(drop(fitted_cate))

  ## function to predict CATE for new data
  return_cate <- function(x)
  {
    unname(drop(cbind(1, x) %*% cate_coefficients))
  }

  list(
    cate_coefficients = cate_coefficients, 
    cate_delta_coefficients = cate_delta_coefficients, 
    cate_function = return_cate,
    fitted_cate = fitted_cate
  )
}
