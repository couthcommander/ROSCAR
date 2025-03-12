#' estimate_cate
#'
#' Estimate conditional average treatment effect with naive/RACER approach.
#'
#' @keywords internal
#' @param X NEEDINFO
#' @param A NEEDINFO
#' @param Y NEEDINFO
#' @param propensity_vec NEEDINFO
#' @param m_of_x_mat NEEDINFO
#' @param normalize_weights NEEDINFO

estimate_cate_with_rct <-
function(X, A, Y, propensity_vec, 
                                m_of_x_mat = NULL, ## regression functions for each trt level
                                normalize_weights = TRUE)
{
  stopifnot("A must be either 0 or 1" = all(sort(unique(A)) == c(0,1)))

  ## inverse propensity score weights
  propens_weights <- ifelse(A == 1, 1/propensity_vec, 1/(1-propensity_vec))

  ## turn binary A into +1 and -1 variable
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
    ## compute shift function
    m_of_x_avg <- rowSums(m_of_x_mat * cbind(propensity_vec, 1 - propensity_vec))
  }

  ## construct working response
  y_tilde <- (A_1_m1 * (Y - m_of_x_avg) * propens_weights) * mean_propens

  ## estimate linear CATE by predicting working response with covariates
  cate_model <- glmnet::cv.glmnet(x = X, y = y_tilde)

  ## estimated coefficients for the linear estimate of the CATE
  cate_coefficients <- stats::predict(cate_model, type = "coef", s = "lambda.min")

  ## function to predict CATE for new data
  return_cate <- function(x)
  {
    unname(drop(stats::predict(cate_model, newx = x, type = "response", s = "lambda.min")))
  }

  list(
    cate_coefficients = cate_coefficients, 
    cate_function = return_cate,
    fitted_cate = return_cate(X)
  )
}
