#' Model Evaluation
#'
#' Compare model estimates
#'
#' @param mod \sQuote{cate_model} object, created with \code{\link{cate_model}}.
#' @param tau True treatment effect for RCT data.
#'
#' @return
#' \item{rmse}{Root mean squared error}
#' \item{rank_corr}{Spearman correlation}
#' \item{accuracy_itr}{Model accuracy with respect to sign}
#' \item{cal_intercept}{Calibration intercept}
#' \item{cal_slope}{Calibration coefficients}
#'
#' @examples
#' data(os, rct, tau)
#' mod <- cate_model(rct, os)
#' model_eval(mod, tau)
#'
#' @export

model_eval <-
function(mod, tau) {
  if(inherits(mod, 'cate_model')) {
    est <- model_pred(mod)
  } else if(inherits(mod, 'cate_est')) {
    est <- mod
  } else {
    stop("'mod' should inherit class 'cate_model'")
  }
  spearman_corr <- function(t, e) {
    tryCatch(stats::cor(t, e, method = "spearman"), warning = function(w) {
      warning(w)
      if(grepl('the standard deviation is zero', w)) {
        return(0)
      }
    })
  }
  t(vapply(est, function(e) {
    cal_coefs <- unname(stats::coef(stats::lm(tau ~ e)))
    c(
      rmse          = sqrt(mean((tau - e)^2)),
      rank_corr     = spearman_corr(tau, e),
      accuracy_itr  = mean(sign(tau) == sign(e)),
      cal_intercept = cal_coefs[1],
      cal_slope     = cal_coefs[2]
    )
  }, numeric(5)))
}
