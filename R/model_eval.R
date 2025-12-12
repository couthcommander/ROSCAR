#' Model Evaluation
#'
#' Compare model estimates
#'
#' @param mod \sQuote{cate_model} object, created with \code{\link{cate_model}}.
#' @param tau NEEDINFO
#'
#' @return
#' \item{rmse}{NEEDINFO}
#' \item{rank_corr}{NEEDINFO}
#' \item{accuracy_itr}{NEEDINFO}
#' \item{cal_intercept}{NEEDINFO}
#' \item{cal_slope}{NEEDINFO}
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
