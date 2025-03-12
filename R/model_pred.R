#' Model Predict
#'
#' CATE estimates with new data
#'
#' @param object \sQuote{cate_model} object, created with \code{\link{model_coef}}.
#' @param newdata NEEDINFO
#'
#' @return A data.frame with estimates from each method.
#'
#' @examples
#' data(os, rct)
#' mod <- model_coef(rct, os)
#' prd <- model_pred(mod)
#' prdNewdata <- model_pred(mod, rct$X)
#'
#' @export

model_pred <-
function(object, newdata) {
  if(!inherits(object, 'cate_model')) {
    stop("'object' should inherit class 'cate_model'")
  }
  if(!missing(newdata)) {
    ests <- lapply(object, function(i) {
      unname(i$cate_function(newdata))
    })
  } else {
    ests <- lapply(object, function(i) {
      unname(i$fitted_cate)
    })
  }
  o <- as.data.frame(ests)
  class(o) <- c('data.frame', 'cate_est')
  o
}
