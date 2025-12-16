#' Model Predict
#'
#' CATE estimates with new data
#'
#' @param object \sQuote{cate_model} object, created with \code{\link{cate_model}}.
#' @param newdata An optional data frame in which to look for
#' variables with which to predict. If omitted, the fitted values are used.
#'
#' @return A data.frame with estimates from each method.
#'
#' @examples
#' data(os, rct)
#' mod <- cate_model(rct, os)
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
    ests <- lapply(cate_function(object), function(f) unname(f(newdata)))
    o <- as.data.frame(ests)
  } else {
    o <- fitted_cate(object)
  }
  class(o) <- c('data.frame', 'cate_est')
  o
}
