#' Print methods
#'
#' Print method for \sQuote{cate_model} objects.
#'
#' Print its argument and return it invisibly.
#'
#' @param x A cate_model object.
#' @param \dots Additional parameters passed to print.data.frame.
#'
#' @examples
#' data(os, rct)
#' mod <- model_coef(rct, os)
#' print(mod, max = 40)
#'
#' @rdname print
#' @export

print.cate_model <-
function(x, ...) {
  cc <- lapply(x, function(i) i[['cate_coefficients']][,1])
  # random forest models don't include coefficients
  d <- as.data.frame(cc[lengths(cc) > 0])
  print(d, ...)
  invisible(x)
}
