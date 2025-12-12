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
#' mod <- cate_model(rct, os)
#' print(mod, max = 40)
#'
#' @rdname print
#' @export

print.cate_model <-
function(x, ...) {
  d <- cate_coefficients(x)
  print(d, ...)
  invisible(x)
}
