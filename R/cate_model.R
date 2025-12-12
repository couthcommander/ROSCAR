#' Compute Model Coefficients
#'
#' CATE estimates
#'
#' @param rct NEEDINFO
#' @param os NEEDINFO
#' @param randomForest NEEDINFO
#' @param sampleSplit NEEDINFO
#'
#' @return
#' \item{estimates by method}{
#' \describe{
#' \item{cate_coefficients}{cate estimate; not available in
#' naive/RACER random forest approach}
#' \item{cate_function}{function to predict CATE for new data}
#' \item{fitted_cate}{fitted CATE values}
#' \item{delta_coef_mat}{difference between coefficient estimates;
#' only available in OSCAR approach}
#' \item{cate_delta_coefficients}{difference between the two CATEs;
#' only available in R-OSCAR approach}
#' \item{cate_model}{CATE model; only available in naive random
#' forest approach}
#' }}
#'
#' @examples
#' data(os, rct)
#' modNoOs <- cate_model(rct)
#' print(modNoOs, max = 40)
#' modOs <- cate_model(rct, os)
#' print(modOs, max = 40)
#'
#' @export

cate_model <-
function(rct, os = NULL,
          randomForest = FALSE,
          sampleSplit = FALSE,
          method = c('naive', 'racer', 'oscar', 'roscar')) {
  args <- list(rct = rct, os = os)
  if(is.null(os)) {
    # remove OSCAR/R-OSCAR if no "OS" data
    method <- setdiff(method, c('oscar', 'roscar'))
  }
  ll <- lapply(method, function(m) switch(m,
    naive = do.call(model_naive, c(args, randomForest = randomForest)),
    racer = do.call(model_racer, c(args, randomForest = randomForest)),
    oscar = do.call(model_oscar, args),
    roscar = do.call(model_roscar, c(args, sampleSplit = sampleSplit)),
    NULL
  ))
  o <- Reduce(modifyList, ll)
  class(o) <- 'cate_model'
  o
}

#' @keywords internal
cate_coefficients <- function(m) {
  if(!inherits(m, 'cate_model')) {
    stop("'m' should inherit class 'cate_model'")
  }
  cc <- lapply(m, function(i) unname(i[['cate_coefficients']][,1]))
  # random forest models don't include coefficients
  as.data.frame(cc[lengths(cc) > 0])
}

#' @keywords internal
cate_function <- function(m) {
  if(!inherits(m, 'cate_model')) {
    stop("'m' should inherit class 'cate_model'")
  }
  lapply(m, function(i) i$cate_function)
}

#' @keywords internal
fitted_cate <- function(m) {
  if(!inherits(m, 'cate_model')) {
    stop("'m' should inherit class 'cate_model'")
  }
  as.data.frame(lapply(m, function(i) i$fitted_cate))
}
