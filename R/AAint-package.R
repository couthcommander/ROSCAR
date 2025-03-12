#' Package
#'
#' what does it do?
#'
#' @docType package
#'
#' @author Amir Asiaeetaheri \email{amir.asiaeetaheri@vumc.org}
#'
#' Maintainer: Amir Asiaeetaheri \email{amir.asiaeetaheri@vumc.org}
#'
#' @importFrom caret createFolds train trainControl
#' @importFrom glmnet cv.glmnet
#' @importFrom stats coef cor lm predict
#' @importFrom utils modifyList
#'
#' @examples
#' data(os, rct, tau)
#' mod <- model_coef(rct, os, randomForest = TRUE)
#' head(model_pred(mod))
#' head(model_pred(mod, rct$X))
#' model_eval(mod, tau)

"_PACKAGE"
