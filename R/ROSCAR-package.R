#' Package
#'
#' R-OSCAR (Robust Observational Studies for CMO-Augmented RCT)
#' is a two-stage estimator that calibrates OS outcome predictions
#' to the RCT population and corrects residual biases through
#' regularized regression. Comparison methods are also provided.
#'
#' @docType package
#'
#' @author Amir Asiaeetaheri \email{amir.asiaeetaheri@vumc.org}
#'
#' Maintainer: Amir Asiaeetaheri \email{amir.asiaeetaheri@vumc.org}
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom methods formalArgs
#' @importFrom stats coef cor lm predict rbinom rnorm runif
#' @importFrom utils modifyList
#'
#' @examples
#' data(os, rct, tau)
#' mod <- cate_model(rct, os, randomForest = TRUE)
#' head(model_pred(mod))
#' head(model_pred(mod, rct$X))
#' model_eval(mod, tau)

"_PACKAGE"
