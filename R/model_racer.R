#' model_racer
#'
#' Estimate conditional average treatment effect.
#'
#' Users will likely not call this function directly, instead
#' relying in behavior in \code{\link{cate_model}}.
#'
#' @export

model_racer <-
function(rct, os = NULL, randomForest = FALSE) {
    rctX <- rct$X
    rctY <- rct$Y
    rctA <- rct$A
    propensity_vec <- rep(mean(rctA), length(rctY))
    stopifnot("A must be either 0 or 1" = all(sort(unique(rctA)) == c(0,1)))

    trt_idx <- tapply(seq_along(rctA), rctA, I)
    if(randomForest) {
        m_of_x_mat <- rf_racer_mu_x(rctX, rctY, trt_idx)
    } else {
        ## outcome regression coefficient matrix -- one column for each treatment
        regression_coef_mat <- vapply(trt_idx, function(i) {
            lasso_outcome_mod <- glmnet::cv.glmnet(x = rctX[i,,drop = FALSE], y = rctY[i])
            unname(drop(as.vector(stats::predict(lasso_outcome_mod, type = "coef", s = "lambda.min"))))
        }, numeric(ncol(rctX)+1))
        m_of_x_mat <- cbind(1, rctX) %*% regression_coef_mat
        colnames(m_of_x_mat) <- colnames(regression_coef_mat)
    }

    ## inverse propensity score weights
    propens_weights <- ifelse(rctA == 1, 1/propensity_vec, 1/(1-propensity_vec))
    ## turn binary A into +1 and -1 variable
    A_1_m1 <- (2 * rctA - 1)
    ## normalize weights
    mean_propens <- mean(propens_weights)
    propens_weights <- propens_weights / mean_propens
    ## compute shift function (CMO)
    m_of_x_avg <- rowSums(m_of_x_mat * cbind(propensity_vec, 1 - propensity_vec))
    ## construct working response
    y_tilde <- (A_1_m1 * (rctY - m_of_x_avg) * propens_weights) * mean_propens

    if(randomForest) {
        ## estimate CATE by predicting working response with covariates using random forest
        cate_model <- caret::train(y = y_tilde, x = X, method = "rf", trControl = train_control,
                            tuneGrid = tGrid, ntree = 50, nodesize = 25)
        cate_coefficients <- NULL
        ## function to predict CATE for new data
        return_cate <- function(x)
        {
          reqNames <- cate_model[['finalModel']][['xNames']]
          if(inherits(x, 'data.frame')) {
            names(x) <- reqNames
          } else if(inherits(x, 'matrix')) {
            colnames(x) <- reqNames
          }
          unname(drop(stats::predict(cate_model, newdata = x)))
        }
        label <- 'racer_RF'
    } else {
        ## estimate linear CATE by predicting working response with covariates
        cate_model <- glmnet::cv.glmnet(x = rctX, y = y_tilde)
        ## estimated coefficients for the linear estimate of the CATE
        cate_coefficients <- stats::predict(cate_model, type = "coef", s = "lambda.min")
        if(sum(abs(as.matrix(cate_coefficients)[-1,]) > 1e-8) == 0) {
            warning('a constant treatment effect estimate was returned')
        }
        ## function to predict CATE for new data
        return_cate <- function(x)
        {
            unname(drop(stats::predict(cate_model, newx = x, type = "response", s = "lambda.min")))
        }
        label <- 'racer'
    }

    o <- list(list(
        cate_coefficients = cate_coefficients,
        cate_function = return_cate,
        fitted_cate = return_cate(rctX),
        delta_coef_mat = NULL,
        cate_delta_coefficients = NULL,
        cate_model = cate_model
    ))
    names(o) <- label
    class(o) <- 'cate_model'
    o
}

#' @keywords internal
rf_racer_mu_x <- function(rctX, rctY, trt_idx) {
    ## Convert X to data frame with column names
    X <- as.data.frame(rctX)
    names(X) <- paste0("Var", seq_len(ncol(X)))
    tGrid <- data.frame(mtry = floor(sqrt(ncol(X))))
    ## define the training control
    train_control <- caret::trainControl(method = "cv", number = 5)

    ## list to store trained models -- one for each treatment
    rf_models <- lapply(trt_idx, function(i) {
        ## train random forest model
        caret::train(
            y = rctY[i],
            x = X[i,,drop = FALSE],
            method = 'rf',
            trControl = train_control,
            tuneGrid = tGrid,
            ntree = 50,
            nodesize = 25
        )
    })
    ## predict mu_x for the training data
    vapply(rf_models, stats::predict, numeric(nrow(X)), newdata = X)
}
