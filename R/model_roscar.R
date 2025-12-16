#' @rdname cate_model
#' @export

model_roscar <-
function(rct, os, sampleSplit = FALSE) {
    rctX <- rct$X
    rctY <- rct$Y
    rctA <- rct$A
    propensity_vec <- rep(mean(rctA), length(rctY))
    stopifnot("A must be either 0 or 1" = all(sort(unique(rctA)) == c(0,1)))
    stopifnot("A must be either 0 or 1" = all(sort(unique(os$A)) == c(0,1)))

    # estimating the OS outcome mean
    trt_idx_OS <- tapply(seq_along(os$A), os$A, I)
    ## outcome regression coefficient matrix -- one column for each treatment
    mu_x_a_coefs <- vapply(trt_idx_OS, function(i) {
        lasso_outcome_mod <- glmnet::cv.glmnet(x = os$X[i,,drop = FALSE], y = os$Y[i])
        unname(drop(as.vector(stats::predict(lasso_outcome_mod, type = "coef", s = "lambda.min"))))
    }, numeric(ncol(os$X)+1))

    if(sampleSplit) {
        return(ss_roscar(rct, propensity_vec, mu_x_a_coefs))
    }

    # calibrate the OS outcome means to the RCT outcomes by estimating
    # a sparse mean discrepancy vector for each treatment arm
    ## calibrate the above regression functions to the RCT data
    trt_idx <- tapply(seq_along(rctA), rctA, I)
    ## delta coefficient matrix
    delta_coef_mat <- vapply(seq_along(trt_idx), function(l) {
        i <- trt_idx[[l]]
        ## get preliminary prediction of outcome function based on prior estimate
        pred_cur <- drop(cbind(1, rctX[i,,drop = FALSE]) %*% mu_x_a_coefs[,l])
        lasso_outcome_mod <- glmnet::cv.glmnet(x = rctX[i,,drop = FALSE], y = rctY[i], offset = pred_cur)
        unname(drop(as.vector(stats::predict(lasso_outcome_mod, type = "coef", s = "lambda.min"))))
    }, numeric(ncol(rctX)+1))
    colnames(delta_coef_mat) <- names(trt_idx)

    # this produces a CMO and a preliminary CATE estimate tailored for the RCT
    ## overall outcome regression coefficients
    overall_regression_coefs <- mu_x_a_coefs + delta_coef_mat
    ## estimated coefficients
    init_cate_coefs <- overall_regression_coefs[,2] - overall_regression_coefs[,1]
    m_of_x_mat <- cbind(1, rctX) %*% overall_regression_coefs
    colnames(m_of_x_mat) <- colnames(overall_regression_coefs)

    ## inverse propensity score weights
    propens_weights <- ifelse(rctA == 1, 1/propensity_vec, 1/(1-propensity_vec))
    ## turn binary A into +1 and -1 variable
    A_1_m1 <- (2 * rctA - 1)
    mean_propens <- mean(propens_weights)
    propens_weights <- propens_weights / mean_propens
    ## compute shift function (CMO)
    m_of_x_avg <- rowSums(m_of_x_mat * cbind(propensity_vec, 1 - propensity_vec))

    ## initial/preliminary cate estimate for the training data
    init_cate_est <- drop(cbind(1, rctX) %*% init_cate_coefs)
    ## construct working response...
    ## here we subtract the preliminary estimate of the CATE and fit the residual
    y_tilde <- (A_1_m1 * (rctY - m_of_x_avg) * propens_weights) * mean_propens - init_cate_est

    ## fit model for the DIFFERENCE between the preliminary CATE estimate
    ## and the final CATE estimate.
    cate_model <- glmnet::cv.glmnet(x = rctX, y = y_tilde)

    ## coefficients for the difference between the two CATES
    cate_delta_coefficients <- stats::predict(cate_model, type = "coef", s = "lambda.min")

    # refine the preliminary CATE by estimating a sparse linear discrepancy function to correct
    # any residual misspecification using the calibrated CMO as the augmentation function
    ## The final CATE estimate is then the prelim est + the difference (delta)
    cate_coefficients <- as.matrix(cate_delta_coefficients + init_cate_coefs)

    ## fitted CATE values for training data
    fitted_cate <- stats::predict(cate_model, newx = rctX, type = "response", s = "lambda.min") + init_cate_est
    fitted_cate <- unname(drop(fitted_cate))

    ## function to predict CATE for new data
    return_cate <- function(x)
    {
        unname(drop(cbind(1, x) %*% cate_coefficients))
    }

    o <- list(list(
        cate_coefficients = cate_coefficients,
        cate_function = return_cate,
        fitted_cate = fitted_cate,
        delta_coef_mat = NULL,
        cate_delta_coefficients = cate_delta_coefficients,
        cate_model = NULL
    ))
    names(o) <- 'r_oscar'
    class(o) <- 'cate_model'
    o
}

#' @keywords internal
ss_roscar <- function(rct, propensity_vec, mu_x_a_coefs) {
    rctX <- rct$X
    rctY <- rct$Y
    rctA <- rct$A
    K <- 5 # make this argument?
    s_cates_coeffs <- matrix(NA, nrow = ncol(rctX)+1, ncol = K) #CATEs will be the columns
    if(!requireNamespace('caret', quietly = TRUE)) {
        stop('package "caret" is required for sample split in R-OSCAR')
    }
    # To have a balance of treatment in each split use AR for creating folds.
    index <- caret::createFolds(rctA,
                        k = K,
                        list = TRUE,
                        returnTrain = FALSE)
    for (i in seq_len(K)) {
        idx_cur <- index[[i]]
        X_cur <- rctX[idx_cur,,drop = FALSE]
        X_i <- rctX[-idx_cur,,drop = FALSE]
        A_i <- rctA[-idx_cur]
        Y_i <- rctY[-idx_cur]
        trt_idx <- tapply(seq_along(A_i), A_i, I)
        ## delta coefficient matrix
        dcm_i <- vapply(seq_along(trt_idx), function(l) {
            i <- trt_idx[[l]]
            pred_cur <- drop(cbind(1, X_i[i,,drop = FALSE]) %*% mu_x_a_coefs[,l])
            lasso_outcome_mod <- glmnet::cv.glmnet(x = X_i[i,,drop = FALSE], y = Y_i[i], offset = pred_cur)
            unname(drop(as.vector(stats::predict(lasso_outcome_mod, type = "coef", s = "lambda.min"))))
        }, numeric(ncol(X_i)+1))
        colnames(dcm_i) <- names(trt_idx)
        ## outcome regression coefficients
        cal_mu_x_a_coefs <- mu_x_a_coefs + dcm_i
        ## estimated coefficients
        cal_cate_coefs <- cal_mu_x_a_coefs[,2] - cal_mu_x_a_coefs[,1]
        cal_m_of_x_mat <- cbind(1, X_cur) %*% cal_mu_x_a_coefs
        colnames(cal_m_of_x_mat) <- colnames(cal_mu_x_a_coefs)

        pv_cur <- propensity_vec[idx_cur]
        A_cur <- rctA[idx_cur]
        Y_cur <- rctY[idx_cur]
        propens_weights <- ifelse(A_cur == 1, 1/pv_cur, 1/(1-pv_cur))
        A_1_m1 <- (2 * A_cur - 1)
        mean_propens <- mean(propens_weights)
        propens_weights <- propens_weights / mean_propens
        m_of_x_avg <- rowSums(cal_m_of_x_mat * cbind(pv_cur, 1 - pv_cur))
        init_cate_est <- drop(cbind(1, X_cur) %*% cal_cate_coefs)
        y_tilde <- (A_1_m1 * (Y_cur - m_of_x_avg) * propens_weights) * mean_propens - init_cate_est
        cate_model <- glmnet::cv.glmnet(x = X_cur, y = y_tilde)
        cate_delta_coefficients <- stats::predict(cate_model, type = "coef", s = "lambda.min")
        split_cc <- as.matrix(cate_delta_coefficients + cal_cate_coefs)
        s_cates_coeffs[, i] <- as.numeric(split_cc)
    }
    cate_coefficients <- as.matrix(rowMeans(s_cates_coeffs))
    ## function to predict CATE for new data
    return_cate <- function(x)
    {
        unname(drop(cbind(1, x) %*% cate_coefficients))
    }
    o <- list(list(
        cate_coefficients = cate_coefficients,
        cate_function = return_cate,
        fitted_cate = return_cate(rctX),
        delta_coef_mat = NULL,
        cate_delta_coefficients = NULL,
        cate_model = NULL
    ))
    names(o) <- 'ss_r_oscar'
    class(o) <- 'cate_model'
    o
}
