#' model_oscar
#'
#' Estimate conditional average treatment effect.
#'
#' Users will likely not call this function directly, instead
#' relying in behavior in \code{\link{cate_model}}.
#'
#' @export

model_oscar <-
function(rct, os) {
    rctX <- rct$X
    rctY <- rct$Y
    rctA <- rct$A
    propensity_vec <- rep(mean(rctA), length(rctY))
    stopifnot("A must be either 0 or 1" = all(sort(unique(rctA)) == c(0,1)))
    stopifnot("A must be either 0 or 1" = all(sort(unique(os$A)) == c(0,1)))

#####
# estimating the OS outcome mean
#
# \hat{\gamma}_a^o =
#                    \argmin_\gamma \frac{1}{n_a^o}\sum_{i: A_i^o = a}
#                        \left[ Y_i^o - \gamma^T X^o_i \right]^2 +
#                        \lambda^o_a\left\|\gamma\right\|_{1}
#####
    trt_idx_OS <- tapply(seq_along(os$A), os$A, I)
    ## outcome regression coefficient matrix -- one column for each treatment
    mu_x_a_coefs <- vapply(trt_idx_OS, function(i) {
        lasso_outcome_mod <- glmnet::cv.glmnet(x = os$X[i,,drop = FALSE], y = os$Y[i])
        unname(drop(as.vector(stats::predict(lasso_outcome_mod, type = "coef", s = "lambda.min"))))
    }, numeric(ncol(os$X)+1))
    ## estimated main effects from obs data
    m_of_x_mat <- cbind(1, rctX) %*% mu_x_a_coefs

    propens_ratio <- ifelse(rctA == 1, (1-propensity_vec)/propensity_vec, propensity_vec/(1-propensity_vec))
    ## compute shift function (CMO)
    m_of_x_avg <- (propens_ratio + 1) * (rctA * m_of_x_mat[,2] +
                                        (1 - rctA) * m_of_x_mat[,1])
    ## inverse propensity score weights
    propens_weights <- ifelse(rctA == 1, 1/propensity_vec, 1/(1-propensity_vec))
    ## construct working response
    y_tilde <- rctY * propens_weights - m_of_x_avg

#####
# instantiate to learn CMO and CATE jointly through the estimation
# of the sparse mean discrepancy cofficient vectors
#
# \hdelta_{a}^j =
#                 \argmin_{d_a} \frac{1}{n^r_a}\sum_{i: A_i^r = a} &
#                     \left[ \left(
#     \frac{a}{\pi^r_{a}} Y_i^r - a\left(1 + \alpha^a\right) \hat{\gamma}_a^{o^T} X^r_i\right) -
#         a \left(1 + \alpha^a\right) d_{a}^T X^r_i \right]^2 +
#                     \lambda_a^j \|d_{a}\|_{1}
#####
    nc <- ncol(rctX)
    trt_idx <- tapply(seq_along(rctA), rctA, I)
    delta_coef_mat <- vapply(trt_idx, function(i) {
        ## need to force an intercept b/c need to mult by factor
        lasso_outcome_mod <- glmnet::cv.glmnet(x = (1 + propens_ratio[i]) * cbind(1, rctX[i,,drop = FALSE]),
                                    y = y_tilde[i],
                                    intercept = FALSE,
                                    ## dont penalize main effect of trtment
                                    penalty.factor = c(0, rep(1, nc))
        )
        unname(drop(as.vector(stats::predict(lasso_outcome_mod, type = "coef", s = "lambda.min"))))[-1]
    }, numeric(nc+1))
    ## estimate cate
    overall_regression_coefs <- mu_x_a_coefs + delta_coef_mat
    ## estimated coefficients for the linear estimate of the CATE
    cate_coefficients <- matrix(overall_regression_coefs[,2] - overall_regression_coefs[,1], ncol = 1)

    ## function to predict CATE for new data
    return_cate <- function(x)
    {
        unname(drop(cbind(1, x) %*% cate_coefficients))
    }

    o <- list(list(
        cate_coefficients = cate_coefficients,
        cate_function = return_cate,
        fitted_cate = return_cate(rctX),
        delta_coef_mat = delta_coef_mat,
        cate_delta_coefficients = NULL,
        cate_model = NULL
    ))
    names(o) <- 'oscar'
    class(o) <- 'cate_model'
    o
}
