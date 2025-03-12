#' estimate_cate
#'
#' Estimate conditional average treatment effect.
#'
#' Users will likely not call this function directly, instead
#' relying in behavior in \code{\link{model_coef}}.
#'
#' @param X NEEDINFO
#' @param A NEEDINFO
#' @param Y NEEDINFO
#' @param propensity_vec NEEDINFO -- is this always `mean(A)`?
#' @param m_of_x_mat NEEDINFO
#' @param normalize_weights NEEDINFO
#' @param mu_x_a_coefs NEEDINFO (required for OSCAR approach)
#' @param init_cate_coefs NEEDINFO (required for R-OSCAR approach)
#' @param K NEEDINFO (required for sample split R-OSCAR approach)
#' @param outcome_model_obs_coefs NEEDINFO (required for sample split R-OSCAR approach)
#' @param useRF NEEDINFO (required for naive random forest approach)
#'
#' @return
#' \item{cate_coefficients}{cate estimate; not available in
#' naive random forest approach}
#' \item{cate_function}{function to predict CATE for new data}
#' \item{fitted_cate}{fitted CATE values}
#' \item{delta_coef_mat}{difference between coefficient estimates;
#' only available in OSCAR approach}
#' \item{cate_delta_coefficients}{difference between the two CATEs;
#' only available in R-OSCAR approach}
#' \item{cate_model}{CATE model; only available in naive random
#' forest approach}
#'
#' @examples
#' data(rct)
#' X <- rct$X;Y <- rct$Y;A <- rct$A
#' estimate_cate(X, A, Y, rep(mean(A), nrow(X)))$cate_coefficients[1:11,1]
#' @export

estimate_cate <-
function(
    X, A, Y, propensity_vec,
    m_of_x_mat = NULL,
    normalize_weights = TRUE,
    mu_x_a_coefs = NULL,
    init_cate_coefs = NULL,
    K = 5,
    outcome_model_obs_coefs = NULL,
    useRF = FALSE
) {
    if(!is.null(outcome_model_obs_coefs)) {
        est <- ss_estimate_cate_r_oscar(X, A, Y, propensity_vec,
            K, outcome_model_obs_coefs, m_of_x_mat, normalize_weights
        )
    } else if(!is.null(init_cate_coefs)) {
        est <- estimate_cate_r_oscar(X, A, Y, propensity_vec,
            init_cate_coefs, m_of_x_mat, normalize_weights
        )
    } else if(!is.null(mu_x_a_coefs)) {
        est <- estimate_cate_oscar(X, A, Y, propensity_vec, mu_x_a_coefs)
    } else if(useRF) {
        est <- estimate_cate_with_rct_rf(X, A, Y, propensity_vec, m_of_x_mat, normalize_weights)
    } else {
        est <- estimate_cate_with_rct(X, A, Y, propensity_vec, m_of_x_mat, normalize_weights)
        if(sum(abs(as.matrix(est$cate_coefficients)[-1,]) > 1e-8) == 0) {
            warning('a constant treatment effect estimate was returned')
        }
    }
    list(
        cate_coefficients = est$cate_coefficients,
        cate_function = est$cate_function,
        fitted_cate = est$fitted_cate,
        delta_coef_mat = est$delta_coef_mat,
        cate_delta_coefficients = est$cate_delta_coefficients,
        cate_model = est$cate_model
    )
}
