#' Build RCT and OS
#'
#' Create necessary RCT and OS components required to
#' run \code{\link{cate_model}}.
#'
#' @param RCT Randomized clincal trial (RCT) data
#' is a list with U (covariates not available in OS),
#' Z (shared covariates with OS), Y (outcome), A (treatment).
#' @param OS Observational study (OS) data is a list
#' with V (covariates not available in RCT),
#' Z (shared covariates with RCT), Y, A.
#' @param dat Instead of \sQuote{RCT} and \sQuote{OS}, a \sQuote{sim_data}
#' object can be provided (see \code{\link{simulate_rct_and_os_data}}).
#' @param RCT_imp_method Imputation function for RCT data;
#' when unspecified \sQuote{Vhat_r} will be imputed with \sQuote{OLS}.
#' Set to NULL to only use shared covariates (no imputation).
#' @param OS_imp_method Imputation function for OS data; when
#' unspecified \sQuote{Uhat_o} will be imputed as zeroes.
#' Set to NULL to only use shared covariates (no imputation).
#' @param has_excl Required for imputation method
#' @param has_shr Required for imputation method
#' @param oth_shr Required for imputation method
#' @param oth_excl Required for imputation method
#'
#' Imputation functions can be provided to impute either with RCT or
#' OS data. In either case, it's expected that the function will have four
#' arguments (has_excl, has_shr, oth_shr, oth_excl). \sQuote{has_excl}
#' and \sQuote{has_shr} refer to the data set receiving imputation.
#' \sQuote{oth_shr} and \sQuote{oth_excl} refer to the data set being
#' used for imputation. In the typical setting of imputing RCT covariates
#' from OS, has_excl = RCT[['U']], has_shr = RCT[['Z']], oth_shr = OS[['Z']],
#' and oth_excl = OS[['V']]. A user-defined function is limited to this
#' information.
#'
#' @return
#' \item{components of RCT and OS data}{
#' \describe{
#' \item{RCT}{covariates (X), outcome (Y), and treatment (A) for RCT data}
#' \item{OS}{covariates (X), outcome (Y), and treatment (A) for OS data}
#' \item{tau}{True treatment effect for RCT data, only available
#' with \sQuote{sim_data} method.}
#' }}
#'
#' @examples
#' data(os, rct)
#' rct1 <- list(Z = rct$X, Y = rct$Y, A = rct$A)
#' os1 <- list(Z = os$X, Y = os$Y, A = os$A)
#' ex1 <- build_data(rct1, os1)
#'
#' # simulate with covariate mismatch
#' x <- simulate_rct_and_os_data(
#'   n_r = 300, n_o = 1000, support_fraction = 1/20,
#'   covariate_effect = c(Z = 2/3, V = 1),
#'   frac_U = 0.3, frac_V = 0.3
#' )
#' rct_os <- build_data(x)
#'
#' @export

build_data <- function(RCT, OS, RCT_imp_method, OS_imp_method) {
    UseMethod('build_data')
}

#' @rdname build_data
#' @export
build_data.default <- function(RCT, OS, RCT_imp_method = impute_OLS, OS_imp_method = impute_0) {
    U_r <- RCT[['U']]
    Z_r <- RCT[['Z']]
    Z_o <- OS[['Z']]
    V_o <- OS[['V']]
    RCT_Y <- RCT[['Y']]
    RCT_A <- RCT[['A']]
    OS_Y <- OS[['Y']]
    OS_A <- OS[['A']]
    stopifnot(
        'RCT list should contain "Y" (outcome) element' = !is.null(RCT_Y),
        'RCT list should contain "A" (treatment) element' = !is.null(RCT_A),
        'RCT list should contain "Z" (and "U")' = !is.null(Z_r),
        'OS list should contain "Y" (outcome) element' = !is.null(OS_Y),
        'OS list should contain "A" (treatment) element' = !is.null(OS_A),
        'OS list should contain "Z" (and "V")' = !is.null(Z_o),
        'same number of shared covariates "Z" required for RCT and OS' = ncol(Z_r) == ncol(Z_o)
    )
    if(is.null(RCT_imp_method) && is.null(OS_imp_method)) {
        RCT_cov <- Z_r
        OS_cov <- Z_o
    } else {
        RCT_cov <- cbind(U_r, Z_r)
        OS_cov <- cbind(Z_o, V_o)
        if(!is.null(V_o) && !is.null(RCT_imp_method)) {
            imp_dim <- c(nrow(RCT_cov), ncol(V_o))
            Vhat_r <- RCT_imp_method(U_r, Z_r, Z_o, V_o)
            stopifnot(identical(imp_dim, dim(Vhat_r)))
            RCT_cov <- cbind(RCT_cov, Vhat_r)
        }
        if(!is.null(U_r) && !is.null(OS_imp_method)) {
            imp_dim <- c(nrow(OS_cov), ncol(U_r))
            Uhat_o <- OS_imp_method(V_o, Z_o, Z_r, U_r)
            stopifnot(identical(imp_dim, dim(Uhat_o)))
            OS_cov <- cbind(Uhat_o, OS_cov)
        }
        if(ncol(RCT_cov) != ncol(OS_cov)) {
            warning(call. = FALSE,
'RCT and OS data have an unequal number of covariates after imputation.
Limiting data to shared covariates.')
            RCT_cov <- Z_r
            OS_cov <- Z_o
        }
    }
    rct <- list(X = RCT_cov, Y = RCT_Y, A = RCT_A)
    os <- list(X = OS_cov, Y = OS_Y, A = OS_A)
    list(RCT = rct, OS = os)
}

#' @rdname build_data
#' @export
build_data.sim_data <- function(dat, RCT_imp_method = impute_OLS, OS_imp_method = impute_0) {
    comp <- dat$components
    if(is.null(comp)) {
        comp <- list(Z_r = dat$X_RCT, Z_o = dat$X_OS)
        RCT_imp_method <- NULL
        OS_imp_method <- NULL
    }
    RCT <- list(U = comp$U, Z = comp$Z_r, Y = dat$y_r$y_observed, A = dat$a_r)
    OS <- list(Z = comp$Z_o, V = comp$V, Y = dat$y_o$y_observed, A = dat$a_o)
    o <- build_data(RCT, OS, RCT_imp_method, OS_imp_method)
    c(o, tau = list(dat$y_r$tau))
}

#' @rdname build_data
#' @export
impute_OLS <- function(has_excl, has_shr, oth_shr, oth_excl) {
    ## Fit V ~ Z on OS (one OLS per column; closed-form matrix works too)
    B_hat <- MASS::ginv(t(oth_shr) %*% oth_shr) %*% (t(oth_shr) %*% oth_excl)
    has_shr %*% B_hat
}

#' @rdname build_data
#' @export
impute_0 <- function(has_excl, has_shr, oth_shr, oth_excl) {
    matrix(0, nrow(has_shr), ncol(oth_excl))
}
