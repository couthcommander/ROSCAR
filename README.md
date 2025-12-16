
Install and load the R package.

``` r
if(!requireNamespace('ROSCAR', quietly = TRUE)) {
  remotes::install_github('couthcommander/ROSCAR')
}
library(ROSCAR)
```

Generate some data (both RCT and OS). In this case, we’ll introduce
covariate mismatch.

``` r
set.seed(2025)
dat <- simulate_rct_and_os_data(
  n_r = 300, n_o = 1000, support_fraction = 1/20,
  covariate_effect = c(Z = 2/3, V = 1),
  frac_U = 0.3, frac_V = 0.3
)
```

Take out necessary components.

``` r
rct <- list(X = dat$X_RCT, Y = dat$y_r$y_observed, A = dat$a_r)
os <- list(X = dat$X_OS, Y = dat$y_o$y_observed, A = dat$a_o)
tau <- dat$y_r$tau
```

Run and evaluate models.

``` r
mod <- cate_model(rct, os)
head(model_pred(mod))
```

    ##        naive      racer      oscar    r_oscar
    ## 1  0.5969415 -0.4252329 -0.9880615 -0.8910070
    ## 2 -0.3047681 -0.5098061 -0.5721287 -0.3644886
    ## 3  2.4979136  2.8464549  2.7485553  3.0285358
    ## 4  1.1565326  1.5171946  1.5858946  1.6751391
    ## 5 -1.1141061 -0.6689769 -0.3313598 -0.4118794
    ## 6  1.0056099  0.9140123  0.5881268  0.7044381

``` r
head(model_pred(mod, rct$X))
```

    ##        naive      racer      oscar    r_oscar
    ## 1  0.5969415 -0.4252329 -0.9880615 -0.8910070
    ## 2 -0.3047681 -0.5098061 -0.5721287 -0.3644886
    ## 3  2.4979136  2.8464549  2.7485553  3.0285358
    ## 4  1.1565326  1.5171946  1.5858946  1.6751391
    ## 5 -1.1141061 -0.6689769 -0.3313598 -0.4118794
    ## 6  1.0056099  0.9140123  0.5881268  0.7044381

``` r
model_eval(mod, tau)
```

    ##             rmse rank_corr accuracy_itr cal_intercept cal_slope
    ## naive   1.207203 0.6072023    0.7266667    -0.1230678  1.080000
    ## racer   1.052028 0.7330681    0.7733333    -0.2038329  1.108677
    ## oscar   1.019060 0.7425060    0.8000000    -0.1032908  1.060239
    ## r_oscar 1.022510 0.7424229    0.8033333    -0.1837056  1.076606

Repeat with shared covariates.

``` r
rct_sh <- rct
rct_sh[['X']] <- do.call(cbind, dat$components[c('U','Z_r')])
os_sh <- os
os_sh[['X']] <- do.call(cbind, dat$components[c('Z_o','V')])
mod_sh <- cate_model(rct_sh, os_sh)
model_eval(mod_sh, tau)
```

    ##             rmse rank_corr accuracy_itr cal_intercept cal_slope
    ## naive   1.197412 0.6130806    0.7366667   -0.09870933 1.0239523
    ## racer   1.042412 0.7352517    0.7800000   -0.15499987 1.0069535
    ## oscar   1.061850 0.7179022    0.7833333   -0.06819577 0.9679672
    ## r_oscar 1.051774 0.7267010    0.7833333   -0.10839042 0.9538993
