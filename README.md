
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
rct_os <- build_data(dat = dat)
```

Run and evaluate models.

``` r
mod <- cate_model(rct_os$RCT, rct_os$OS)
head(model_pred(mod))
```

    ##        naive     racer     oscar   r_oscar
    ## 1 -0.5629661 -1.079727 -1.470010 -1.784453
    ## 2  0.7013522  1.745146  2.370394  2.409035
    ## 3  2.8328786  2.588514  2.996943  2.873513
    ## 4 -1.6335755 -1.183861 -1.538120 -1.931878
    ## 5 -1.2256149 -1.768328 -1.155722 -1.030950
    ## 6  1.5433424  1.384840  1.976091  2.388333

``` r
model_eval(mod, rct_os$tau)
```

    ##             rmse rank_corr accuracy_itr cal_intercept cal_slope
    ## naive   1.511083 0.7071487    0.7966667   0.008642961 1.1356766
    ## racer   1.385921 0.7641734    0.8066667  -0.143827935 1.1413424
    ## oscar   1.367726 0.7677134    0.8066667  -0.161038795 1.0111305
    ## r_oscar 1.366732 0.7693170    0.8000000  -0.078221907 0.9613386

Repeat with only shared covariates.

``` r
rct_os <- build_data(dat = dat, RCT_imp_method = NULL, OS_imp_method = NULL)
mod_sh <- cate_model(rct_os$RCT, rct_os$OS)
model_eval(mod_sh, rct_os$tau)
```

    ##             rmse rank_corr accuracy_itr cal_intercept cal_slope
    ## naive   1.757951 0.5436260    0.6833333    0.02921390  1.002343
    ## racer   1.785060 0.5313090    0.7000000   -0.06055680  1.168131
    ## oscar   1.736864 0.5691734    0.7166667   -0.22671934  1.085618
    ## r_oscar 1.730254 0.5679761    0.7300000   -0.01768978  1.047507

Try a quick experiment.

``` r
iter <- 10
set.seed(iter)
l <- vector('list', iter)
for(i in seq_along(l)) {
  dat <- simulate_rct_and_os_data(n_r=250, n_o=500)
  rct_os <- build_data(dat = dat)
  mod <- cate_model(rct_os$RCT, rct_os$OS)
  l[[i]] <- model_eval(mod, dat$y_r$tau)
}
# mean RMSE for each method
rowMeans(vapply(l, \(i) i[,'rmse'], numeric(4)))
```

    ##     naive     racer     oscar   r_oscar 
    ## 1.3424593 0.3153759 0.2095290 0.2080282
