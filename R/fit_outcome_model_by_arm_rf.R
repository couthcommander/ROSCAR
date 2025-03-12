#' Fit Outcome Model with Random Forest
#'
#' Fit outcome regression model for each treatment arm.
#'
#' @keywords internal
#' @param X NEEDINFO
#' @param A NEEDINFO
#' @param Y NEEDINFO

fit_outcome_model_by_arm_rf <-
function(X, A, Y)
{
  ## Convert X to data frame with column names
  X <- as.data.frame(X)
  names(X) <- paste0("Var", seq_len(ncol(X)))
  tGrid <- data.frame(mtry = floor(sqrt(ncol(X))))

  trt_idx <- tapply(seq_along(A), A, I)
  ## list to store trained models -- one for each treatment
  rf_models <- lapply(trt_idx, function(i) {
    ## define the training control
    train_control <- caret::trainControl(method = "cv", number = 5)

    ## train random forest model
    caret::train(
      y = Y[i],
      x = X[i,,drop = FALSE],
      method = 'rf',
      trControl = train_control,
      tuneGrid = tGrid,
      ntree = 50,
      nodesize = 25
    )
  })

  ## function to predict mu_x for new data
  return_regression_mat <- function(x)
  {
    vapply(rf_models, stats::predict, numeric(nrow(x)), newdata = x)
  }
  
  ## predict mu_x for the training data
  mu_x_mat <- return_regression_mat(X)
  
  list(mu_x = mu_x_mat, 
       outcome_models = rf_models, 
       outcome_function = return_regression_mat)
}
