#' estimate_cate
#'
#' Estimate conditional average treatment effect with naive/RACER random forest approach.
#'
#' @keywords internal
#' @param X NEEDINFO
#' @param A NEEDINFO
#' @param Y NEEDINFO
#' @param propensity_vec NEEDINFO
#' @param m_of_x_mat NEEDINFO
#' @param normalize_weights NEEDINFO

estimate_cate_with_rct_rf <-
function(X, A, Y, propensity_vec, 
                             m_of_x_mat = NULL, ## regression functions for each trt level
                             normalize_weights = TRUE)
{
  stopifnot("A must be either 0 or 1" = all(sort(unique(A)) == c(0,1)))

  ## inverse propensity score weights
  propens_weights <- ifelse(A == 1, 1/propensity_vec, 1/(1-propensity_vec))

  ## turn binary A into +1 and -1 variable
  A_1_m1 <- (2 * A - 1)

  if (normalize_weights)
  {
    mean_propens <- mean(propens_weights)
    propens_weights <- propens_weights / mean_propens
  } else
  {
    mean_propens <- 1
  }

  if (is.null(m_of_x_mat))
  {
    m_of_x_avg <- numeric(NROW(X))
  } else
  {
    ## compute shift function
    m_of_x_avg <- rowSums(m_of_x_mat * cbind(propensity_vec, 1 - propensity_vec))
  }
  ## Convert X to data frame with column names
  X <- as.data.frame(X)
  names(X) <- paste0("Var", seq_len(ncol(X)))

  ## construct working response
  y_tilde <- (A_1_m1 * (Y - m_of_x_avg) * propens_weights) * mean_propens

  ## define the training control
  train_control <- caret::trainControl(method = "cv", number = 5)

  ## estimate CATE by predicting working response with covariates using random forest
  cate_model <- caret::train(y = y_tilde, x = X, method = "rf", trControl = train_control,
                      tuneGrid = data.frame(mtry = floor(sqrt(ncol(X)))), ntree = 50, nodesize = 25)

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

  list(
    cate_model = cate_model, 
    cate_function = return_cate, 
    fitted_cate = return_cate(X)
  )
}
