#Perturb m elements of a vector by amount
#' @keywords internal
perturb <- function(vector, ind, amount){
  m <- length(ind)
  vector[ind] <- vector[ind] + sample(c(-1,1), replace=T, size=m) * stats::runif(m, min=.5*amount, max=1*amount)
  return(vector)
}

#' @keywords internal
linearFun <- function(x, b) {
  x %*% b
}

#' @keywords internal
quadraticFun <- function(x, b) {
  p <- length(b)
  m <- stats::runif(p, min = 0.25, max = 0.5)
  m[b == 0] <- 0
  x %*% b + x^2 %*% m
}

#' @keywords internal
sineFun <- function(x, b) {
  p <- length(b)
  m <- stats::runif(p, min = 0.25, max = 0.5)
  m[b == 0] <- 0
  x %*% b + sin(x) %*% m
}

# Simulate outcome model for potential and observed outcomes and compute individual treatment effect
#' @keywords internal
outcome_model <- function(x, b0, b1, sigma, treatment_assigned, FUN = linearFun) {
  n <- nrow(x)
  stopifnot('user provided FUN must have arguments "x" and "b"' = identical(c('x','b'), methods::formalArgs(FUN)))
  y0 <- drop(FUN(x, b0))
  y1 <- drop(FUN(x, b1))
  tau <- y1 - y0
  y1 <- y1 + stats::rnorm(n, sd=.5*sigma)
  y0 <- y0 + stats::rnorm(n, sd=.5*sigma)

  y_potential <- cbind(y0, y1)
  colnames(y_potential) <- c("untreated", "treated")
  y_observed <- y_potential[cbind(seq_len(n), treatment_assigned+1)]
  list(y_potential=y_potential, y_observed=y_observed, tau=tau)
}
