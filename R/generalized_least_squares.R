#' Generalized least squares
#'
#' @param V residual among-species covariance matrix
#' @param X model matrix
#' @param y the observations
#' @param y the parameter vector
#'
#' @return the log likelihood value
#' @export
#'
#' @examples
#'
#' n_obs <- 10
#'
#' V <- crossprod(matrix(rnorm(n_obs^2),n_obs,n_obs))
#' X <- matrix(1,n_obs,1)
#' y <- rnorm(n_obs, mean = 5, sd = 1)
#' betas <- seq(3, 7, length.out = 20)
#'
#' logls <- sapply(betas, function(beta) generalized_least_squares(V, X, y, beta))
#'
#' plot(betas, logls)
generalized_least_squares <- function(V, X, y, beta1){
  # Generalised least square method (aka variance-covariance matrix method)
  ntip <- nrow(V)
  L = t(chol(V)) # lower triangular matrix

  log_det_V <- 2*sum(log(diag(L)))

  r = forwardsolve(L, y - X %*% beta1) # what does de-correlated residuals mean?

  dot_r <- (t(r) %*% r)[1]
  logl = -0.5 * ntip * log(2*pi) -0.5 * log_det_V - 0.5 * dot_r
  return(logl)
}
