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
#' V <- crossprod(matrix(rnorm(16),4,4))
#' X <- matrix(1,4,1)
#' y <- c(0.1, 0.3, -0.3, 0.0)
#' beta1 <- 0.2
#'
#' generalized_least_squares(V, X, y, beta1)
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
