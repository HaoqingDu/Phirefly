#' Generalized least squares
#'
#' @param V residual among-species covariance matrix
#' @param X model matrix
#'
#' @return the log likelihood value
#' @export
#'
#' @examples
#'
#' V <- crossprod(matrix(rnorm(16),4,4))
#' X <- matrix(1,4,1)
#'
#' generalized_least_squares(V, X)
generalized_least_squares <- function(V, X){
  ntip <- nrow(V)
  L = t(chol(V)) # lower triangular matrix

  log_det_V <- 2*sum(log(diag(L)))

  r = solve(L) %*% y - solve(L) %*% X * theta # what does de-correlated residuals mean?

  logl = -0.5 * ntip * log(2*pi) -0.5 * log_det_V - 0.5 * dot(r, r)
  return(logl)
}
