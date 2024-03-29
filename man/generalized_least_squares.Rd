% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generalized_least_squares.R
\name{generalized_least_squares}
\alias{generalized_least_squares}
\title{Generalized least squares}
\usage{
generalized_least_squares(V, X, y, beta1)
}
\arguments{
\item{V}{residual among-species covariance matrix}

\item{X}{model matrix}

\item{y}{the parameter vector}
}
\value{
the log likelihood value
}
\description{
Generalized least squares
}
\examples{

n_obs <- 10

V <- crossprod(matrix(rnorm(n_obs^2),n_obs,n_obs))
X <- matrix(1,n_obs,1)
y <- rnorm(n_obs, mean = 5, sd = 1)
betas <- seq(3, 7, length.out = 20)

logls <- sapply(betas, function(beta) generalized_least_squares(V, X, y, beta))

plot(betas, logls)
}
