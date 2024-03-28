# Author: Haoqing Du
# Latest Editing Time: 26/03/2024

## function is to calculate the likelihood:
## - with the Brownian Motion model
## - under the situation that (i) single trait; (ii) single sample per taxa.
## - using vcv methods


## Likelihood calculation
## the result gives out the log-value of the likelihood

#' Title
#'
#' @param sig2
#' @param td
#' @param trait.names
#' @param mu
#' @param model
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples

loglik_vcv <- function(td, trait_names,
                   mu, sig2,
                   model = C("BM", "OU"),
                   alpha = NULL) {

  phy <- td@phylo
  if(class(phy) != "phylo") {stop(td," does not have a \" phylo \".")}

  ntaxa <- length(phy$tip.label)
  chr.values <- td@data[[trait_names]][1:ntaxa]

  if(length(model) > 1) {
    warning("Please specify the model! (\"BM\"/\"OU\")")
    print("The given result is under the BM model")
    model <- "BM"
  }

  m <- model

  V <- vcv.matrix(td,
                  model = m,
                  alpha)
  v1 <- matrix(1, ntaxa, 1)

  # loglikelihood = -1/2 * (X - mu)^T (sigma2 C)^-1 (X - mu)
  #                 -1/2 * n * log(2pi) _ 1/2 log(det(sigma2 C))
  log.likelihood <- -1/2 * t(chr.values-mu * v1) %*%
    solve(sig2*V) %*% (chr.values-mu * v1) -
    1/2*ntaxa*log(2*pi) - 1/2*log(det(sig2*V))

  return(as.numeric(log.likelihood))
  }
