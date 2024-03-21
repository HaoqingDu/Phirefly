# Author: Haoqing Du
# Latest Editing Time: 21/03/2024

## function is to calculate the likelihood:
## - with the Brownian Motion model
## - under the situation that (i) single trait; (ii) single sample per taxa.


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

loglik <- function(td, trait.names,
                   mu, sig2,
                   model = C("BM", "OU"),
                   alpha = NULL) {

  if(class(td@phylo) != "phylo") {stop(td," does not have a \" phylo \".")}

  ntaxa <- length(td@phylo$tip.labels)
  chr.values <- td@data$trait.names

  C <- vcv.matrix(td,
                  model = model,
                  alpha)
  v1 <- matrix(rep(1, ntaxa), ncol = 1)
  log.likelihood <- -1/2 * t(chr.values-mu %x% v1) %*%
    solve(sig2*C) %*% (chr.values-mu %x% v1) -
    1/2*ntaxa*log(2*pi) - 1/2*det(sig2*C)
  return(as.numeric(log.likelihood))
  }
