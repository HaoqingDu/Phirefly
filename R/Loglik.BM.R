# Author: Haoqing Du
# Latest Editing Time: 20/03/2024

## function is to calculate the likelihood:
## - with the Brownian Motion model
## - under the situation that (i) single trait; (ii) single sample per taxa.


## Likelihood calculation
## the result gives out the log-value of the likelihood

#' Title
#'
#' @param phy
#' @param chr.values
#' @param anc.values
#' @param sig2
#'
#' @return
#' @export
#'
#' @examples

loglik.BM <- function(phy, chr.values, anc.values, sig2) {
  C <- vcv.BM(phy)
  v1 <- matrix(rep(1, ncol(C)), ncol = 1)
  likelihood <- exp(-1/2 * t(chr.values-anc.values %x% v1) %*%
                      solve(sig2*C) %*%
                      (chr.values-anc.values %x% v1)) / sqrt((2*pi)^ncol(C) * det(sig2*C))
  return(as.numeric(log(likelihood)))
}
