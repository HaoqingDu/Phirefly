## Variance-Covariance matrix for a phylogeny
## Model: Brownian Motion

vcv.BM <- function(phy) {
  if(class(phy) != "phylo") stop(phy," is not a \" phylo \".")
  C <- matrix(NA, nrow = length(phy$tip.label), ncol = length(phy$tip.label))
  for (n in 1:ncol(C)) {
    for (m in 1:n) {
      if (m!=n) {
        C[m,n] <- from.root(phy, common.ancestor(phy, m, n))
        C[n,m] <- C[m,n]
      } else {
        C[m,n] <- from.root(phy, phy$edge[phy$edge[,2] == m, 1]) + phy$edge.length[phy$edge[,2] == m]
      }
    }
  }
  return(C)
}

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
