# Author: Haoqing Du
# Latest Editing Time: 20/03/2024

## Under the BM model, the maximum likelihood estimates of mu (root values) and sigma2 (variance per unit time)
## can actually be calculated analytically (O’Meara et al. 2006).


## Root value estimates

root.values <- function(phy, chr.values) {
  C <- vcv.BM(phy)
  v1 <- rep(1, ncol(C))
  a <- solve(t(v1) %*% solve(C) %*% v1) %*% (t(v1) %*% solve(C) %*% chr.values)
  return(a)
}


## Maximum likelihood estimate of variance per unit time (sigma2)

sigma2.BM <- function(phy, chr.values, anc.values) {
  C <- vcv.BM(phy)
  v1 <- matrix(rep(1, ncol(C)), ncol = 1)
  sig2.hat <- t(chr.values-anc.values %x% v1) %*% solve(C) %*%
    (chr.values-anc.values %x% v1)/ncol(C)
  return(sig2.hat)
}
