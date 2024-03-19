# Author: Haoqing Du
# Latest Editing Time: 19/03/2024

## function is to calculate the likelihood:
## - with the Brownian Motion model
## - under the situation that (i) single trait; (ii) single sample per taxa.


## Ancestral value estimates
## the maximum likelihood estimate of the ancestral values can be calculated (O’Meara et al. 2006)

anc.values <- function(phy, chr.values) {
  C <- vcv.BM(phy)
  v1 <- rep(1, ncol(C))
  a <- solve(t(v1) %*% solve(C) %*% v1) %*% (t(v1) %*% solve(C) %*% chr.values)
  print(a)
}


## Likelihood calculation
## the result gives out the log-value of the likelihood

Loglik.BM <- function(phy, chr.values, anc.values, sig2) {
  C <- vcv.BM(phy)
  v1 <- matrix(rep(1, ncol(C)), ncol = 1)
  likelihood <- exp(-1/2 * t(chr.values-anc.values %x% v1) %*%
                      solve(sig2*C) %*%
                      (chr.values-anc.values %x% v1)) / sqrt((2*pi)^ncol(C) * det(sig2*C))
  print(log(likelihood))
}


## Maximum likelihood estimate of variance per unit time (sigma2)
## (O’Meara et al. 2006)

sigma2.BM <- function(phy, chr.values, anc.values) {
  C <- vcv.BM(phy)
  v1 <- matrix(rep(1, ncol(C)), ncol = 1)
  sig2.hat <- t(chr.values-anc.values %x% v1) %*% solve(C) %*%
    (chr.values-anc.values %x% v1)/ncol(C)
  print(sig2.hat)
}


## Maximum likelihood

ML.BM <- function(phy, chr.values) {
  z_0 <- anc.values(phy, chr.values) # ML estimates of ancestral values
  sigma2.hat <- as.numeric(sigma2.BM(phy, chr.values, z_0)) # ML estimates of variance

  ML <- Loglik.BM(phy, chr.values, z_0, sigma2.hat)

  print(ML)
}



