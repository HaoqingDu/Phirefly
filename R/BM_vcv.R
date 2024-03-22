# Author: Haoqing Du
# Latest Editing Time: 20/03/2024

## function is to calculate the likelihood:
## - with the Brownian Motion model
## - under the situation that (i) single trait; (ii) single sample per taxa.

## Likelihood calculation
## the result gives out the log-value of the likelihood

#' Calculate the likelihood of a univariate BM process using the vcv method
#'
#' @inherit logl_OU_fitzjohn
#'
#' @export
#'
#' @examples
#' data("artiodactyla")
#' logl_BM_vcv(artiodactyla, 0.5, 0.5, "brain_mass_g_log_mean")
logl_BM_vcv <- function(td, mu, sigma2, trait_name) {
  phy <- td@phylo
  ntip <- length(phy$tip.label)
  chr.values <- td@data[[trait_name]][1:ntip]

  C <- BM_vcv(phy)
  v1 <- matrix(rep(1, ncol(C)), ncol = 1)
  likelihood <- exp(-1/2 * t(chr.values - mu %x% v1) %*%
                      solve(sigma2*C) %*%
                      (chr.values-mu %x% v1)) / sqrt((2*pi)^ncol(C) * det(sigma2*C))
  return(as.numeric(log(likelihood)))
}

## Variance-Covariance matrix for a phylogeny
## Model: Brownian Motion

BM_vcv <- function(phy) {
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



