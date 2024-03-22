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
#' @inherit logl_BM_fitzjohn
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

  V <- BM_vcv(phy, sigma2)
  X <- matrix(1, ntip, 1)
  beta1 <- matrix(mu, 1,1)

  logl <- generalized_least_squares(V, X, chr.values, beta1)
  return(logl)
}

## Variance-Covariance matrix for a phylogeny
## Model: Brownian Motion

BM_vcv <- function(phy, sigma2) {
  if(class(phy) != "phylo") stop(phy," is not a \" phylo \".")
  V <- matrix(NA, nrow = length(phy$tip.label), ncol = length(phy$tip.label))
  for (n in 1:ncol(V)) {
    for (m in 1:n) {
      if (m!=n) {
        V[m,n] <- from.root(phy, common.ancestor(phy, m, n))
        V[n,m] <- V[m,n]
      } else {
        V[m,n] <- from.root(phy, phy$edge[phy$edge[,2] == m, 1]) + phy$edge.length[phy$edge[,2] == m]
      }
    }
  }
  V <- V * sigma2
  return(V)
}



