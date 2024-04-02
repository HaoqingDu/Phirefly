# Author: Haoqing Du
# Latest Editing Time: 26/03/2024

## function is to calculate the likelihood:
## - under the situation that (i) single trait; (ii) single sample per taxa.
## - using vcv methods

## vcv matrix
## appliable to both BM and OU

vcv.matrix <- function(td,
                       model = c("BM", "OU"),
                       alpha = NULL) {
  phy <- td@phylo
  if(class(phy) != "phylo") stop(td," does not have a \" phylo \".")

  ntaxa <- length(phy$tip.label)
  V <- matrix(NA, nrow = ntaxa, ncol = ntaxa)

  if(length(model) > 1) {
    warning("Please specify the model! (\"BM\"/\"OU\")")
    print("The given result is under the BM model")
    model <- "BM"
  }

  if(model == "OU") {
    # alpha is required under OU process
    if(!is.null(alpha)) {
      for (n in 1:ntaxa) {
        for (m in 1:m) {
          # C[i,j] = C[j,i] = 1/(2alpha)*exp[-alpha*t_ij]*[1-exp[-2alpha*t_ra]]
          # t_ij = distance between i & j, t_ra = distance between root and mrca
          V[m,n] <- 1/(2*alpha)*
            exp(-alpha*dist(phy, m, n))*
            (1-exp(-2*alpha*from.root(phy, common.ancestor(phy, m, n))))
          V[n,m] <- C[m,n]
        }
      }
    } else {
      warning("alpha is required for the OU model!")
      # when alpha is missing under the OU model, the function switch to BM automatically
      print("The given result is under the BM model.")
      vcv.matrix(td,
                 model <- "BM")
    }
  }

  if(model == "BM") {
    # when model is not specified, use BM.
    for (n in 1:ntaxa) {
      for (m in 1:n) {
        if (m!=n) {
          V[m,n] <- from.root(phy, common.ancestor(phy, m, n))
          V[n,m] <- V[m,n]
        } else {
          V[m,n] <- from.root(phy, phy$edge[phy$edge[,2] == m, 1]) + phy$edge.length[phy$edge[,2] == m]
        }
      }
    }
  }

  return(V)
}

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
