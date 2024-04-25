# Author: Haoqing Du
# Latest Editing Time: 28/03/2024

## function is to calculate the likelihood:
## - with the Brownian Motion model
## - under the situation that (i) multiple independent traits; (ii) single sample per taxa.
## - using pruning methods


## pruning
#' Title
#'
#' @param td
#' @param trait_names
#' @param sigma2
#' @param node_index
#'
#' @return
#' @export
#'
#' @examples
pruning_BM <- function(td, trait_names, sigma2, node_index) {

  phy <- td@phylo
  # if(class(phy) != "phylo") {stop(td," does not have a \" phylo \".")}

  ntaxa <- length(phy$tip.label)
  nchr <- length(trait_names)
  characters <- as.matrix(td@data[trait_names][1:ntaxa,])
  descendants <- phy$edge[phy$edge[,1] == node_index,][,2]

  if (descendants[1] <= ntaxa) { # if left descendant is a tip
    edge_left <- phy$edge.length[phy$edge[,2] == descendants[1]]
    v_left <- sigma2 * edge_left
    chr_left <- characters[descendants[1],]
    loglik_left <- 0
  }
  else { # if left offspring node is not a tip
    recursive <- pruning_BM(td, trait_names, sigma2, descendants[1])

    edge_left <- phy$edge.length[phy$edge[,2] == descendants[1]]+
      recursive$extended.edge
    v_left <- sigma2 * edge_left

    chr_left <- recursive$chr.values
    loglik_left <- recursive$loglik
  }

  if (descendants[2] <= ntaxa) { # if right offspring node is a tip
   edge_right <- phy$edge.length[phy$edge[,2] == descendants[2]]
  #  v_right <- sigma2 * edge_right
    chr_right <- characters[descendants[2],]
    loglik_right <- 0
  }
  else { # right offspring note is not a tip
    recursive <- pruning_BM(td, trait_names, sigma2, descendants[2])

    edge_right <- phy$edge.length[phy$edge[,2] == descendants[2]]+
      recursive$extended.edge
    v_right <- sigma2 * edge_right

    chr_right <- recursive$chr.values
    loglik_right <- recursive$loglik
  }

  #extended branch length
  edge <- (edge_left * edge_right) / (edge_left + edge_right)

  # character value for the node is the weighted average of two descendants
  chr.values <- (edge_left*chr_right + edge_right*chr_left) / (edge_left + edge_right)

  # log-likellihood = sum of left & right descendants' likelihood and the likelihood of itself
  log.likelihood <- loglik_left + loglik_right
  log.likelihood <- log.likelihood - 1/2 * nchr * log(2*pi) -
    1/2 * log(det(v_left + v_right)) -
    1/2 * t(chr_left - chr_right) %*% solve(v_left + v_right) %*% (chr_left - chr_right)

  # log.likelihood <- log.likelihood - 1/2 * log(2*pi*sigma2*(edge_left + edge_right)) -
  #   1/2 * (chr_left - chr_right)^2 / (sigma2 * (edge_left + edge_right))


  return(list(extended.edge = edge,
              chr.values = chr.values,
              loglik = log.likelihood))
}

## log likelihood printing

#' Title
#'
#' @param td
#' @param trait_names
#' @param mu
#' @param sigma2
#'
#' @return
#' @export
#'
#' @examples
loglik_BM_prun <- function(td, trait_names, mu, sigma2) {

  phy <- td@phylo
  if(class(phy) != "phylo") {stop(td," does not have a \" phylo \".")}

  ntaxa <- length(phy$tip.label)
  nchr <- length(trait_names)
  # characters <- td@data[[trait_names]][1:ntaxa]
  root_index <- ntaxa + 1

  root <- pruning_BM(td, trait_names, sigma2, root_index)
  log.likelihood <- root$loglik - 1/2 * nchr * log(2*pi) -
    1/2 * log(det(sigma2 * root$extended.edge)) -
    1/2 * t(root$chr.values - mu) %*% solve(sigma2 * root$extended.edge) %*% (root$chr.values - mu)

  return(log.likelihood)
}
