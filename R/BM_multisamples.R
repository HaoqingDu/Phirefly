# Author: Haoqing Du
# Latest Editing Time: 03/05/2024

## function is to calculate the likelihood:
## - with the Brownian Motion model
## - under the situation that (i) multiple independent traits; (ii) multi-samples per taxa.
## - using pruning methods

pruning_BM_multi <- function(td, trait_names, sigma2, tau, nsamples, node_index) {

  phy <- td@phylo
  # if(class(phy) != "phylo") {stop(td," does not have a \" phylo \".")}

  ntaxa <- length(phy$tip.label)
  nchr <- length(trait_names)
  characters <- as.matrix(td@data[trait_names][1:ntaxa,])
  descendants <- phy$edge[phy$edge[,1] == node_index,][,2]

  # inv_sigma2 <- solve(sigma2)

  if (descendants[1] <= ntaxa) { # if left descendant is a tip
    edge_left <- phy$edge.length[phy$edge[,2] == descendants[1]] +
      tau[descendants[1],]/nsamples[descendants[1],]
    # edge_left <- matrix(phy$edge.length[phy$edge[,2] == descendants[1]], nchr, nchr) +
    #   diag(tau[descendants[1],]/nsamples[descendants[1],], nchr, nchr, names = T)
    v_left <- sigma2 * edge_left

    chr_left <- characters[descendants[1],]
    # likelihood of sample mean
    loglik_left <- - 1/2 * nchr * log(2*pi) -
      1/2 * log(det(sigma2 * tau[descendants[1],]))
    # 1/2 * log(det(sigma2 * diag(tau[descendants[1],], nchr, nchr, names = T)))
    # - 1/2 * t(chr_left - chr_right) %*% solve(v_left + v_right) %*% (chr_left - chr_right)
    # loglik_left <- sum(dnorm(sample, mean, sd ,log = T))
  }
  else { # if left offspring node is not a tip
    recursive <- pruning_BM_multi(td, trait_names, sigma2, tau, nsamples, descendants[1])

    edge_left <- phy$edge.length[phy$edge[,2] == descendants[1]]+
      recursive$extended.edge
    v_left <- sigma2 * edge_left

    chr_left <- recursive$chr.values
    loglik_left <- recursive$loglik
  }

  if (descendants[2] <= ntaxa) { # if right offspring node is a tip
    edge_right <- phy$edge.length[phy$edge[,2] == descendants[2]] +
      tau[descendants[2],]/nsamples[descendants[2],]
    # edge_right <- matrix(phy$edge.length[phy$edge[,2] == descendants[2]], nchr, nchr) +
    #   diag(tau[descendants[2],]/nsamples[descendants[2],], nchr, nchr, names = T)
    v_right <- sigma2 * edge_right

    chr_right <- characters[descendants[2],]
    # likelihood of sample mean
    loglik_right <- - 1/2 * nchr * log(2*pi) -
      1/2 * log(det(sigma2 * tau[descendants[2],]))
    # 1/2 * log(det(sigma2 * diag(tau[descendants[2],], nchr, nchr, names = T)))
    # loglik_left <- sum(dnorm(sample, mean, sd ,log = T))
  }
  else { # right offspring note is not a tip
    recursive <- pruning_BM_multi(td, trait_names, sigma2, tau, nsamples, descendants[2])

    edge_right <- phy$edge.length[phy$edge[,2] == descendants[2]]+
      recursive$extended.edge
    v_right <- sigma2 * edge_right

    chr_right <- recursive$chr.values
    loglik_right <- recursive$loglik
  }

  #extended branch length and variance
  edge <- (edge_left * edge_right) / (edge_left + edge_right)
  # print(v_left)
  # print(v_right)
  # v_node <- (v_left * v_right) / (v_left + v_right)
  # print(v_node)

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



loglik_BM_multi <- function(td, trait_names, mu, sigma2, tau, nsamples) {

  phy <- td@phylo
  if(class(phy) != "phylo") {stop(td," does not have a \" phylo \".")}

  ntaxa <- length(phy$tip.label)
  nchr <- length(trait_names)
  # characters <- td@data[[trait_names]][1:ntaxa]
  root_index <- ntaxa + 1

  root <- pruning_BM_multi(td, trait_names, sigma2, tau, nsamples, root_index)
  log.likelihood <- root$loglik - 1/2 * nchr * log(2*pi) -
    1/2 * log(det(sigma2 * root$extended.edge)) -
    1/2 * t(root$chr.values - mu) %*% solve(sigma2 * root$extended.edge) %*% (root$chr.values - mu)

  return(log.likelihood)
}
