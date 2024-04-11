# Author: Haoqing Du
# Latest Editing Time: 08/04/2024

## function is to calculate the likelihood:
## - with the Brownian Motion model
## - under the situation that (i) multiple independent traits; (ii) multi-samples per taxa.
## - using pruning methods

pruning_multi <- function(td, trait_names, sigma2, tau2, nsamples, node_index) {

  phy <- td@phylo
  # if(class(phy) != "phylo") {stop(td," does not have a \" phylo \".")}

  ntaxa <- length(phy$tip.label)
  nchr <- length(trait_names)
  characters <- as.matrix(td@data[trait_names][1:ntaxa,])
  descendants <- phy$edge[phy$edge[,1] == node_index,][,2]

  if (descendants[1] <= ntaxa) { # if left descendant is a tip
    v_left <- sigma2 * phy$edge.length[phy$edge[,2] == descendants[1]] +
      diag(tau2[descendants[1],], names = T)
    # add extra variance
    # do we need to think about covariance here? Or is it kind of redundancy?
    chr_left <- characters[descendants[1],]
    # likelihood of sample mean
    loglik_left <- sum(dnorm(characters[descendants[1],],
                         mean = characters[descendants[1],],
                         sd = tau2[descendants[1],]/nsamples[descendants[1],],
                         log = T))
    # loglik_left <- sum(dnorm(sample, mean, sd ,log = T))
  }
  else { # if left offspring node is not a tip
    recursive <- pruning_multi(td, trait_names, sigma2, tau2, nsamples, descendants[1])

    v_left <- sigma2 * phy$edge.length[phy$edge[,2] == descendants[1]]+
      recursive$variance
    chr_left <- recursive$chr.values
    loglik_left <- recursive$loglik
  }

  if (descendants[2] <= ntaxa) { # if right offspring node is a tip
    v_right <- sigma2 * phy$edge.length[phy$edge[,2] == descendants[2]] +
      diag(tau2[descendants[2],], names = T)
    # same for this epsilon
    chr_right <- characters[descendants[2],]
    # likelihood of sample mean
    loglik_right <- sum(dnorm(characters[descendants[2],],
                             mean = characters[descendants[2],],
                             sd = tau2[descendants[2],]/nsamples[descendants[2],],
                             log = T))
    # loglik_left <- sum(dnorm(sample, mean, sd ,log = T))
  }
  else { # right offspring note is not a tip
    recursive <- pruning_multi(td, trait_names, sigma2, tau2, nsamples, descendants[2])

    v_right <- sigma2 * phy$edge.length[phy$edge[,2] == descendants[2]]+
      recursive$variance
    chr_right <- recursive$chr.values
    loglik_right <- recursive$loglik
  }

  #extended branch length and variance
  v_node <- (v_left * v_right) / (v_left + v_right)

  # character value for the node is the weighted average of two descendants
  chr.values <- (v_left*chr_right + v_right*chr_left) / (v_left + v_right)

  # log-likellihood = sum of left & right descendants' likelihood and the likelihood of itself
  log.likelihood <- loglik_left + loglik_right
  log.likelihood <- log.likelihood - 1/2 * nchr * log(2*pi) -
    1/2 * log(det(v_right + v_left)) -
    1/2 * t(chr_left - chr_right) %*% solve(v_left + v_right) %*% (chr_left - chr_right)

  # log.likelihood <- log.likelihood - 1/2 * log(2*pi*sigma2*(edge_left + edge_right)) -
  #   1/2 * (chr_left - chr_right)^2 / (sigma2 * (edge_left + edge_right))


  return(list(variance = v_node,
              chr.values = chr.values,
              loglik = log.likelihood))
}



loglik_BM_multi <- function(td, trait_names, mu, sigma2, tau2, nsamples) {

  phy <- td@phylo
  if(class(phy) != "phylo") {stop(td," does not have a \" phylo \".")}

  ntaxa <- length(phy$tip.label)
  nchr <- length(trait_names)
  # characters <- td@data[[trait_names]][1:ntaxa]
  root_index <- ntaxa + 1

  root <- pruning_multi(td, trait_names, sigma2, tau2, nsamples, root_index)
  log.likelihood <- root$loglik - 1/2 * nchr * log(2*pi) -
    1/2 * log(det(root$variance)) -
    1/2 * t(root$chr.values - mu) %*% solve(root$variance) %*% (root$chr.values - mu)

  return(log.likelihood)
}
