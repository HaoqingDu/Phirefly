# Author: Haoqing Du
# Latest Editing Time: 26/04/2024

## function is to calculate the likelihood:
## - with the Ornstein-Uhlenbeck Process
## - under the situation that (i) single trait; (ii) multi-samples per taxa.
## - using pruning methods

pruning_OU <- function(td, trait_names, sigma2, alpha, theta, node_index) {
  phy <- td@phylo
  # if(class(phy) != "phylo") {stop(td," does not have a \" phylo \".")}

  ntaxa <- length(phy$tip.label)
  nchr <- length(trait_names)
  characters <- as.matrix(td@data[trait_names][1:ntaxa,])
  descendants <- phy$edge[phy$edge[,1] == node_index,][,2]

  if (descendants[1] <= ntaxa) { # if left descendant is a tip
    edge_left <- phy$edge.length[phy$edge[,2] == descendants[1]]
    v_left <- sigma2 / (2 * alpha) * expm1(2 * alpha * edge_left)
    chr_left <- characters[descendants[1],]
    mean_left <- exp(alpha * edge_left) * (chr_left - theta) + theta
    loglik_left <- 0
  }
  else { # if left offspring node is not a tip
    recursive <- pruning_OU(td, trait_names, sigma2, alpha, theta, descendants[1])

    edge_left <- phy$edge.length[phy$edge[,2] == descendants[1]]
    v_left <- sigma2 / (2 * alpha) * expm1(2 * alpha * edge_left) +
      recursive$Var * exp(2 * alpha * edge_left)

    chr_left <- recursive$chr.values
    mean_left <- exp(alpha * edge_left) * (chr_left - theta) + theta
    loglik_left <- recursive$loglik
  }

  if (descendants[2] <= ntaxa) { # if right offspring node is a tip
    edge_right <- phy$edge.length[phy$edge[,2] == descendants[2]]
    v_right <- sigma2 / (2 * alpha) * expm1(2 * alpha * edge_right)
    chr_right <- characters[descendants[2],]
    mean_right <- exp(alpha * edge_right) * (chr_right - theta) + theta
    loglik_right <- 0
  }
  else { # right offspring note is not a tip
    recursive <- pruning_OU(td, trait_names, sigma2, alpha, theta, descendants[2])

    edge_right <- phy$edge.length[phy$edge[,2] == descendants[2]]
    v_right <- sigma2 / (2 * alpha) * expm1(2 * alpha * edge_right) +
      recursive$Var * exp(2 * alpha * edge_right)

    chr_right <- recursive$chr.values
    mean_right <- exp(alpha * edge_right) * (chr_right - theta) + theta
    loglik_right <- recursive$loglik
  }

  # Variance of the node
  var_node <- v_left * v_right / (v_left + v_right)

  # character values of the node
  chr_node <- (mean_left * v_right + mean_right * v_left) / (v_left + v_right)
  #print(chr_node)

  # log-likelihood
  log.likelihood <- loglik_left + loglik_right
  log.likelihood <- log.likelihood + alpha * (edge_left + edge_right)
  log.likelihood <- log.likelihood - 1/2 * log(2*pi) -
    1/2 * log(v_left + v_right) -
    1/2 * (mean_left - mean_right)^2 / (v_left + v_right)

  # output
  return(list(Var = var_node,
              chr.values = chr_node,
              loglik = log.likelihood))
}


loglik_OU_prun <- function(td, trait_names, sigma2, alpha, theta) {

  phy <- td@phylo
  if(class(phy) != "phylo") {stop(td," does not have a \" phylo \".")}

  ntaxa <- length(phy$tip.label)
  nchr <- length(trait_names)
  # characters <- td@data[[trait_names]][1:ntaxa]
  root_index <- ntaxa + 1

  root <- pruning_OU(td, trait_names, sigma2, alpha, theta, root_index)

  log.likelihood <- root$loglik - 1/2 * log(2*pi) -
    1/2 * log(root$Var) -
    1/2 * (root$chr.values - theta)^2 / (root$Var)

  return(log.likelihood)
}
