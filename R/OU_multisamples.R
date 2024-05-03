# Author: Haoqing Du
# Latest Editing Time: 03/05/2024

## function is to calculate the likelihood:
## - with the Ornstein-Uhlenbeck Process
## - under the situation that (i) multiple traits; (ii) multi-samples per taxa.
## - using pruning methods

pruning_OU_MVN <- function(td, trait_names, sigma2, alpha, theta, tau, nsamples, node_index) {
  phy <- td@phylo
  # if(class(phy) != "phylo") {stop(td," does not have a \" phylo \".")}

  ntaxa <- length(phy$tip.label)
  nchr <- length(trait_names)
  characters <- as.matrix(td@data[trait_names][1:ntaxa,])
  descendants <- phy$edge[phy$edge[,1] == node_index,][,2]

  if (descendants[1] <= ntaxa) { # if left descendant is a tip
    edge_left <- phy$edge.length[phy$edge[,2] == descendants[1]] +
      tau[descendants[1],]/nsamples[descendants[1],]
    t_left <- 1 / (2 * alpha) * expm1(2 * alpha * edge_left)
    v_left <- sigma2 * t_left
    chr_left <- characters[descendants[1],]
    mean_left <- exp(alpha * edge_left) * (chr_left - theta) + theta
    loglik_left <- 1/2 * nchr * log(2*pi) -
      1/2 * log(det(sigma2 * tau[descendants[1],]))
  }
  else { # if left offspring node is not a tip
    recursive <- pruning_OU_MVN(td, trait_names, sigma2, alpha, theta, descendants[1])

    edge_left <- phy$edge.length[phy$edge[,2] == descendants[1]]
    t_left <- 1 / (2 * alpha) * expm1(2 * alpha * edge_left) +
      recursive$extended * exp(2.0 * alpha * edge_left)
    v_left <- sigma2 * t_left

    chr_left <- recursive$chr.values
    mean_left <- exp(alpha * edge_left) * (chr_left - theta) + theta
    loglik_left <- recursive$loglik
  }

  if (descendants[2] <= ntaxa) { # if right offspring node is a tip
    edge_right <- phy$edge.length[phy$edge[,2] == descendants[2]] +
      tau[descendants[2],]/nsamples[descendants[2],]
    t_right <- 1 / (2 * alpha) * expm1(2 * alpha * edge_right)
    v_right <- sigma2 * t_right
    chr_right <- characters[descendants[2],]
    mean_right <- exp(alpha * edge_right) * (chr_right - theta) + theta
    loglik_right <- - 1/2 * nchr * log(2*pi) -
      1/2 * log(det(sigma2 * tau[descendants[2],]))
  }
  else { # right offspring note is not a tip
    recursive <- pruning_OU_MVN(td, trait_names, sigma2, alpha, theta, descendants[2])

    edge_right <- phy$edge.length[phy$edge[,2] == descendants[2]]
    t_right <- 1 / (2 * alpha) * expm1(2 * alpha * edge_right) +
      recursive$extended * exp(2.0 * alpha * edge_right)
    v_right <- sigma2 * t_right

    chr_right <- recursive$chr.values
    mean_right <- exp(alpha * edge_right) * (chr_right - theta) + theta
    loglik_right <- recursive$loglik
  }

  # Variance of the node
  t_node <- t_left * t_right / (t_left + t_right)
  var_node <- sigma2 * t_node

  # character values of the node
  chr_node <- (mean_left * t_right + mean_right * t_left) / (t_left + t_right)

  # log-likelihood
  log.likelihood <- loglik_left + loglik_right
  log.likelihood <- log.likelihood - 1/2 * nchr * log(2*pi) -
    1/2 * log(det(v_left + v_right)) -
    1/2 * t(mean_left - mean_right) %*% solve(v_left + v_right) %*% (mean_left - mean_right)

  # output
  return(list(extended = t_node,
              Var = var_node,
              chr.values = chr_node,
              loglik = log.likelihood))
}


loglik_OU_MVN <- function(td, trait_names, sigma2, alpha, theta) {

  phy <- td@phylo
  if(class(phy) != "phylo") {stop(td," does not have a \" phylo \".")}

  ntaxa <- length(phy$tip.label)
  nchr <- length(trait_names)
  # characters <- td@data[[trait_names]][1:ntaxa]
  root_index <- ntaxa + 1

  root <- pruning_OU_MVN(td, trait_names, sigma2, alpha, theta, root_index)
  log.likelihood <- root$loglik - 1/2 * nchr * log(2*pi) -
    1/2 * log(det(root$Var)) -
    1/2 * t(root$chr.values - theta) %*% solve(root$Var) %*% (root$chr.values - theta)

  return(log.likelihood)
}
