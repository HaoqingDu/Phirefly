# Author: Haoqing Du
# Latest Editing Time: 25/04/2024

## function is to calculate the likelihood:
## - with the Ornstein-Uhlenbeck Process
## - under the situation that (i) multiple independent traits; (ii) multi-samples per taxa.
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
    v_left <- sigma2 * edge_left
    chr_left <- characters[descendants[1],]
    loglik_left <- 0
  }
}
