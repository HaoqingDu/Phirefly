# Author: Haoqing Du
# Latest Editing Time: 15/03/2024

## Variance-Covariance matrix for a phylogeny


vcv.BM <- function(phy){
  if(class(phy) != "phylo") stop(phy," is not a \" phylo \".")
  common.ancestor <- function(i,j){
    if (phy$edge[phy$edge[,2]==i,1]
        == phy$edge[tree$edge[,2]==j,1]) {
      return(tree$edge[tree$edge[,2]==i,1])
    }
    else if (tree$edge[tree$edge[,2]==i,1]>tree$edge[tree$edge[,2]==j,1]) {
      return(common.ancestor(tree$edge[tree$edge[,2]==i,1],j))
    }
    else {
      return(common.ancestor(i,tree$edge[tree$edge[,2]==j,1]))
    }
  }

  from.root <- function(node){
    if (node == length(tree$tip.label)+1) {
      return(0)
    } else{
      return(tree$edge.length[tree$edge[,2]==node]+from.root(tree$edge[tree$edge[,2]==node,1]))
    }
  }
  for (i in 1:length(tree$tip.label)) {
    tree$edge[tree$edge[,2]==i,1]
  }

  C <- matrix(NA, nrow = length(tree$tip.label), ncol = length(tree$tip.label))

  for (n in 1:ncol(C)) {
    for (m in 1:n) {
      if (m!=n) {
        C[m,n] <- from.root(common.ancestor(m,n))
        C[n,m] <- C[m,n]
      } else {
        C[m,n] <- from.root(tree$edge[tree$edge[,2] == m,1]) + tree$edge.length[tree$edge[,2] == m]
      }
    }
  }
  return(C)
}
