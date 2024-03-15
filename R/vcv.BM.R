# Author: Haoqing Du
# Latest Editing Time: 15/03/2024

## Variance-Covariance matrix for a phylogeny


vcv.BM <- function(phy){
  if(class(phy) != "phylo") stop(phy," is not a \" phylo \".")
  common.ancestor <- function(i,j){ # function to find CA of taxa i & j
    if (phy$edge[phy$edge[,2]==i,1]
        == phy$edge[phy$edge[,2]==j,1]) {
      return(phy$edge[phy$edge[,2]==i,1])
    }
    else if (phy$edge[phy$edge[,2]==i,1]>phy$edge[phy$edge[,2]==j,1]) {
      return(common.ancestor(phy$edge[phy$edge[,2]==i,1],j))
    }
    else {
      return(common.ancestor(i,phy$edge[phy$edge[,2]==j,1]))
    }
  }

  from.root <- function(node){ # the total branch length from the root to certain node
    if (node == length(phy$tip.label)+1) {
      return(0)
    } else{
      return(phy$edge.length[phy$edge[,2]==node]+from.root(phy$edge[phy$edge[,2]==node,1]))
    }
  }
  for (i in 1:length(phy$tip.label)) {
    phy$edge[phy$edge[,2]==i,1]
  }

  C <- matrix(NA, nrow = length(phy$tip.label), ncol = length(phy$tip.label))

  for (n in 1:ncol(C)) {
    for (m in 1:n) {
      if (m!=n) {
        C[m,n] <- from.root(common.ancestor(m,n))
        C[n,m] <- C[m,n]
      } else {
        C[m,n] <- from.root(phy$edge[phy$edge[,2] == m,1]) + phy$edge.length[phy$edge[,2] == m]
      }
    }
  }
  return(C)
}
