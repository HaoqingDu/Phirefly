# Author: Haoqing Du
# Latest Editing Time: 15/03/2024

## Most recent common ancestor of two taxa

common.ancestor <- function(phy, tip1, tip2){
  if (sum(phy$edge[,2]==tip1) == 0 | sum(phy$edge[,2]==tip2) == 0) {
    stop(tip1, "\ or\ ", tip2, ": not a valid tip label in this phylo")
    }
  if (phy$edge[phy$edge[,2]==tip1,1]
      == phy$edge[phy$edge[,2]==tip2,1]) {
    return(phy$edge[phy$edge[,2]==tip1,1])
  }
  else if (phy$edge[phy$edge[,2]==tip1,1]>phy$edge[phy$edge[,2]==tip2,1]) {
    return(common.ancestor(phy, phy$edge[phy$edge[,2]==tip1,1], tip2))
  }
  else {
    return(common.ancestor(phy, tip1, phy$edge[phy$edge[,2]==tip2,1]))
  }
  }



## The branch length from the root to this node
## If the node is the mrca of two taxa, the result is their shared branch length


from.root <- function(phy, node){
  if (node == length(phy$tip.label)+1) {
    return(0)
  } else{
    return(phy$edge.length[phy$edge[,2]==node]+from.root(phy$edge[phy$edge[,2]==node,1]))
  }
}
for (i in 1:length(phy$tip.label)) {
  phy$edge[phy$edge[,2]==i,1]
}


## The distance between two taxa

dist <- function(phy, taxa1, taxa2){
  mrca <- common.ancestor(phy, taxa1, taxa2)

  edgelength <- function(phy, tip, node) {
    if (phy$edge[phy$edge[,2] == tip, 1] == node){
      return(phy$edge.length[phy$edge[,2] == tip])
    } else {
      return(phy$edge.length[phy$edge[,2] == tip] +
               edgelength(phy, phy$edge[phy$edge[,2] == tip, 1], node))
    }
  }

  print(edgelength(phy, taxa1, mrca) + edgelength(phy, taxa2, mrca))
}


## Variance-Covariance matrix for a phylogeny
## Model: Brownian Motion

vcv.BM <- function(phy){
  if(class(phy) != "phylo") stop(phy," is not a \" phylo \".")
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


## vcv matrix for Ornstein-Uhlenbeck process

vcv.OU <- function(phy, sigma2, adpt.rate){
  if(class(phy) != "phylo") stop(phy," is not a \" phylo \".")
  C <- matrix(NA, nrow = length(phy$tip.label), ncol = length(phy$tip.label))
  for (n in 1:ncol(C)) {
    for (m in 1:n) {
      C[m,n] <- sigma2/(2*adpt.rate)*
        exp(-adpt.rate*dist(phy, m, n))*
        (1-exp(-2*adpt.rate*from.root(common.ancestor(m,n))))
      C[n,m] <- C[m,n]
    }
  }
  return(C)
}
