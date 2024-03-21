# Author: Haoqing Du
# Latest Editing Time: 21/03/2024

## Most recent common ancestor of two taxa

common.ancestor <- function(phy, tip1, tip2) {
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


from.root <- function(phy, node) {
  if (node == length(phy$tip.label)+1) {
    return(0)
  } else{
    return(phy$edge.length[phy$edge[,2]==node]+from.root(phy, phy$edge[phy$edge[,2]==node,1]))
  }
}


## The distance between two taxa

dist <- function(phy, taxa1, taxa2) {
  mrca <- common.ancestor(phy, taxa1, taxa2)

  edgelength <- function(phy, tip, node) {
    if (phy$edge[phy$edge[,2] == tip, 1] == node){
      return(phy$edge.length[phy$edge[,2] == tip])
    } else {
      return(phy$edge.length[phy$edge[,2] == tip] +
               edgelength(phy, phy$edge[phy$edge[,2] == tip, 1], node))
    }
  }

  return(edgelength(phy, taxa1, mrca) + edgelength(phy, taxa2, mrca))
}


## Variance-Covariance matrix for a phylogeny
## Model: Brownian Motion

vcv.BM <- function(td) {
  if(class(td@phylo) != "phylo") stop(td," does not have a \" phylo \".")
  C <- matrix(NA, nrow = length(td@phylo$tip.label), ncol = length(td@phylo$tip.label))
  for (n in 1:ncol(C)) {
    for (m in 1:n) {
      if (m!=n) {
        C[m,n] <- from.root(td@phylo, common.ancestor(phy, m, n))
        C[n,m] <- C[m,n]
      } else {
        C[m,n] <- from.root(td@phylo, td@phylo$edge[td@phylo$edge[,2] == m, 1]) +
          td@phylo$edge.length[td@phylo$edge[,2] == m]
      }
    }
  }
  return(C)
}


## vcv matrix for Ornstein-Uhlenbeck process

vcv.OU <- function(td, alpha) {
  if(class(td@phylo) != "phylo") stop(td," does not have a \" phylo \".")
  C <- matrix(NA, nrow = length(td@phylo$tip.label), ncol = length(td@phylo$tip.label))
  for (n in 1:ncol(C)) {
    for (m in 1:n) {
      # Cov[i, j] = sigma2/2alpha*exp(-alpha*tij)*[1-exp(-2alpha*tra)]
      C[m,n] <- 1/(2*alpha)*
        exp(-alpha*dist(td@phylo, m, n))*
        (1-exp(-2*alpha*from.root(td@phylo, common.ancestor(td@phylo, m, n))))
      C[n,m] <- C[m,n]
    }
  }
  return(C)
}


## merge two functions into one vcv.matrix function

vcv.matrix <- function(td,
                       model = c("BM", "OU"),
                       alpha = NULL) {

  if(class(td@phylo) != "phylo") stop(td," does not have a \" phylo \".")

  ntaxa <- length(td@phylo$tip.label)
  C <- matrix(NA, nrow = ntaxa, ncol = ntaxa)

  if(length(model) > 1) {
    warning("Please specify the model!")
    print("The given result is under the BM model")
  }

  if(model == "OU") {
    # alpha is required under OU process
    if(!is.null(alpha)) {
      for (n in 1:ntaxa) {
        for (m in 1:m) {
          # C[i,j] = C[j,i] = 1/(2alpha)*exp[-alpha*t_ij]*[1-exp[-2alpha*t_ra]]
          # t_ij = distance between i & j, t_ra = distance between root and mrca
          C[m,n] <- 1/(2*alpha)*
            exp(-alpha*dist(td@phylo, m, n))*
            (1-exp(-2*alpha*from.root(td@phylo, common.ancestor(td@phylo, m, n))))
          C[n,m] <- C[m,n]
        }
      }
    } else {
      warning("alpha is required for the OU model!")
      # when alpha is missing under the OU model, the function switch to BM automatically
      print("The given result is under the BM model.")
      vcv.matrix(td,
                 model = "BM")
    }
  }

  else {
    # when model is not specified, use BM.
    for (n in 1:ncol(C)) {
      for (m in 1:n) {
        C[m,n] <- from.root(td@phylo, common.ancestor(phy, m, n))
        C[n,m] <- C[m,n]
      }
    }
  }
}


