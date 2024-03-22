## vcv matrix for Ornstein-Uhlenbeck process

vcv.OU <- function(phy, sigma2, adpt.rate) {
  if(class(phy) != "phylo") stop(phy," is not a \" phylo \".")
  C <- matrix(NA, nrow = length(phy$tip.label), ncol = length(phy$tip.label))
  for (n in 1:ncol(C)) {
    for (m in 1:n) {
      C[m,n] <- sigma2/(2*adpt.rate)*
        exp(-adpt.rate*dist(phy, m, n))*
        (1-exp(-2*adpt.rate*from.root(phy, common.ancestor(m,n))))
      C[n,m] <- C[m,n]
    }
  }
  return(C)
}
