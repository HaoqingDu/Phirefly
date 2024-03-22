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


# Generalised least square method (aka variance-covariance matrix method)
simple_ou_vcv <- function(tree, continuousChar, alpha, sigma2, theta){
  ntip <- length(tree$tip.label)
  mrca1 <- ape::mrca(tree) # get the ancestral node label for each pair of tips
  times <- ape::node.depth.edgelength(tree) # get time at each node from root
  ta <- matrix(times[mrca1], nrow=ntip, dimnames = list(tree$tip.label, tree$tip.label)) # get time of divergence for each pair of tips
  T.term <- times[1:ntip] # get time at tips
  tia <- times[1:ntip] - ta
  tja <- t(tia)
  tij <- tja + tia # distance in time unit between two tips

  vy = sigma2 / (2*alpha)

  #V = vy * (1 - exp(-2 * alpha * ta)) * exp(-alpha * tij)
  V = vy * -1 * expm1(-2 * alpha * ta) * exp(-alpha * tij) ### ta = time tgt; tij = time not tgt (sum of two branches)

  X = matrix(1, ntip)

  C = chol(V) # upper triangular matrix
  L = t(C) # lower triangular matrix
  log_det_V = 0
  for (i in 1:ntip){
    log_det_V = log_det_V + log(L[i,i])
  }
  log_det_V = log_det_V *2.0 # equals to julia implementation to 12 sig. fig.

  y = NULL
  for (species in tree$tip.label){
    y[species] = as.numeric(continuousChar[species])
  }

  r = solve(L) %*% y - solve(L) %*% X * theta # what does de-correlated residuals mean?

  # res = - (n/2) * log(2*pi) - 0.5 * log_det_V - 0.5 * dot(r, r)
  #     = exp(-n/2)^(2*pi) * exp(-0.5)^det_V * exp(-0.5)^dot(r, r) ?
  res = 0.0
  res = res - 0.5 * ntip * log(2*pi)
  res = res - 0.5 * log_det_V
  res = res - 0.5 * dot(r, r) # what is r and what is  dot product of r?

  return(res)
}
