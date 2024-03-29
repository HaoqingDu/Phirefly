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

#' Calculate the likelihood of an univariate OU process using the vcv method
#'
#' @inherit logl_OU_fitzjohn
#'
#' @export
#'
#' @examples
#'
#' data("artiodactyla")
#'
#' logl_OU_vcv(artiodactyla, 0.5, 0.5, 5.04, "brain_mass_g_log_mean")
logl_OU_vcv <- function(td, alpha, sigma2, theta, trait_name){
  tree <- td@phylo
  ntip <- length(tree$tip.label)
  continuous_trait <- td@data[[trait_name]][1:ntip]

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

  beta1 <- matrix(theta, 1,1)

  logl <- generalized_least_squares(V, X, continuous_trait, beta1)

  return(logl)
}
