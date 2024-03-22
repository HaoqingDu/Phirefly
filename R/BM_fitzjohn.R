#' Calculate the likelihood of a univariate BM process using the FitzJohn algorithm
#'
#' @inherit logl_OU_fitzjohn
#'
#' @param mu the parameter controlling the mean, or the ancestral value of the BM process
#' @export
#'
#' @examples
#'
#' data("artiodactyla")
#'
#' logl_BM_fitzjohn(artiodactyla, 0.1, 5.0, "brain_mass_g_log_mean")
logl_BM_fitzjohn <- function(td, sigma2, mu, trait_name){

  ntip <- length(td@phylo$tip.label)
  descendants <- make_descendants(td@phylo)

  root_index <- ntip+1
  l <- logl_BM_fitzjohn_po(td, sigma2, root_index, descendants, ntip, trait_name)


  lnl <- dnorm(l$x, mean = mu, sd = sqrt(l$v), log = TRUE)
  lnl <- lnl + l$log_norm_factor

  return(lnl)
}

logl_BM_fitzjohn_po <- function(td, sigma2, node_index, descendants, ntip, trait_name){
  is_internal <- node_index > ntip

  if(is_internal){
    desc <- descendants[[node_index]]
    left_edge_idx <- desc[[1]]
    right_edge_idx <- desc[[2]]

    left_node_idx <- td@phylo$edge[left_edge_idx,2]
    right_node_idx <- td@phylo$edge[right_edge_idx,2]

    left_bl <- td@phylo$edge.length[left_edge_idx]
    right_bl <- td@phylo$edge.length[right_edge_idx]

    left <- logl_BM_fitzjohn_po(td, sigma2, left_node_idx, descendants, ntip, trait_name)
    right <- logl_BM_fitzjohn_po(td, sigma2, right_node_idx, descendants, ntip, trait_name)

    x_left <- left$x
    x_right <- right$x

    v_left <- left$v + sigma2*left_bl
    v_right <- right$v + sigma2*right_bl

    log_norm_factor_left <- left$log_norm_factor
    log_norm_factor_right <- right$log_norm_factor

    x_node <- (v_left*x_right + v_right*x_left)/(v_left+v_right)
    v_node <- v_left*v_right/(v_left+v_right)

    contrast <- x_left-x_right
    ## norm_factor <- exp(-(x_left-x_right)^2/(2*(v_left+v_right))) / (sqrt(2*pi*(v_left+v_right)))
    log_norm_factor <- log_norm_factor_left + log_norm_factor_right
    log_norm_factor <- log_norm_factor -(x_left-x_right)^2/(2*(v_left+v_right))
    log_norm_factor <- log_norm_factor - 0.5*log(2*pi*(v_left+v_right))

  }else{
    x_node <- td@data[[trait_name]][node_index] ## does not work unless td@data is sorted by node
    v_node <- 0.0
    log_norm_factor <- 0.0
  }

  l <- list(
    "x" = x_node,
    "v" = v_node,
    "log_norm_factor" = log_norm_factor
    )
  return(l)
}
