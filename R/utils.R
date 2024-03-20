#' Title
#'
#' @param an object of class "phylo"
#'
#' @return a list
#' @export
#'
#' @examples
#'
#' data("artiodactyla")
#' descendants <- make_descendants(artiodactyla@phylo)
make_descendants <- function(phy){
  ## iterate over internal nodes

  ntips <- length(phy$tip.label)
  root_index <- ntips+1

  max_node_index <- max(phy$edge)

  internal_node_indices <- root_index:max_node_index

  descendants <- list()

  for (internal_node in internal_node_indices){
    edge_indices <- which(phy$edge[,1] == internal_node)

    #child_node_indices <- phy$edge[edge_indices,2]
    descendants[[internal_node]] <- edge_indices
  }
  return(descendants)
}
