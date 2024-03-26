# Author: Haoqing Du
# Latest Editing Time: 26/03/2024

## function is to calculate the likelihood:
## - with the Brownian Motion model
## - under the situation that (i) single trait; (ii) single sample per taxa.
## - using pruning methods


## pruning
pruning <- function(td, trait_names, node_index) {

  phy <- td@phylo
  if(class(phy) != "phylo") {stop(td," does not have a \" phylo \".")}

  offspring_node <- tree$edge[tree$edge[,1] == node.label,][,2]
  if (sum(tree$edge[,1] ==offspring_node[1]) == 0) { # if left offspring node is a tip
    edge_left <- tree$edge.length[which(tree$edge[,1] == node.label)[1]]
    chr_left <- character.value[offspring_node[1],]
    if (sum(tree$edge[,1] ==offspring_node[2] ) == 0) { # if right offspring node is a tip
      edge_right <- tree$edge.length[which(tree$edge[,1] == node.label)[2]]
      chr_right <- character.value[offspring_node[2],]
      #footprint <<- c(footprint, node.label)
      contrast <<- rbind(contrast, chr_left-chr_right)
      edgelength <<- rbind(edgelength, c(edge_left, edge_right))
      return(list(chr.values = rbind(chr_left, chr_right), edge = c(edge_left, edge_right)))
    } else { # right offspring note is not a tip
      recursive <- pruning2(character.value, tree, offspring_node[2])
      edge_right <- tree$edge.length[which(tree$edge[,1] == node.label)[2]]+
        prod(recursive$edge)/
        sum(recursive$edge)
      x_left <- recursive$chr.values[1,]
      x_right <- recursive$chr.values[2,]
      v_left <- recursive$edge[1]
      v_right <- recursive$edge[2]
      chr_right <- (x_left*v_right+x_right*v_left)/
        (v_left+v_right)
      contrast <<- rbind(contrast, chr_left-chr_right)
      edgelength <<- rbind(edgelength, c(edge_left, edge_right))
      #footprint <<- c(footprint, node.label)
      return(list(chr.values = rbind(chr_left, chr_right), edge = c(edge_left, edge_right)))
    }
  } else { # if left offspring node is not a tip
    recursive <- pruning2(character.value, tree, offspring_node[1])
    edge_left <- tree$edge.length[which(tree$edge[,1] == node.label)[1]]+
      prod(recursive$edge)/
      sum(recursive$edge)
    x_left <- recursive$chr.values[1,]
    x_right <- recursive$chr.values[2,]
    v_left <- recursive$edge[1]
    v_right <- recursive$edge[2]
    chr_left <- (x_left*v_right+x_right*v_left)/
      (v_left+v_right)
    if (sum(tree$edge[,1] ==offspring_node[2] ) == 0) { # if right offspring node is a tip
      edge_right <- tree$edge.length[which(tree$edge[,1] == node.label)[2]]
      chr_right <- character.value[offspring_node[2],]
      contrast <<- rbind(contrast, chr_left-chr_right)
      edgelength <<- rbind(edgelength, c(edge_left, edge_right))
      #footprint <<- c(footprint, node.label)
      return(list(chr.values = rbind(chr_left, chr_right), edge = c(edge_left, edge_right)))
    } else { # right offspring note is not a tip
      recursive <- pruning2(character.value, tree, offspring_node[2])
      edge_right <- tree$edge.length[which(tree$edge[,1] == node.label)[2]]+
        prod(recursive$edge)/
        sum(recursive$edge)
      x_left <- recursive$chr.values[1,]
      x_right <- recursive$chr.values[2,]
      v_left <- recursive$edge[1]
      v_right <- recursive$edge[2]
      chr_right <- (x_left*v_right+x_right*v_left)/
        (v_left+v_right)
      contrast <<- rbind(contrast, chr_left-chr_right)
      edgelength <<- rbind(edgelength, c(edge_left, edge_right))
      #footprint <<- c(footprint, node.label)
      return(list(chr.values = rbind(chr_left, chr_right), edge = c(edge_left, edge_right)))
    }
  }
}
