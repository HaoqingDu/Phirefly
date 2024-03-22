# Pruning algorithm

# Postorder function
postorder <- function(node_index, edge, tree, continuousChar,
                      mu, V, log_norm_factor, branch_lengths, alpha, sigma2, theta){
  ntip = length(tree$tip.label)

  # if is internal node
  if (node_index > ntip){

    left_edge  = which(edge[,1] == node_index)[1] # index of left child edge
    right_edge = which(edge[,1] == node_index)[2] # index of right child edge
    left = edge[left_edge,2] # index of left child node
    right = edge[right_edge,2] # index of right child node

    output_left <- postorder(left, edge, tree, continuousChar,
                             mu, V, log_norm_factor, branch_lengths, alpha, sigma2, theta)
    mu <- output_left[[1]]
    V <- output_left[[2]]
    log_norm_factor <- output_left[[3]]

    output_right <- postorder(right, edge, tree, continuousChar,
                              mu, V, log_norm_factor, branch_lengths, alpha, sigma2, theta)
    mu <- output_right[[1]]
    V <- output_right[[2]]
    log_norm_factor <- output_right[[3]]


    bl_left = branch_lengths[left_edge] # all branch of left child edge
    bl_right = branch_lengths[right_edge] # all branch of right child edge

    # 1) variance of the normal variable: this branch (v_left) and the subtree (V[left])

    v_left = sigma2/(2*alpha) *expm1(2.0*alpha*bl_left)
    var_left = v_left + V[left] * exp(2.0 * alpha * bl_left)

    v_right = sigma2/(2*alpha) *expm1(2.0*alpha*bl_right)
    var_right = v_right + V[right] * exp(2.0 * alpha * bl_right)

    # 2) mean of the normal variable
    mean_left = exp(alpha*bl_left)*(mu[left] - theta) + theta
    mean_right = exp(alpha*bl_right)*(mu[right] - theta) + theta

    ## compute the mean and variance of the node
    mean_ancestor = (mean_left * var_right + mean_right * var_left) / (var_left + var_right)
    mu[node_index] = mean_ancestor
    var_node = (var_left * var_right) / (var_left + var_right)
    V[node_index] = var_node

    ## compute the normalizing factor, the left-hand side of the pdf of the normal variable
    log_nf_left = bl_left * alpha
    log_nf_right = bl_right * alpha

    contrast = mean_left - mean_right
    a = -(contrast*contrast / (2*(var_left+var_right)))
    b = log(2*pi*(var_left+var_right))/2.0
    log_nf = log_nf_left + log_nf_right + a - b
    log_norm_factor[node_index] = log_nf


    return(list(mu, V, log_norm_factor))
  }


  # if is tip
  else{
    species = tree$tip.label[node_index]

    mu[node_index] = as.numeric(continuousChar[[which(names(continuousChar) == species)]])
    V[node_index] = 0.0 ## if there is no observation error

    return(list(mu, V, log_norm_factor))
  }
}

## log-likelihood + root treatment
simple_ou_pruning <- function(tree, continuousChar, alpha, sigma2, theta){
  ntip = length(tree$tip.label) # number of tips
  edge = tree$edge # equals tree[:edge] in Julia
  n_edges = length(edge[,1]) # number of edges
  max_node_index = max(tree$edge) # total number of nodes

  V = numeric(max_node_index)
  mu = numeric(max_node_index)
  log_norm_factor = numeric(max_node_index)

  branch_lengths = tree$edge.length

  root_index = ntip + 1

  output <- postorder(root_index, edge, tree, continuousChar,
                      mu, V, log_norm_factor, branch_lengths, alpha, sigma2, theta)
  mu <- output[[1]]
  V <- output[[2]]
  log_norm_factor <- output[[3]]

  ## assume root value equal to theta
  mu_root = mu[root_index]
  v_root = V[root_index]
  lnl = dnorm(theta, mean = mu_root, sd = sqrt(v_root), log = TRUE) # are \theta and \mu in correct positions?

  ## add norm factor
  for (log_nf in log_norm_factor){
    lnl = lnl + log_nf
  }
  return(lnl)
}


# multiple continuous characters
# continuousChar. alpha, sigma2, and theta should be vectors of same length
#multiChar.ou_pruning <- function(tree, continuousChar, alpha, sigma2, theta){
#  sum_likelihood = 0
#  for (i in 1:length(continuousChar)){
#    sum_likelihood = sum_likelihood + simple_ou_pruning(tree, continuousChar[i], alpha[i], sigma2[i], theta[i])
#  }
#  return(sum_likelihood)
#}


#multiSample.postorder <- function(node_index, edge, tree, continuousChar,
#                      mu, V, log_norm_factor, branch_lengths, alpha, sigma2, theta){
#  ntip = length(tree$tip.label)
#
#  # if is internal node
#  if (node_index > ntip){
#
#    left_edge  = which(edge[,1] == node_index)[1] # index of left child edge
#    right_edge = which(edge[,1] == node_index)[2] # index of right child edge
#    left = edge[left_edge,2] # index of left child node
#    right = edge[right_edge,2] # index of right child node
#
#    output_left <- postorder(left, edge, tree, continuousChar,
#                             mu, V, log_norm_factor, branch_lengths, alpha, sigma2, theta)
#    mu <- output_left[[1]]
#    V <- output_left[[2]]
#    log_norm_factor <- output_left[[3]]
#
#    output_right <- postorder(right, edge, tree, continuousChar,
#                              mu, V, log_norm_factor, branch_lengths, alpha, sigma2, theta)
#    mu <- output_right[[1]]
#    V <- output_right[[2]]
#    log_norm_factor <- output_right[[3]]
#
#
#    bl_left = branch_lengths[left_edge] # all branch of left child edge
#    bl_right = branch_lengths[right_edge] # all branch of right child edge
#
#    # 1) variance of the normal variable: this branch (v_left) and the subtree (V[left])
#
#    v_left = sigma2/(2*alpha) *expm1(2.0*alpha*bl_left)
#    var_left = v_left + V[left] * exp(2.0 * alpha * bl_left)
#
#    v_right = sigma2/(2*alpha) *expm1(2.0*alpha*bl_right)
#    var_right = v_right + V[right] * exp(2.0 * alpha * bl_right)
#
#    # 2) mean of the normal variable
#    mean_left = exp(alpha*bl_left)*(mu[left] - theta) + theta
#    mean_right = exp(alpha*bl_right)*(mu[right] - theta) + theta
#
#    ## compute the mean and variance of the node
#    mean_ancestor = (mean_left * var_right + mean_right * var_left) / (var_left + var_right)
#    mu[node_index] = mean_ancestor
#    var_node = (var_left * var_right) / (var_left + var_right)
#    V[node_index] = var_node
#
#    ## compute the normalizing factor, the left-hand side of the pdf of the normal variable
#    log_nf_left = bl_left * alpha
#    log_nf_right = bl_right * alpha
#
#    contrast = mean_left - mean_right
#    a = -(contrast*contrast / (2*(var_left+var_right)))
#    b = log(2*pi*(var_left+var_right))/2.0
#    #b = log(2*pi)/2.0 + log(var_left+var_right)/2.0
#    log_nf = log_nf_left + log_nf_right + a - b
#    log_norm_factor[node_index] = log_nf
#
#
#    return(list(mu, V, log_norm_factor))
#  }
#
#
#  # if is tip
#  else{
#    species = tree$tip.label[node_index]
#
#    mu[node_index] = as.numeric(continuousChar[[which(names(continuousChar) == species)]])
#    V[node_index] = 0.0 ## if there is no observation error
#
#    return(list(mu, V, log_norm_factor))
#  }
#}
#multiSample.ou_pruning <-
