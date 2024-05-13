# Author: Haoqing Du
# Latest Editing Time: 06/05/2024

## Markov Chain Monte Carlo
## The function is designed specially for the Phirefly package.

# tools for conversion between a var-cov matrix and a vector
from.matrix.to.vector <- function(mat){
  n <- ncol(mat)

  vec <- mat[lower.tri(mat, diag = T)]

  return(vec)
}

from.vector.to.matrix <- function(vec) {
  n <- (sqrt(1 + 8 * length(vec)) - 1) / 2

  if ((n %% 1) != 0) {
    stop("The vector cannot be converted into a var-cov matrix.")
  }

  mat <- matrix(0, n, n)
  mat[lower.tri(mat, diag = T)] <- vec
  mat <- t(mat) + mat - diag(diag(mat))

  return(mat)
}

# mcmc
mcmc.chain <- function(LikelihoodFunction,
                       td,
                       trait_names,
                       priors, # prior probability distribution function, returns log prior probbility
                       parameters, # a list containing all the parameters
                       logTransforms,
                       delta # a vector of step lengths
                       ) {


  for (j in sample(1:length(para_vec))) { # randomize the order?
    if (logTransforms[j]) {
      if (para_vec[j] == 0) {
        stop("cannot log-transform 0")
      }
      logvalue <- log(parameters[j]) # propose a new value for parameter[j]
      logvalue_new <- logvalue + rnorm(1, 0, delta[j])
      new_value <- exp(logvalue_new)
      para_vec[j] <- new_value

      new_prior     <- 0.0
      for ( k in 1:length(parameters) ) {
        new_prior <- new_prior + priors[[k]](parameters[k])
      }
      if ( is.finite(new_prior) ) {
        new_lnl <- do.all(LikelihoodFunction(parameters))
        new_pp <- new_lnl + new_prior
      }

      # accept / reject
      if ( is.finite(new_pp) &&  new_pp-pp > log(runif(1,0,1)) ) {
        pp <- new_pp
        lnl <- new_lnl
        ln_prior <- new_prior
      }
    } else { # logTransformes = F
      value <- parameters[j] # propose a new value for parameter[j]
      new_value <- value + rnorm(1,0,delta[j])
      parameters[j] <- new_value

      new_prior <- 0.0
      for ( k in 1:length(parameters) ) {
        new_prior <- new_prior + priors[[k]](parameters[k])
      }
      if ( is.finite(new_prior) ) {
        new_lnl <- do.all(LikelihoodFunction(parameters))
        new_pp <- new_lnl + new_prior
      }

      # accept / reject
      if ( is.finite(new_pp) &&  new_pp-pp > log(runif(1,0,1)) ) {
        pp <- new_pp
        lnl <- new_lnl
        ln_prior <- new_prior
      }
    }
  }

  # output
  return(list(Parameters = parameters,
              log.prior = ln_prior,
              log.lokelihood = lnl))
}






mcmc <- function(LikelihoodFunction,
                 priors, # prior probability distribution function, returns log prior probbility
                 parameters,
                 logTransforms,
                 delta, # a vector of step lengths
                 iterations,
                 burnin = round(interations/3), # Default: burn-in first 25%
                 thining = 1, # sample every k steps
                 print.info = FALSE,
                 save.file

){
  # preparation
  ntraits <- length(trait_names)

  if (is.matrix(parameters$sigma2)) {
    parameters$sigma2 <- from.matrix.to.vector(parameters$sigma2)
  }

  para_vec <- unlist(parameters)

  # sample list to save parameter values
  chain <- array(dim = c(ceiling(iterations/thining), length(parameters)+2))

  # pre-compute current posterior probability
  lnl <- do.call(LikelihoodFunction(td, trait_names, parameters)) # do.call
  ln_prior <- 0
  for (j in 1:length(para_vec)) {
    ln_prior <- ln_prior + priors[[j]](para_vec[j])
  }
  pp <- lnl + ln_prior



  ## burn-in
  if (burnin > 0) {
    if (print.info) {
      cat("Start burn-in")
      cat("0--------25--------50--------75--------100\n")
      bar <- txtProgressBar(style=1,width=42)
    }
  } else { # burnin < 0
    chain[1,] <- c(lnl, ln_prior, parameters)
  }

  for (i in 1:(burnin+iterations)) {
    if (print.info) {
      if ( i <= burnin ) {
        setTxtProgressBar(bar,i/burnin)
      } else if (i == (burnin+1) ) {
        cat("\nFinished burnin period!\n\n")
        cat("Running the chain ...\n")
        cat("0--------25--------50--------75--------100\n")
        bar <- txtProgressBar(style=1,width=42)
      } else {
        setTxtProgressBar(bar,(i-burnin)/iterations)
      }
    }
    }



}
