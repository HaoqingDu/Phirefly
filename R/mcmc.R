# Author: Haoqing Du
# Latest Editing Time: 06/05/2024

## Markov Chain Monte Carlo
## The function is designed specially for the Phirefly package.

mcmc.chain <- function(LikelihoodFunction,
                       priors, # prior probability distribution function, returns log prior probbility
                       parameters,
                       logTransforms,
                       delta # a vector of step lengths
                       ) {

  # pre-compute current posterior probability
  lnl <- do.all(LikelihoodFunction(parameters)) # do.call
  ln_prior <- 0
  for (j in 1:length(parameters)) {
    ln_prior <- ln_prior + priors[[j]](parameters[j])
  }
  pp <- lnl + ln_prior

  for (j in sample(1:length(parameters))) { # randomize the order?
    if (logTransforms[j]) {
      if (parameters[j] == 0) {
        stop("cannot log-transform 0")
      }
      logvalue <- log(parameters[j]) # propose a new value for parameter[j]
      logvalue_new <- logvalue + rnorm(1, 0, delta[j])
      new_value <- exp(logvalue_new)
      parameters[j] <- new_value

      new_prior     <- 0.0
      for ( k in 1:length(parameters) ) {
        new_prior <- new_prior + priors[[k]](parameters[k])
      }
      if ( is.finite(new_prior) ) {
        new_lnl <- do.all(LikelihoodFunction(parameters))
        new_pp <- new_lnl + new_prior
      }

      # accept / reject
      if ( is.finite(new_pp) && is.finite(hr) &&  new_pp-pp > log(runif(1,0,1)) ) {
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
      if ( is.finite(new_pp) && is.finite(hr) &&  new_pp-pp > log(runif(1,0,1)) ) {
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
  # sample list to save parameter values
  chain <- array(dim = c(ceiling(iterations/thining), length(parameters)+2))



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


  }
}
