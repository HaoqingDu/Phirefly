# Author: Haoqing Du
# Latest Editing Time: 06/05/2024

## Markov Chain Monte Carlo
## The function is designed specially for the Phirefly package.

mcmc <- function(LikelihoodFunction,
                 priors, # prior probability distribution function, returns log prior probbility
                 parameters,
                 logTransforms,
                 delta, # a vector of step lengths
                 iterations,
                 burnin = round(interations/3), # reason?
                 thining = 1, # sample every k steps
                 print.info = FALSE,
                 save.file,

){
  # sample list to save parameter values
  chain <- array(dim = c(ceiling(iterations/thining), length(parameters)+2))

  # pre-compute current posterior probability
  lnl <- LikelihoodFunction(parameters) # do.call
  ln_prior <- 0
  for (i in 1:length(parameters)) {
    ln_prior <- ln_prior + priors[[i]](parameters[i])
  }
  pp <- lnl + ln_prior

  if (burnin > 0) {
    if (print.info) {
      cat("Start burn-in")
      bar <- txtProgressBar(style=1,width=42)
    }

    delta
  }
  # if (print.info) {
  #   cat("Burning-in the chain ...\n")
  #   cat("0--------25--------50--------75--------100\n")
  #   bar <- txtProgressBar(style=1,width=42)
  # }
  #
  # if (burnin == 0) {
  #   chain[1,] <- c(lnl, ln_prior, parameters)
  # }
  #
  # for (i in 1:(burnin + iterations)) {
  #   if (print.info) {
  #     if (i <= burnin) {
  #
  #     }
  #   }
  # }

}
