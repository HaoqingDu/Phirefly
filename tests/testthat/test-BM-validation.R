####################################
##
##  the pruning likelihood
##
####################################

data("artiodactyla")

sigma2 <- 0.1
mu <- 5.0

pruning_likelihood <- logl_BM_fitzjohn(artiodactyla, sigma2, mu, "brain_mass_g_log_mean")

####################################
##
##  the vcv likelihood algorithm
##
####################################

## vcv_likelihood <- vcv.loglik(artiodactyla, sigma2, mu, trait_name)
vcv_likelihood <- Loglik.BM(artiodactyla@phylo, artiodactyla@data$brain_mass_g_log_mean[1:43] , mu, sigma2)

####################################
##
##  the vcv likelihood algorithm
##
####################################

test_that("BM vcv and pruning likelihoods are the same", {
  expect_equal(
    vcv_likelihood,
    pruning_likelihood
  )
})
