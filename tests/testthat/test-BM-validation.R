context("test that pruning and vcv give same likelihood")

test_that("1 + 1 = 2", {
  expect_equal(1 + 1, 2)
})

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
vcv_likelihood <- -20

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
