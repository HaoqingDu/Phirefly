####################################
##
##  the pruning likelihood
##
####################################

data("artiodactyla")

sigma2 <- 0.5
alpha <- 0.5
theta <- 5.0
pruning_likelihood <- logl_OU_fitzjohn(artiodactyla, alpha, sigma2, theta, "brain_mass_g_log_mean")

####################################
##
##  the vcv likelihood algorithm
##
####################################

vcv_likelihood <- logl_OU_vcv(artiodactyla, alpha, sigma2, theta, "brain_mass_g_log_mean")

####################################
##
##  the vcv likelihood algorithm
##
####################################

test_that("OU vcv and pruning likelihoods are the same", {
  expect_equal(
    vcv_likelihood,
    pruning_likelihood
  )
})
