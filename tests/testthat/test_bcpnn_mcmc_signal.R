context("bcpnn_mcmc_signal() checks")




test_that("bcpnn_mcmc_signal() correctly throws errors for bad input", {
  
  ### zeros ok for BCPNN method
  # expect_warning(
  #   bcpnn_mcmc_signal(0,2,1,0),
  #   "Zero counts detected"
  # )
  
  expect_error(
    bcpnn_mcmc_signal(-1,2,3,4),
    "Supplied values must be non-negative"
  )
  expect_error(
    bcpnn_mcmc_signal(1,2,3,NA),
    "Supplied values must not be or contain `NA` or infinite values"
  )
  expect_error(
    bcpnn_mcmc_signal(1,Inf,3,7),
    "Supplied values must not be or contain `NA` or infinite values"
  )
  expect_error(
    bcpnn_mcmc_signal(1,10,NULL,7),
    "Supplied values must not be `NULL`"
  )
  
})


test_that("bcpnn_mcmc_signal() correctly calculates known values", {
  
  
  
  # These examples are taken from:
  # see Noren et al. (2006). Stats in Med, 25:3740–3757. (Table III, page 3749)
  
  n11 <- 25
  n1. <- 1126
  n.1 <- 87
  n.. <- 572573
  set.seed(123456)
  stratum1 <- bcpnn_mcmc_signal(n11, n1. - n11, n.1 - n11, n.. - n11 - (n1. - n11) - (n.1 - n11))
  # stratum1$est - 5.25 # tab III values (will be slightly different because of sampling)
  # stratum1$ci_lo - 4.64 # tab III values (will be slightly different because of sampling)
  
  expect_equal(
    stratum1$est,
    5.247846,
    tol = 1e-6
  )
  expect_equal(
    stratum1$ci_lo,
    4.649270,
    tol = 1e-6
  )
  
  
  n11 <- 203
  n1. <- 30068
  n.1 <- 508
  n.. <- 155209
  set.seed(123456)
  stratum3 <- bcpnn_mcmc_signal(n11, n1. - n11, n.1 - n11, n.. - n11 - (n1. - n11) - (n.1 - n11))
  # stratum3$est - 1.04 # tab III values (will be slightly different because of sampling)
  # stratum3$ci_lo - (0.87) # tab III values (will be slightly different because of sampling)
  
  expect_equal(
    stratum3$est,
    1.0408,
    tol = 1e-6
  )
  expect_equal(
    stratum3$ci_lo,
    0.8803391,
    tol = 1e-6
  )
  
  n11 <- 0
  n1. <- 5232
  n.1 <- 3
  n.. <- 80140
  set.seed(123456)
  stratum4 <- bcpnn_mcmc_signal(n11, n1. - n11, n.1 - n11, n.. - n11 - (n1. - n11) - (n.1 - n11))
  # stratum4$est - (-0.48) # tab III values (will be slightly different because of sampling)
  # stratum4$ci_lo - (-11.10) # tab III values (will be slightly different because of sampling)
  
  expect_equal(
    stratum4$est,
    -0.4768594,
    tol = 1e-6
  )
  expect_equal(
    stratum4$ci_lo,
    -10.3099,
    tol = 1e-6
  )
  
  
  
  n11 <- 0
  n1. <- 10
  n.1 <- 0
  n.. <- 453481
  set.seed(123456)
  stratum7 <- bcpnn_mcmc_signal(n11, n1. - n11, n.1 - n11, n.. - n11 - (n1. - n11) - (n.1 - n11))
  # stratum7$est - 0.00 # tab III values (will be slightly different because of sampling)
  # stratum7$ci_lo - (-10.66) # tab III values (will be slightly different because of sampling)
  
  expect_equal(
    stratum7$est,
    1.590577e-06,
    tol = 1e-6
  )
  expect_equal(
    stratum7$ci_lo,
    -9.90467,
    tol = 1e-6
  )
  
  
  
  
})



test_that("bcpnn_mcmc_signal() correctly calculates known __vectors__ of values", {
  
  
  
  # These examples are taken from:
  # see Noren et al. (2006). Stats in Med, 25:3740–3757. (Table III, page 3749)
  
  n11 <- 25
  n1. <- 1126
  n.1 <- 87
  n.. <- 572573
  set.seed(123456)
  stratum1 <- bcpnn_mcmc_signal(n11, n1. - n11, n.1 - n11, n.. - n11 - (n1. - n11) - (n.1 - n11))
  
  n11 <- 203
  n1. <- 30068
  n.1 <- 508
  n.. <- 155209
  set.seed(123456)
  stratum3 <- bcpnn_mcmc_signal(n11, n1. - n11, n.1 - n11, n.. - n11 - (n1. - n11) - (n.1 - n11))
  
  
  n11 <- c(25, 203)
  n1. <- c(1126, 30068)
  n.1 <- c(87, 508)
  n.. <- c(572573, 155209)
  set.seed(123456)
  stratum_1_and_3 <- bcpnn_mcmc_signal(n11, n1. - n11, n.1 - n11, n.. - n11 - (n1. - n11) - (n.1 - n11))
  
  
  # stratum1
  # stratum3
  # stratum_1_and_3
  
  expect_equal(
    rbind(stratum1, stratum3)$est,
    stratum_1_and_3$est,
    tolerance = 0.001, # won't be the same as different seeds
    scale = rbind(stratum1, stratum3)$est
  )
  expect_equal(
    rbind(stratum1, stratum3)$ci_lo,
    stratum_1_and_3$ci_lo,
    tolerance = 0.001, 
    scale = rbind(stratum1, stratum3)$ci_lo
  )
  expect_equal(
    rbind(stratum1, stratum3)$ci_hi,
    stratum_1_and_3$ci_hi,
    tolerance = 0.001, 
    scale = rbind(stratum1, stratum3)$ci_hi
  )

  
  
})




