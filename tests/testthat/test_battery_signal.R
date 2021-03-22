context("battery_signal() checks")




test_that("battery_signal() correctly throws errors for bad input", {
  
  expect_warning(
    battery_signal(0,2,1,0),
    "Zero counts detected"
  )
  expect_error(
    battery_signal(-1,2,3,4),
    "Supplied values must be non-negative"
  )
  expect_error(
    battery_signal(1,2,3,NA),
    "Supplied values must not be or contain `NA` or infinite values"
  )
  expect_error(
    battery_signal(1,Inf,3,7),
    "Supplied values must not be or contain `NA` or infinite values"
  )
  expect_error(
    battery_signal(1,10,NULL,7),
    "Supplied values must not be `NULL`"
  )
  
})


test_that("battery_signal() is equiv to ror, prr etc for scalar values", {
  
  
  
  # These examples are taken from:
  # see Noren et al. (2006). Stats in Med, 25:3740–3757. (Table III, page 3749)
  
  n11 <- 25
  n1. <- 1126
  n.1 <- 87
  n.. <- 572573
  
  set.seed(123456)
  str_1_bat <- battery_signal(n11, n1. - n11, n.1 - n11, n.. - n11 - (n1. - n11) - (n.1 - n11))
  
  set.seed(123456)
  str_1_ror <- ror_signal(n11, n1. - n11, n.1 - n11, n.. - n11 - (n1. - n11) - (n.1 - n11))
  str_1_prr <- prr_signal(n11, n1. - n11, n.1 - n11, n.. - n11 - (n1. - n11) - (n.1 - n11))
  str_1_bno <- bcpnn_norm_signal(n11, n1. - n11, n.1 - n11, n.. - n11 - (n1. - n11) - (n.1 - n11))
  str_1_bmc <- bcpnn_mcmc_signal(n11, n1. - n11, n.1 - n11, n.. - n11 - (n1. - n11) - (n.1 - n11))
  
  
  expect_equal(
    str_1_bat[, !(colnames(str_1_bat) %in% c("lab", "analysis"))],
    dplyr::as_tibble(dplyr::bind_rows(str_1_ror, str_1_prr, str_1_bno, str_1_bmc))[, -1]
  )

  
})



test_that("battery_signal() battery_signal() is equiv to ror, prr etc for __vectors__ of values", {
  
  
  
  # These examples are taken from:
  # see Noren et al. (2006). Stats in Med, 25:3740–3757. (Table III, page 3749)
  
  n11 <- c(25, 203)
  n1. <- c(1126, 30068)
  n.1 <- c(87, 508)
  n.. <- c(572573, 155209)
  
  set.seed(123456)
  vec_bat <- battery_signal(n11, n1. - n11, n.1 - n11, n.. - n11 - (n1. - n11) - (n.1 - n11))
  
  set.seed(123456)
  vec_ror <- ror_signal(n11, n1. - n11, n.1 - n11, n.. - n11 - (n1. - n11) - (n.1 - n11))
  vec_prr <- prr_signal(n11, n1. - n11, n.1 - n11, n.. - n11 - (n1. - n11) - (n.1 - n11))
  vec_bno <- bcpnn_norm_signal(n11, n1. - n11, n.1 - n11, n.. - n11 - (n1. - n11) - (n.1 - n11))
  vec_bmc <- bcpnn_mcmc_signal(n11, n1. - n11, n.1 - n11, n.. - n11 - (n1. - n11) - (n.1 - n11))
  
  
  expect_equal(
    vec_bat[, !(colnames(vec_bat) %in% c("lab", "analysis"))],
    dplyr::as_tibble(dplyr::bind_rows(vec_ror, vec_prr, vec_bno, vec_bmc))[, -1],
    tolerance = 1e-3 # mcmc estimates have different seeds
  )
  
  
  
  
})




