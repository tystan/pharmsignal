context("bcpnn_norm_signal() checks")




test_that("bcpnn_norm_signal() correctly throws errors for bad input", {
  
  ### zeros ok for BCPNN method
  # expect_warning(
  #   bcpnn_norm_signal(0,2,1,0),
  #   "Zero counts detected"
  # )
  
  expect_error(
    bcpnn_norm_signal(-1,2,3,4),
    "Supplied values must be non-negative"
  )
  expect_error(
    bcpnn_norm_signal(1,2,3,NA),
    "Supplied values must not be or contain `NA` or infinite values"
  )
  expect_error(
    bcpnn_norm_signal(1,Inf,3,7),
    "Supplied values must not be or contain `NA` or infinite values"
  )
  expect_error(
    bcpnn_norm_signal(1,10,NULL,7),
    "Supplied values must not be `NULL`"
  )
  
})




test_that("bcpnn_norm_signal() correctly calculates known values", {
  
  
  # These examples are taken from:
  # see Gould (2003). Pharmacoepidemiology and Drug Safety, 12: 559–574. (Table 1, page 562)
  #
  # NOTE: (#1) our implementation is the 'exact WHO approach' using the posterior distribution digamma and trigamma
  #   equations (A3) and (A4) on page 572 of the appendix
  # NOTE: (#2) our implementation uses slightly different a priori values than Gould (2003), so the values are 
  #   slightly different than those in the table 1 (but close)
  
  
  # Headache

  n11 <- 1614
  n1. <- 85304
  n.1 <- 71209
  n.. <- 4640498
  
  headache <- bcpnn_norm_signal(n11, n1. - n11, n.1 - n11, n.. - n11 - (n1. - n11) - (n.1 - n11))
  
  expect_equal(
    headache$est,
    0.301503,
    tol = 1e-6
  )
  expect_equal(
    headache$var_est,
    0.001341455,
    tol = 1e-6
  )
  expect_equal(
    headache$sd_est,
    0.03662588,
    tol = 1e-6
  )
  expect_equal(
    headache$ci_lo,
    0.2297176,
    tol = 1e-6
  )
  expect_equal(
    headache$ci_hi,
    0.3732884,
    tol = 1e-6
  )
  
  
  
  bcpnn_norm_signal(328, 85304 - 328, 3001 - 328, 4640498 - 328 - (85304 - 328) - (3001 - 328))
  
  
  # Akathisia
  
  n11 <- 328
  n1. <- 85304
  n.1 <- 3001
  n.. <- 4640498

  akathisia <- bcpnn_norm_signal(n11, n1. - n11, n.1 - n11, n.. - n11 - (n1. - n11) - (n.1 - n11))
  
  expect_equal(
    akathisia$est,
    2.547888,
    tol = 1e-6
  )
  expect_equal(
    akathisia$var_est,
    0.007052478,
    tol = 1e-6
  )
  expect_equal(
    akathisia$sd_est,
    0.08397904,
    tol = 1e-6
  )
  expect_equal(
    akathisia$ci_lo,
    2.383292,
    tol = 1e-6
  )
  expect_equal(
    akathisia$ci_hi,
    2.712484,
    tol = 1e-6
  )
  
  
  
  
})




test_that("bcpnn_norm_signal() correctly calculates known __vectors__ of values", {
  
  
  
  # These examples are taken from:
  # see Noren et al. (2006). Stats in Med, 25:3740–3757. (Table III, page 3749)
  
  n11 <- 25
  n1. <- 1126
  n.1 <- 87
  n.. <- 572573
  set.seed(123456)
  stratum1 <- bcpnn_norm_signal(n11, n1. - n11, n.1 - n11, n.. - n11 - (n1. - n11) - (n.1 - n11))
  
  n11 <- 203
  n1. <- 30068
  n.1 <- 508
  n.. <- 155209
  set.seed(123456)
  stratum3 <- bcpnn_norm_signal(n11, n1. - n11, n.1 - n11, n.. - n11 - (n1. - n11) - (n.1 - n11))
  
  
  n11 <- c(25, 203)
  n1. <- c(1126, 30068)
  n.1 <- c(87, 508)
  n.. <- c(572573, 155209)
  set.seed(123456)
  stratum_1_and_3 <- bcpnn_norm_signal(n11, n1. - n11, n.1 - n11, n.. - n11 - (n1. - n11) - (n.1 - n11))
  
  
  # stratum1
  # stratum3
  # stratum_1_and_3
  
  expect_equal(
    rbind(stratum1, stratum3)[, -1],
    stratum_1_and_3[, -1]
  )

  
  
  
})











