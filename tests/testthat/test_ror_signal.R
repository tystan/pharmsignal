context("ror_signal() checks")




test_that("ror_signal() correctly throws errors for bad input", {
  
  expect_warning(
    ror_signal(0,2,1,0),
    "Zero counts detected"
  )
  expect_error(
    ror_signal(-1,2,3,4),
    "Supplied values must be non-negative"
  )
  expect_error(
    ror_signal(1,2,3,NA),
    "Supplied values must not be or contain `NA` or infinite values"
  )
  expect_error(
    ror_signal(1,Inf,3,7),
    "Supplied values must not be or contain `NA` or infinite values"
  )
  expect_error(
    ror_signal(1,10,NULL,7),
    "Supplied values must not be `NULL`"
  )
  
})




test_that("ror_signal() correctly calculates known values", {
  
  # Singhal et al. p409 table 16. Int J Pharm Pharm Sci, Vol 7, Issue 6, 405-411
  Singhal1 <- ror_signal(28, 942, 17, 31435)
  Singhal2 <- ror_signal(122, 1320, 381, 31341)
  Singhal3 <- ror_signal(c(28, 122), c(942, 1320), c(17, 381), c(31435, 31341))
  
  expect_equal(
    rbind(Singhal1, Singhal2)[, -1], # get rid of "analysis" column
    Singhal3[, -1]
  )
  
  # manual calcs
  ror_est <- (28 / 942) / (17 / 31435)
  ror_se <- sqrt(1 / 28 + 1 / 942 + 1 / 17 + 1 / 31435)
  # function version
  ror_df <- ror_signal(28, 942, 17, 31435)
  
  expect_equal(
    ror_df$est,
    ror_est
  )
  expect_equal(
    ror_df$ci_lo,
    exp(log(ror_est) - qnorm(0.975) * ror_se)
  )
  expect_equal(
    ror_df$ci_hi,
    exp(log(ror_est) + qnorm(0.975) * ror_se)
  )
  
  expect_equal(
    ror_df$sd_est,
    ror_se
  )
  
})




