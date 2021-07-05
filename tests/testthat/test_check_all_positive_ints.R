context("check_all_positive_ints() checks")







test_that("check_all_positive_ints() correctly throws errors for bad input", {
  

  # these will throw errors
  expect_error(
    check_all_positive_ints(1, character(0)), 
    "Supplied values must be numeric"
  )
  expect_error(
    check_all_positive_ints(1, numeric(0)), 
    "Supplied values must have length 1 \\(scalar\\) or greater \\(vector\\)"
  )
  expect_error(
    check_all_positive_ints(-11, 2, 1:2), 
    "Supplied values must be non-negative"
  )
  expect_error(
    check_all_positive_ints(-11, NA), 
    "Supplied values must not be or contain `NA` or infinite values"
  )
  expect_error(
    check_all_positive_ints(-11, c(1, NA)), 
    "Supplied values must not be or contain `NA` or infinite values"
  )
  expect_error(
    check_all_positive_ints(-11.00001), 
    "Supplied values must be \\(coercible to\\) integer values"
  )
  expect_error(
    check_all_positive_ints(-11, NULL), 
    "Supplied values must not be `NULL`"
  )
  expect_error(
    check_all_positive_ints(c(11.00001, 12.00000000001), 1:2), 
    "Supplied values must be \\(coercible to\\) integer values"
  )
  expect_error(
    check_all_positive_ints(11.00001, 1:2), 
    "Supplied values must be \\(coercible to\\) integer values"
  )
  
})



test_that("check_all_positive_ints() correctly throws warnings for zero counts", {
  
  # these will throw warnings and return FALSE
  expect_warning(check_all_positive_ints(2, 1:2, 0), "Zero counts detected")
  expect_warning(check_all_positive_ints(2, 1:2, 0:9), "Zero counts detected")
  
})



test_that("check_all_positive_ints() correctly returns FALSE", {
  
  # these return FALSE
  expect_false(check_all_positive_ints(2, 1:2, 0, warn_zeros = FALSE))
  expect_false(check_all_positive_ints(2, 1:2, 0:9, warn_zeros = FALSE))
  
})



test_that("check_all_positive_ints() correctly returns TRUE", {
  
  # these return FALSE# these return TRUE
  expect_true(check_all_positive_ints(1, 2, 3e12))
  expect_true(check_all_positive_ints(1, 2, 1:2))
  expect_true(check_all_positive_ints(12.0000000000000000001, 1:2))
  
  
})












