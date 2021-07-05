
#' Test whether input arguments are non-negative scalar
#'
#' @author Ty Stanford <tystan@gmail.com>
#' @description Test whether input arguments are non-negative scalar
#' @param ... arbitrary number of inputs (each to be checked whether they are non-negative scalar)
#' @param warn_zeros \code{TRUE} will provide a warning if any 0 values, (\code{FALSE}, default) will ignore 0 counts
#'
#' @export
#' @details
#' Five checks are made for each input:
#' \itemize{
#'   \item input must not be \code{NULL}
#'   \item input must not be or contain \code{NA} or infinite values
#'   \item input must be integer/numeric
#'   \item input must have a length of 1
#'   \item input must be 0 or greater in value
#' }
#'
#' Each argument can therefore be either a numeric/integer scalar.
#'
#' Returns TRUE if all passed arguments are non-negative scalar (pass the above tests).
#' An error is thrown if all passed arguments are NOT non-negative scalar
#' @examples
#' # these return TRUE
#' check_all_non_neg_scalar(1, 2, 3e12)
#' check_all_non_neg_scalar(0, 2, 3e12)
#'
#' # these will throw warnings and return TRUE
#' # check_all_non_neg_scalar(2, 0, warn_zeros = TRUE)
#' 
#' # these will throw errors
#' # check_all_non_neg_scalar(12.0000000000000000001, 1:2)
#' # check_all_non_neg_scalar(1, 2, 1:2)
#' # check_all_non_neg_scalar(1, character(0))
#' # check_all_non_neg_scalar(1, numeric(0))
#' # check_all_non_neg_scalar(-11, 2, 1:2)
#' # check_all_non_neg_scalar(-11, NA)
#' # check_all_non_neg_scalar(-11, c(1, NA))
#' # check_all_non_neg_scalar(-11.00001)
#' # check_all_non_neg_scalar(-11, NULL)
#' # check_all_non_neg_scalar(c(11.00001, 12.00000000001), 1:2)
#' # check_all_non_neg_scalar(11.00001, 1:2)




check_all_non_neg_scalar <- function(..., warn_zeros = FALSE) {
  
  dots <- list(...)
  ndots <- length(dots)
  approx_zero <- 1e-12
  
  for (i in 1:ndots) {
    if (TRUE %in% is.null(dots[[i]])) {
      stop("Supplied values must not be `NULL`")
    }
  }
  
  for (i in 1:ndots) {
    if (!all(is.finite(dots[[i]]) %in% TRUE)) {
      stop("Supplied values must not be or contain `NA` or infinite values")
    }
  }
  
  for (i in 1:ndots) {
    if (!all(is.numeric(dots[[i]]))) {
      stop("Supplied values must be numeric")
    }
  }
  
  for (i in 1:ndots) {
    if (!all(length(dots[[i]]) == 1)) {
      stop("Supplied values must have length 1 (scalar)")
    }
  }
  
  for (i in 1:ndots) {
    if (!all(dots[[i]] >= 0)) {
      stop("Supplied values must be non-negative")
    }
  }
  
  for (i in 1:ndots) {
    if (!all(dots[[i]] > 0)) {
      if (warn_zeros) {
        warning("Zero values detected in supplied non-negative scalars")
      }
      # return(FALSE)
    }
  }
  
  return(TRUE)
  
}




