
#' Test whether input arguments are all integers
#'
#' @author Ty Stanford <tystan@gmail.com>
#' @description Test whether input arguments are all integers
#' @param ... arbitrary number of inputs (that may be vectors themselves)
#' @param warn_zeros (\code{TRUE}, default) will provide a warning if any counts are 0, FALSE will ignore 0 counts
#'
#' @export
#' @details
#' Six checks are made for each input:
#' \itemize{
#'   \item input must not be \code{NULL}
#'   \item input must not be or contain \code{NA} or infinite values
#'   \item input must be integer/numeric
#'   \item input must have a length of 1 or more (so numeric(0) values can't sneak past)
#'   \item input must be (coercible to) integer values (i.e., not further than \code{1e-12} from a whole number)
#'   \item input must be 0 or greater in value
#' }
#'
#' Each argument can therefore be either a numeric/integer scalar or vector.
#'
#' Returns TRUE if all passed arguments are non-negative integers (pass the above tests).
#' Returns FALSE if 0 values are detected.
#' An error is thrown if all passed arguments are NOT non-negative integers
#' @examples
#' # these return TRUE
#' check_all_positive_ints(1, 2, 3e12)
#' check_all_positive_ints(1, 2, 1:2)
#' check_all_positive_ints(12.0000000000000000001, 1:2)
#'
#' # these return FALSE
#' check_all_positive_ints(2, 1:2, 0, warn_zeros = FALSE)
#' check_all_positive_ints(2, 1:2, 0:9, warn_zeros = FALSE)
#'
#' # these will throw warnings and return FALSE
#' # check_all_positive_ints(2, 1:2, 0)
#' # check_all_positive_ints(2, 1:2, 0:9)

#' # these will throw errors
#' # check_all_positive_ints(1, character(0))
#' # check_all_positive_ints(1, numeric(0))
#' # check_all_positive_ints(-11, 2, 1:2)
#' # check_all_positive_ints(-11, NA)
#' # check_all_positive_ints(-11, c(1, NA))
#' # check_all_positive_ints(-11.00001)
#' # check_all_positive_ints(-11, NULL)
#' # check_all_positive_ints(c(11.00001, 12.00000000001), 1:2)
#' # check_all_positive_ints(11.00001, 1:2)




check_all_positive_ints <- function(..., warn_zeros = TRUE) {

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
    if (!all(length(dots[[i]]) > 0)) {
      stop("Supplied values must have length 1 (scalar) or greater (vector)")
    }
  }

  for (i in 1:ndots) {
    if (!all(abs(dots[[i]] - floor(dots[[i]])) <= approx_zero)) {
      stop("Supplied values must be (coercible to) integer values")
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
        warning("Zero counts detected")
      }
      return(FALSE)
    }
  }

  return(TRUE)

}




