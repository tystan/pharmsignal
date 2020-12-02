


#' Test whether multiplication of two integers will cause machine overflow
#'
#' @author Ty Stanford <tystan@gmail.com>
#' @description Test whether multiplication of two integers will cause machine overflow (if stored as an integer).
#' @param int1 the first integer value (may be a numeric/double)
#' @param int2 the second integer value (may be a numeric/double)
#' @export
#' @details
#' Two numbers are supplied and the result is boolean/logical answering the question,
#' "if we multiply the two numbers supplied, can the result be stored as an integer?"
#' That is, the supplied values need to be numeric or integer, but the function is not picky about type
#' beyond being able to multiple/divide values
#'
#' The logic used is as follows:
#'
#' if \code{a * b > max_int} where \code{a,b} integers ==> overflow, then....
#'
#' \code{a > max_int / b} ==> overflow.
#'
#' The maximum integer value (max_int above) that can be stored in your R session can be found using \code{.Machine$integer.max}.




is_mult_int_overflow <- function(int1, int2) {
  # note if a * b > .Machine$integer.max where a,b integers ==> overflow
  # then....
  # a > .Machine$integer.max / b ==> overflow

  check_all_positive_ints(int1, int2) # sanity check for positive integers

  # however, need to test int2 not == 0,
  # then we won't have int overflow (and cannot divide with)
  if (all(int2 != 0)) {
    return(any(int1 > (.Machine$integer.max / int2)))
  } else {
    return(FALSE)
  }
}


