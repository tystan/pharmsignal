
#' Test whether input arguments vectors of the same length
#'
#' @author Ty Stanford <tystan@gmail.com>
#' @description Test whether input arguments vectors of the same length
#' @param ... arbitrary number of vectors
#'
#' @export
#' @details
#' It is assumed the input has been checked using \code{check_all_positive_ints(...)} before using this function.
#' Each argument can therefore be either a numeric/integer scalar or vector.
#'
#' Returns TRUE if all vectors supplied are of the same length. Error thrown if not.
#' @examples
#' ### these return TRUE
#' check_lengths_equal(1)
#' check_lengths_equal(1:3)
#' check_lengths_equal(1, 2, 3e12)
#' check_lengths_equal(1:3, 2:4, 3:5)
#' check_lengths_equal(12.0000000000000000000001, 2)
#'
#' ### these will throw errors
#' # check_lengths_equal(2, 1:2, 0)
#' # check_lengths_equal(1:9, 0:9)






check_lengths_equal <- function(...) {

  dots <- list(...)
  ndots <- length(dots)
  first_len <- length(dots[[1]])

  for (i in 1:ndots) {
    if (length(dots[[i]]) != first_len) {
      stop("Supplied vectors are not the same length")
    }
  }

  return(TRUE)

}




