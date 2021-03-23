











#' Proportional reporting ratio (PRR) of a \code{2x2} contingency table
#'
#' @author Ty Stanford <tystan@gmail.com>
#' @description Proportional reporting ratio (PRR) of a \code{2x2} contingency table
#' @param a also referred to as \eqn{n_{11}}{n11} as this is the count of event of interest under exposure of interest
#' @param b also referred to as \eqn{n_{10}}{n10} as this is the count of \emph{not} event of interest under exposure of interest
#' @param c also referred to as \eqn{n_{01}}{n01} as this is the count of event of interest under \emph{not} exposure of interest
#' @param d also referred to as \eqn{n_{00}}{n00} as this is the count of \emph{not} event of interest under \emph{not} exposure of interest
#' @param alpha for construction of the \code{100*(1-alpha)\%} confidence intervals
#' @export
#' @details
#' Note that the \code{a},\code{b},\code{c},\code{d} inputs may be vectors of equal length, for which the function
#' will perform the calculations for each individual \code{(a,b,c,d)}-tuple moving across the vectors. 
#' 
#' It is assumed that the contingency table under consideration has drugs/exposures in the rows and outcomes/events in the columns.
#'
#' Let's assume we are interested in drug Y and outcome Z, the contingency table will look like this
#'
#' \tabular{ccc}{
#' . \tab outcome Z \tab not outcome Z\cr
#' drug Y \tab \code{a}=\eqn{n_{11}}{n11} \tab \code{b}=\eqn{n_{10}}{n10} \cr
#' not drug Y \tab \code{c}=\eqn{n_{01}}{n01}  \tab \code{d}=\eqn{n_{00}}{n00} \cr
#' }
#'
#' @return 
#' The function returns a \code{data.frame} with the following columns:
#'
#' \itemize{
#'   \item \code{analysis}: names from the supplied vector \code{a} if present, the values 1:length(a) (with 0 padding) if not
#'   \item \code{n11}: drug Y and outcome Z count (\code{a})
#'   \item \code{drug margin}: margin for drug Y (\code{a + b})
#'   \item \code{event margin}: margin for outcome Z (\code{a + c})
#'   \item \code{n..}: total count from contingency table (\code{a + b + c + d})
#'   \item \code{E_n11}: expected value of \code{n11} based on margins (\code{(a + b) * (a + c) / (a + b + c + d)})
#'   \item \code{est_name}: \code{"prr"}
#'   \item \code{est_scale}: \code{"orig scale"}
#'   \item \code{est}: ROR estimate on original scale (not log_e)
#'   \item \code{var_scale}: \code{"ln"}
#'   \item \code{var_est}: variance of the log_e(PRR) estimate
#'   \item \code{sd_est}: standard deviation of the log_e(PRR) estimate
#'   \item \code{alpha}: alpha level supplied for the \code{100*(1-alpha)\%} confidence intervals
#'   \item \code{ci_lo}: lower bound of the \code{100*(1-alpha)\%} confidence interval
#'   \item \code{ci_hi}: upper bound of the \code{100*(1-alpha)\%} confidence interval
#' }
#'
#' @examples
#' # Singhal et al. p409 table 16. Int J Pharm Pharm Sci, Vol 7, Issue 6, 405-411
#' prr_signal(28, 942, 17, 31435)
#' prr_signal(c("Cisplatin-Ototoxicity" = 28), 942, 17, 31435)
#' prr_signal(122, 1320, 381, 31341)
#' prr_signal(
#'   c("Cisplatin-Ototoxicity" = 28, "Carboplatin-Pruritis" = 122), 
#'   c(942, 1320), 
#'   c(17, 381), 
#'   c(31435, 31341)
#' )



prr_signal <- function(a, b, c, d, alpha = 0.05) {
  
  # make sure values are positive integers
  check_all_positive_ints(a, b, c, d)
  # make sure we have equal number of elements in each of a, b, c, d
  check_lengths_equal(a, b, c, d)
  
  # inherit names of `a` vector if they exist
  analysis_names <- names(a)
  if (is.null(analysis_names)) {
    n_a <- length(a)
    n_a_digits <- nchar(as.character(n_a))
    analysis_names <- sprintf(paste0("%0", n_a_digits,".0f"), 1:n_a)
  }
  
  n.. <- a + b + c + d
  n1. <- a + b
  n.1 <- a + c
  E_n11 <- (n1. / n..) * n.1
  
  n11 <- a
  n10 <- b
  n01 <- c
  n00 <- d
  
  log_prr <- log((n11 / (n11 + n10)) / (n01 / (n01 + n00)))
  var_log_prr <- 1 / n11 - 1 / (n11 + n10) + 1 / n01 - 1 / (n01 + n00)
  
  
  log_lb <- qnorm(alpha / 2, log_prr, sqrt(var_log_prr))
  log_ub <- qnorm(1 - alpha / 2, log_prr, sqrt(var_log_prr))
  
  res_df <-
    data.frame(
      analysis_names,
      a, n1., n.1, n.., E_n11,
      "prr", "orig scale", exp(log_prr),
      "ln", var_log_prr, sqrt(var_log_prr),
      alpha, exp(log_lb), exp(log_ub),
      stringsAsFactors = FALSE # for R versions < 4.0
    )
  
  colnames(res_df) <-
    c(
      "analysis",
      "n11", "drug margin", "event margin", "n..", "E_n11",
      "est_name", "est_scale", "est",
      "var_scale", "var_est", "sd_est",
      "alpha", "ci_lo", "ci_hi"
    )
  
  rownames(res_df) <- NULL
  
  return(res_df)
  
}



