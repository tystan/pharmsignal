

#' BCPNN IC using normal approximation CIs for a \code{2x2} contingency table
#'
#' @author Ty Stanford <tystan@gmail.com>
#' @description
#' Bayesian Confidence Propagation Neural Network Information Component (BCPNN IC) using
#' normal approximation CIs for a \code{2x2} contingency table
#' @param a also referred to as \eqn{n_{11}} as this is the count of event of interest under exposure of interest
#' @param b also referred to as \eqn{n_{10}} as this is the count of \emph{not} event of interest under exposure of interest
#' @param c also referred to as \eqn{n_{01}} as this is the count of event of interest under \emph{not} exposure of interest
#' @param d also referred to as \eqn{n_{00}} as this is the count of \emph{not} event of interest under \emph{not} exposure of interest
#' @param alpha for construction of the \code{100*(1-alpha)\%} confidence intervals
#' @export
#' @details
#' It is assumed that the contingency table under consideration has drugs/exposures in the rows and outcomes/events in the columns.
#'
#' Let's assume we are interested in drug Y and outcome Z, the contingency table will look like this
#'
#' \tabular{ccc}{
#' . \tab outcome Z \tab not outcome Z\cr
#' drug Y \tab \code{a}=\eqn{n_{11}} \tab \code{b}=\eqn{n_{10}} \cr
#' not drug Y \tab \code{c}=\eqn{n_{01}}  \tab \code{d}=\eqn{n_{00}} \cr
#' }
#'
#' The function returns a \code{data.frame} with the following columns:
#'
#' \itemize{
#'   \item \code{analysis}: names from the supplied vector \code{a} if present, the values 1:length(a) (with 0 padding) if not
#'   \item \code{n11}: drug Y and outcome Z count (\code{a})
#'   \item \code{drug margin}: margin for drug Y (\code{a + b})
#'   \item \code{event margin}: margin for outcome Z (\code{a + c})
#'   \item \code{n..}: total count from contingency table (\code{a + b + c + d})
#'   \item \code{E_n11}: expected value of \code{n11} based on margins (\code{(a + b) * (a + c) / (a + b + c + d)})
#'   \item \code{est_name}: \code{"ror"}
#'   \item \code{est_scale}: \code{"orig scale"}
#'   \item \code{est}: ROR estimate on original scale (not log_e)
#'   \item \code{var_scale}: \code{"ln"}
#'   \item \code{var_est}: variance of the log_e(ROR) estimate
#'   \item \code{sd_est}: standard deviation of the log_e(ROR) estimate
#'   \item \code{alpha}: alpha level supplied for the \code{100*(1-alpha)\%} confidence intervals
#'   \item \code{ci_lo}: lower bound of the \code{100*(1-alpha)\%} confidence interval
#'   \item \code{ci_hi}: upper bound of the \code{100*(1-alpha)\%} confidence interval
#' }
#'
#' @examples
#' # Singhal et al. p409 table 16. Int J Pharm Pharm Sci, Vol 7, Issue 6, 405-411
#' bcpnn_norm_signal(28, 942, 17, 31435)
#' bcpnn_norm_signal(
#'   c("Cisplatin-Ototoxicity" = 28), 942, 17, 31435
#' )
#' bcpnn_norm_signal(122, 1320, 381, 31341)
#' bcpnn_norm_signal(
#'   c("Cisplatin-Ototoxicity" = 28, "Carboplatin-Pruritis" = 122), 
#'   c(942, 1320), 
#'   c(17, 381), 
#'   c(31435, 31341)
#' )



bcpnn_norm_signal <- function(a, b, c, d, alpha = 0.05) {
  
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
  
  m_obs <- length(a)
  na_numeric <- as.numeric(NA) # NAs by default are logical data type
  
  n.. <- a + b + c + d
  n11 <- a
  n1. <- a + b
  n.1 <- a + c
  E_n11 <- (n1. / n..) * n.1
  
  p1 <- 1 + n1.
  p2 <- 1 + n.. - n1.
  q1 <- 1 + n.1
  q2 <- 1 + n.. - n.1
  r1 <- 1 + n11
  r2b <- n.. - n11 - 1 + (2 + n..)^2 / (q1 * p1)
  
  # digamma and trigamma return the first and second derivatives of the logarithm of the gamma function
  # digamma(x) = psigamma(x, deriv = 1) = d/dx{ln Gam(x)} = Gam'(x) / Gam(x)
  E_ic <- 
    log(2)^(-1) * (
      digamma(r1) - digamma(r1 + r2b) - (
        digamma(p1) - digamma(p1 + p2) + digamma(q1) - digamma(q1 + q2)
      )
    )
  var_ic <- 
    log(2)^(-2) * (
      trigamma(r1) - trigamma(r1 + r2b) + (
        trigamma(p1) - trigamma(p1 + p2) + trigamma(q1) - trigamma(q1 + q2)
      )
    )
  
  log2_lb <- qnorm(alpha / 2, E_ic, sqrt(var_ic))
  log2_ub <- qnorm(1 - alpha / 2, E_ic, sqrt(var_ic))
  
  res_df <-
    data.frame(
      analysis_names,
      n11, n1., n.1,  n.., E_n11,
      "bcpnn_norm", "log2", E_ic,
      "log2", var_ic, sqrt(var_ic),
      alpha, log2_lb, log2_ub,
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
