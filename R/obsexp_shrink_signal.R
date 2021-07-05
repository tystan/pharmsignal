

#' Observed-to-expected ratios with shrinkage for a \code{2x2} contingency table
#'
#' @author Ty Stanford <tystan@gmail.com>
#' @description
#' Observed-to-expected ratios with simple shrinkage of Norén (2013) for a \code{2x2} contingency table
#' @param a also referred to as \eqn{n_{11}}{n11} as this is the count of event of interest under exposure of interest
#' @param b also referred to as \eqn{n_{10}}{n10} as this is the count of \emph{not} event of interest under exposure of interest
#' @param c also referred to as \eqn{n_{01}}{n01} as this is the count of event of interest under \emph{not} exposure of interest
#' @param d also referred to as \eqn{n_{00}}{n00} as this is the count of \emph{not} event of interest under \emph{not} exposure of interest
#' @param alpha1 (default \code{0.5}) the numerator shrinkage parameter \code{>=0}. See details.
#' @param alpha2 (default \code{0.5}) the denominator shrinkage parameter \code{>=0}. See details.
#' @param alpha for construction of the \code{100*(1-alpha)\%} confidence intervals
#' @export
#' @details
#' 
#' The observed to expected (OE) ratio with approximate \code{100*(1-alpha)\%} confidence intervals
#' are constructed on the log2 scale as outlined in Norén et al. (2013).
#' 
#' The OE ratio with shrinkage estimates is calculated as \code{(O + alpha1) / (E + alpha2)}.
#' 
#' If \code{(O + alpha1)} < \code{1}, then the exact uncertainty limits should be used. 
#' That is the \code{100*(1-alpha)\%} confidence intervals as implemented in \code{bcpnn_mcmc_signal()} (Norén et al., 2013).
#' 
#' \code{log2(OE)} approximates the Bayesian confidence propagation neural network information component (IC) 
#' with reasonable accuracy when \code{alpha1 = alpha2 = 0.5} (Norén et al., 2013).
#' 
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
#'   \item \code{est_name}: \code{"obs exp ratio"}
#'   \item \code{est_scale}: \code{"log2"}
#'   \item \code{est}: OE ratio (with shrinkage) estimate on log2 scale 
#'   \item \code{var_scale}: \code{NA}
#'   \item \code{var_est}: \code{NA}
#'   \item \code{sd_est}: \code{NA}
#'   \item \code{alpha}: alpha level supplied for the \code{100*(1-alpha)\%} confidence intervals
#'   \item \code{ci_lo}: lower bound of the \code{100*(1-alpha)\%} confidence interval
#'   \item \code{ci_hi}: upper bound of the \code{100*(1-alpha)\%} confidence interval
#' }
#'
#' @references 
#' Norén GN, Hopstadius J, Bate A. 
#' Shrinkage observed-to-expected ratios for robust and transparent large-scale pattern discovery. 
#' Statistical methods in medical research. 2013 Feb;22(1):57-69.
#' 
#' @seealso 
#' \code{\link{bcpnn_mcmc_signal}}, 
#' \code{\link{bcpnn_norm_signal}}, 
#' \code{\link{ror_signal}}, 
#' \code{\link{prr_signal}}
#' 
#' @examples
#' # Singhal et al. p409 table 16. Int J Pharm Pharm Sci, Vol 7, Issue 6, 405-411
#' obsexp_shrink_signal(28, 942, 17, 31435)
#' obsexp_shrink_signal(
#'   c("Cisplatin-Ototoxicity" = 28), 942, 17, 31435
#' )
#' obsexp_shrink_signal(122, 1320, 381, 31341)
#' obsexp_shrink_signal(
#'   c("Cisplatin-Ototoxicity" = 28, "Carboplatin-Pruritis" = 122), 
#'   c(942, 1320), 
#'   c(17, 381), 
#'   c(31435, 31341)
#' )
#' # case when (O + alpha1) < 1 [see details]
#' set.seed(1234)
#' obsexp_shrink_signal(0, 12, 49, 300) # CI the same as below for bcpnn_mcmc
#' set.seed(1234)
#' bcpnn_mcmc_signal(0, 12, 49, 300) 
#' 



obsexp_shrink_signal <- function(a, b, c, d, alpha1 = 0.5, alpha2 = 0.5, alpha = 0.05) {
  
  # make sure values are positive integers
  check_all_positive_ints(a, b, c, d, warn_zeros = FALSE) # can ignore 0 counts for BCPNN
  # make sure we have equal number of elements in each of a, b, c, d
  check_lengths_equal(a, b, c, d)
  # check alpha1 and alpha2 shrinkage params are not negative values or vectors
  check_all_non_neg_scalar(alpha1, alpha2)
  
  # inherit names of `a` vector if they exist
  analysis_names <- names(a)
  if (is.null(analysis_names)) {
    n_a <- length(a)
    n_a_digits <- nchar(as.character(n_a))
    analysis_names <- sprintf(paste0("%0", n_a_digits,".0f"), 1:n_a)
  }
  
  m_obs <- length(a)
  na_numeric <- as.numeric(NA) # NAs by default are logical data type
  na_char <- as.character(NA) # NAs by default are logical data type
  
  n.. <- a + b + c + d
  n11 <- a
  n1. <- a + b
  n.1 <- a + c
  E_n11 <- (n1. / n..) * n.1
  
  log2_oe <- log2((n11 + alpha1) / (E_n11 + alpha2))
  log2_lb <- log2_oe - 3.3 * (n11 + alpha1) ^ (-1 / 2) - 2 * (n11 + alpha1) ^ (-3 / 2)
  log2_ub <- log2_oe + 2.4 * (n11 + alpha1) ^ (-1 / 2) - 0.5 * (n11 + alpha1) ^ (-3 / 2)
    
  need_exact_lims <- (a + alpha1) < 1
  if (any(need_exact_lims)) {
    
    cat("NOTE: there is", sum(need_exact_lims), "case(s) of (O + alpha1) < 1.\n")
    cat("      These rows will use the CIs from the BCPNN IC (using MCMC estimation).\n")
    
    # use default mcmc draws
    a_star <- a[need_exact_lims]
    b_star <- b[need_exact_lims]
    c_star <- c[need_exact_lims]
    d_star <- d[need_exact_lims]
    ic <- bcpnn_mcmc_signal(a_star, b_star, c_star, d_star, alpha = alpha)
    # replace the CIs for those with (O + alpha1) < 1
    log2_lb[need_exact_lims] <- ic[["ci_lo"]]
    log2_ub[need_exact_lims] <- ic[["ci_hi"]]
    
  }
  
  res_df <-
    data.frame(
      analysis_names,
      n11, n1., n.1,  n.., E_n11,
      "obsexp_shrink", "log2", log2_oe,
      na_char, na_numeric, na_numeric,
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
