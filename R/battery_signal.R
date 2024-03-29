
globalVariables(c(
  "ci_hi", "ci_lo", "est", "est_scale"
))


#' ROR, PRR, BCPNN (MCMC and norm), obs-exp ratio of a \code{2x2} contingency table
#'
#' @author Ty Stanford <tystan@gmail.com>
#' @description ROR, PRR, BCPNN (MCMC and norm), obs-exp ratio of a \code{2x2} contingency table
#' @param a also referred to as \eqn{n_{11}}{n11} as this is the count of event of interest under exposure of interest
#' @param b also referred to as \eqn{n_{10}}{n10} as this is the count of \emph{not} event of interest under exposure of interest
#' @param c also referred to as \eqn{n_{01}}{n01} as this is the count of event of interest under \emph{not} exposure of interest
#' @param d also referred to as \eqn{n_{00}}{n00} as this is the count of \emph{not} event of interest under \emph{not} exposure of interest
#' @param orig_scale \code{FALSE} (default), do you want the estimates to be on the ratio scale (e.g., BCPNN is on the log2 scale)
#' @param alpha for construction of the \code{100*(1-alpha)\%} confidence intervals
#' @param n_mcmc number of MCMC simulations per \code{(a,b,c,d)}-tuple (for BCPNN MCMC estimation only)
#' @param alpha1 (default \code{0.5}) the numerator shrinkage parameter \code{>=0} (for observed-expected ratio estimation only)
#' @param alpha2 (default \code{0.5}) the denominator shrinkage parameter \code{>=0} (for observed-expected ratio estimation only)
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
#' @seealso 
#' \code{\link{bcpnn_mcmc_signal}}, 
#' \code{\link{bcpnn_norm_signal}}, 
#' \code{\link{ror_signal}}, 
#' \code{\link{prr_signal}},
#' \code{\link{obsexp_shrink_signal}}
#' 
#' @return 
#' The function returns a \code{data.frame} with the following columns and rows corresponding to 
#'\code{ror_signal}, \code{prr_signal}, \code{bcpnn_norm_signal}, and \code{bcpnn_mcmc_signal} function calls:
#'
#' \itemize{
#'   \item \code{analysis}
#'   \item \code{n11}
#'   \item \code{drug margin}
#'   \item \code{event margin}
#'   \item \code{n..}
#'   \item \code{E_n11}
#'   \item \code{est_name}
#'   \item \code{est_scale}
#'   \item \code{est}
#'   \item \code{var_scale}
#'   \item \code{var_est}
#'   \item \code{sd_est}
#'   \item \code{alpha}
#'   \item \code{ci_lo}
#'   \item \code{ci_hi}
#' }

battery_signal <- function(
  a, b, c, d, 
  orig_scale = FALSE, 
  alpha = 0.05, n_mcmc = 1e5,
  alpha1 = 0.5, alpha2 = 0.5
) {
  
  res <- 
    as_tibble(rbind(
      ror_signal(a, b, c, d, alpha = alpha),
      prr_signal(a, b, c, d, alpha = alpha),
      bcpnn_norm_signal(a, b, c, d, alpha = alpha),
      bcpnn_mcmc_signal(a, b, c, d, alpha = alpha, n_mcmc = n_mcmc),
      obsexp_shrink_signal(a, b, c, d, alpha1 = alpha1, alpha2 = alpha2, alpha = alpha)
    ))
  
  if (orig_scale) {
    res <-
      res %>%
      mutate(
        ci_lo = ifelse(est_scale == "log2", 2 ^ ci_lo, ci_lo),
        ci_hi = ifelse(est_scale == "log2", 2 ^ ci_hi, ci_hi),
        est   = ifelse(est_scale == "log2", 2 ^ est  , est  ),
        est_scale = ifelse(est_scale == "log2", "orig scale", est_scale)
      )
  }
  
  res
  
}


