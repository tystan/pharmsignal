
globalVariables(c(
  "ci_hi", "ci_lo", "est", "est_scale"
))


#' ROR, PRR, BCPNN (MCMC and norm) of a \code{2x2} contingency table
#'
#' @author Ty Stanford <tystan@gmail.com>
#' @description ROR, PRR, BCPNN (MCMC and norm) of a \code{2x2} contingency table
#' @param a also referred to as \eqn{n_{11}} as this is the count of event of interest under exposure of interest
#' @param b also referred to as \eqn{n_{10}} as this is the count of \emph{not} event of interest under exposure of interest
#' @param c also referred to as \eqn{n_{01}} as this is the count of event of interest under \emph{not} exposure of interest
#' @param d also referred to as \eqn{n_{00}} as this is the count of \emph{not} event of interest under \emph{not} exposure of interest
#' @param label \code{"exposure v outcome"} (default), string to append as column to data.frame describing data 
#' @param orig_scale \code{FALSE} (default), do you want the estimates to be on the ratio scale (e.g., BCPNN is on the log2 scale)
#' @param alpha for construction of the \code{100*(1-alpha)\%} confidence intervals
#' @param n_mcmc number of MCMC simulations per \code{(a,b,c,d)}-tuple
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
  label = "exposure v outcome", orig_scale = FALSE, 
  alpha = 0.05, n_mcmc = 1e5
) {
  
  res <- rbind(
    ror_signal(a, b, c, d, alpha = alpha),
    prr_signal(a, b, c, d, alpha = alpha),
    bcpnn_norm_signal(a, b, c, d, alpha = alpha),
    bcpnn_mcmc_signal(a, b, c, d, alpha = alpha, n_mcmc = n_mcmc)
  )
  res <- as_tibble(res)
  res$lab <- label
  
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


