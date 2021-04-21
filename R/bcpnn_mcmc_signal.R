

#' BCPNN IC using MCMC for a \code{2x2} contingency table
#'
#' @author Ty Stanford <tystan@gmail.com>
#' @description
#' Bayesian Confidence Propagation Neural Network Information Critereon (BCPNN IC) using
#' Markov Chain Monte Carlo (MCMC) for a \code{2x2} contingency table
#' @param a also referred to as \eqn{n_{11}}{n11} as this is the count of event of interest under exposure of interest
#' @param b also referred to as \eqn{n_{10}}{n10} as this is the count of \emph{not} event of interest under exposure of interest
#' @param c also referred to as \eqn{n_{01}}{n01} as this is the count of event of interest under \emph{not} exposure of interest
#' @param d also referred to as \eqn{n_{00}}{n00} as this is the count of \emph{not} event of interest under \emph{not} exposure of interest
#' @param alpha for construction of the \code{100*(1-alpha)\%} confidence intervals
#' @param n_mcmc number of MCMC simulations per \code{(a,b,c,d)}-tuple
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
#' bcpnn_mcmc_signal(28, 942, 17, 31435)
#' bcpnn_mcmc_signal(
#'   c("Cisplatin-Ototoxicity" = 28), 942, 17, 31435
#' )
#' bcpnn_mcmc_signal(122, 1320, 381, 31341)
#' bcpnn_mcmc_signal(
#'   c("Cisplatin-Ototoxicity" = 28, "Carboplatin-Pruritis" = 122), 
#'   c(942, 1320), 
#'   c(17, 381), 
#'   c(31435, 31341)
#' )



bcpnn_mcmc_signal <- function(a, b, c, d, alpha = 0.05, n_mcmc = 1e5) {

  # make sure values are positive integers
  check_all_positive_ints(a, b, c, d, warn_zeros = FALSE) # can ignore 0 counts for BCPNN
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
  na_char <- as.character(NA) # NAs by default are logical data type

  n.. <- a + b + c + d
  n11 <- a
  n10 <- b
  n01 <- c
  n00 <- d
  n1. <- a + b
  n.1 <- a + c
  E_n11 <- (n1. / n..) * n.1

  # hyper-parameter setup etc
  q1. <- (n1. + 0.5) / (n.. + 1)
  q.1 <- (n.1 + 0.5) / (n.. + 1)
  q.0 <- (n.. - n.1 + 0.5) / (n.. + 1)
  q0. <- (n.. - n1. + 0.5) / (n.. + 1)
  a.. <- 0.5 / (q1. * q.1)
  a11 <- q1. * q.1 * a..
  a10 <- q1. * q.0 * a..
  a01 <- q0. * q.1 * a..
  a00 <- q0. * q.0 * a..
  g11 <- a11 + n11
  g10 <- a10 + n10
  g01 <- a01 + n01
  g00 <- a00 + n00
  g1. <- g11 + g10
  g.1 <- g11 + g01

  # note from Soren (2006):
  # p11 ~ Be(g11, g10 + g01 + g00) ==> E[p11] = g11 / (g11 + g10 + g01 + g00)
  # p1. ~ Be(g11 + g10, g01 + g00) ==> E[p1.] = (g11 + g10) / (g11 + g10 + g01 + g00)
  # p.1 ~ Be(g11 + g01, g10 + g00) ==> E[p.1] = (g11 + g01) / (g11 + g10 + g01 + g00)
  # therefore:
  # map_ic =
  #   log2(E[p11] / (E[p1.] * E[p.1])) =
  #   log2(g11 * (g11 + g10 + g01 + g00) / ((g11 + g10) * (g11 + g01)))
  # as E[X] = alpha / (alpha + beta) for X ~ Beta(alpha, beta)

  map_ic <-
    log2(
      g11 * (g11 + g10 + g01 + g00) / ((g11 + g10) * (g11 + g01))
    )

  m <- NULL # initialise (otherwise package check complains not globally available)
  mcmc_res <-
    foreach(m = 1:m_obs, .combine = rbind) %do% {

      p <- rdirichlet(n_mcmc, c(g11[m], g10[m], g01[m], g00[m]))
      p11 <- p[, 1]
      p1. <- p11 + p[, 2]
      p.1 <- p11 + p[, 3]
      ic_monte <- log2(p11 / (p1. * p.1))

      # posterior distribution quantiles: (0.025, 0.5, 0.975) for alpha = 0.05
      qs <- stats::quantile(ic_monte, c(alpha / 2, 0.5, 1 - alpha / 2))


      ## this is where the m.a.p. estimate is likely to differ from mean (Noren, 2006)
      # if (n11[m] <= 10) {
      #   cat("Note the empirical MCMC IC median is", qs[2], "and the m.a.p. IC est is", map_ic[m])
      #   cat(" resulting in a difference of", round(qs[2] - map_ic[m], 3), "\n")
      # }

      data.frame(
        anlys = m,
        map_ic = map_ic[m],
        empirical_median = qs[2],
        lo = qs[1],
        hi = qs[3]
      )

    }

  # not guaranteed foreach par returns in order of iteration
  mcmc_res <- mcmc_res[order(mcmc_res$anlys), ]


  res_df <-
    data.frame(
      analysis_names,
      n11, n1., n.1,  n.., E_n11,
      "bcpnn_mcmc", "log2", mcmc_res$map_ic,
      na_char, na_numeric, na_numeric,
      alpha, mcmc_res$lo, mcmc_res$hi,
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
