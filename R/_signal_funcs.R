
library("testthat")
library("dplyr")
library("MCMCpack")
library("foreach")
# won't load PhViD but poackage is required for test below
# library("PhViD")
# won't load DescTools but poackage is required for test below
# library("DescTools")

# ---- bcpnn_norm ----

bcpnn_norm <- function(a, b, c, d, RR0 = 1, alpha = 0.05) {

  N <- a + b + c + d
  n11 <- a
  n1. <- a + b
  n.1 <- a + c
  n10 <- n1. - n11
  n01 <- n.1 - n11
  n00 <- N - (n11 + n10 + n01)
  E <- (n1. / N) * n.1

  Nb.Cell <- length(n11)
  post.H0 <- matrix(nrow = Nb.Cell, ncol = length(RR0))
  p1 <- 1 + n1.
  p2 <- 1 + N - n1.
  q1 <- 1 + n.1
  q2 <- 1 + N - n.1
  r1 <- 1 + n11
  r2b <- N - n11 - 1 + (2 + N)^2/(q1 * p1)
  # digamma and trigamma return the first and second derivatives of the logarithm of the gamma function
  # digamma(x) = psigamma(x, deriv = 1) = d/dx{ln Gam(x)} = Gam'(x) / Gam(x)
  EICb <- 
    log(2)^(-1) * (
      digamma(r1) - digamma(r1 + r2b) - (
        digamma(p1) - digamma(p1 + p2) + digamma(q1) - digamma(q1 + q2)
      )
    )
  VICb <- 
    log(2)^(-2) * (
      trigamma(r1) - trigamma(r1 + r2b) + (
        trigamma(p1) - trigamma(p1 + p2) + trigamma(q1) - trigamma(q1 + q2)
      )
    )
  post.H0 <- pnorm(log(RR0), EICb, sqrt(VICb))
  LB <- qnorm(alpha / 2, EICb, sqrt(VICb))
  UB <- qnorm(1 - alpha / 2, EICb, sqrt(VICb))

  RES <- 
    data.frame(
      n11, n1., n.1, N, E, 
      "bcpnn_norm", "log2",
      EICb, VICb, sqrt(VICb),
      alpha, LB, UB,
      stringsAsFactors = FALSE
    )
  colnames(RES) <- 
    c("n11", 
      "drug margin", "event margin", "n..",
      "E_n11", 
      "est type", "est scale",
      "est", "V_est", "sd_est", "alpha",
      "ci_lo", "ci_hi"
    )

  attr(RES, "INPUT.PARAM") <- data.frame(RR0, alpha)
    
  RES
}


# tests

# Singhal et al. p409 table 15. Int J Pharm Pharm Sci, Vol 7, Issue 6, 405-411
bcpnn_norm(28, 942, 17, 31435)
bcpnn_norm(122, 1320, 381, 31341)
bcpnn_norm(c(28, 122), c(942, 1320), c(17, 381), c(31435, 31341))

# see Statist. Med. 2006; 25:3740–3757
bcpnn_norm(25, 1126 - 25, 87 - 25, 572573 - 25 - (1126 - 25) - (87 - 25))

# AL Gould et al. Pharmacoepidemiology and Drug Safety, 2003; 12: 559–574
bcpnn_norm(1614, 85304 - 1614, 71209 - 1614, 4640498 - 1614 - (85304 - 1614) - (71209 - 1614))
# bcpnn_norm(1614, 85304 - 1614, 71209 - 1614, 4864480 - 1614 - (85304 - 1614) - (71209 - 1614))


attr(bcpnn_norm(c(28, 122), c(942, 1320), c(17, 381), c(31435, 31341)), "INPUT.PARAM")
attr(rbind(
  bcpnn_norm(28, 942, 17, 31435),
  bcpnn_norm(122, 1320, 381, 31341)
), "INPUT.PARAM")

expect_equal(
  rbind(
    bcpnn_norm(28, 942, 17, 31435),
    bcpnn_norm(122, 1320, 381, 31341)
  ),
  bcpnn_norm(c(28, 122), c(942, 1320), c(17, 381), c(31435, 31341))
)

outcome_interest <- "Ototoxicity"
drug_interest <- "Cisplatin"

(phv_test <- data.frame(
  exposure = c("Cisplatin", "Cisplatin", "not Cisplatin", "not Cisplatin"),
  outcome = c("Ototoxicity", "not Ototoxicity", "Ototoxicity", "not Ototoxicity"),
  n = c(28, 942, 17, 31435)
))

phvid_dat <- PhViD::as.PhViD(phv_test) # to look at structure: str(phvid_dat)
bcpnn_sym <- 
  PhViD::BCPNN(phvid_dat, MIN.n11 = 3, DECISION = 3, RANKSTAT = 2, MC = FALSE)$SIGNALS  %>% 
  dplyr::filter(`event effect` == outcome_interest, `drug code` == drug_interest)
  
expect_equal(
  bcpnn_norm(28, 942, 17, 31435)$`ci_lo`,
  bcpnn_sym$`Q_0.025(log(IC))`
)
  
  

# ---- bcpnn_mcmc ----

bcpnn_mcmc <- function(a, b, c, d, RR0 = 1, alpha = 0.05, NB.MC = 10000) {

  N <- a + b + c + d
  n11 <- a
  n1. <- a + b
  n.1 <- a + c
  n10 <- n1. - n11
  n01 <- n.1 - n11
  n00 <- N - (n11 + n10 + n01)
  E <- (n1. / N) * n.1

  # n1. <- n11 + n10
  # n.1 <- n11 + n01
  Nb_Obs <- length(n11)
  q1. <- (n1. + 0.5) / (N + 1)
  q.1 <- (n.1 + 0.5) / (N + 1)
  q.0 <- (N - n.1 + 0.5) / (N + 1)
  q0. <- (N - n1. + 0.5) / (N + 1)
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
  post.H0 <- vector(length = Nb_Obs)
  LB <- UB <- EIC <- vector(length = Nb_Obs)
  for (m in 1:Nb_Obs) {
    p <- rdirichlet(NB.MC, c(g11[m], g10[m], g01[m], g00[m]))
    p11 <- p[, 1]
    p1. <- p11 + p[, 2]
    p.1 <- p11 + p[, 3]
    IC_monte <- log2(p11 / (p1. * p.1))
    temp <- IC_monte < log2(RR0)
    post.H0[m] <- sum(temp) / NB.MC
    # LB[m] <- sort(IC_monte)[round(NB.MC * (alpha / 2))]
    # UB[m] <- sort(IC_monte)[round(NB.MC * (1 - alpha / 2))]
    # EIC[m] <- sort(IC_monte)[round(NB.MC / 2)]
    qs <- stats::quantile(IC_monte, c(alpha / 2, 0.5, 1 - alpha / 2))
    # print(qs)
    LB[m] <- qs[1]
    EIC[m] <- qs[2]
    UB[m] <- qs[3]
    
    # note from Soren (2006):
    # p11 ~ Be(g11, g10 + g01 + g00) ==> E[p11] = g11 / (g11 + g10 + g01 + g00)
    # p1. ~ Be(g11 + g10, g01 + g00) ==> E[p1.] = (g11 + g10) / (g11 + g10 + g01 + g00)
    # p.1 ~ Be(g11 + g01, g10 + g00) ==> E[p.1] = (g11 + g01) / (g11 + g10 + g01 + g00)
    # therefore:
    # map_IC = log2(E[p11] / (E[p1.] * E[p.1])) =
    # log2(g11 * (g11 + g10 + g01 + g00) / ((g11 + g10) * (g11 + g01)))
    
    map_IC <- log2(g11 * (g11 + g10 + g01 + g00) / ((g11 + g10) * (g11 + g01)))
    
    cat("Note the empirical MCMC IC mean is", qs[2], "and the m.a.p. IC est is", map_IC)
    cat(" resulting in a difference of", round(qs[2] - map_IC, 3), "\n")
    
    
    EIC[m] <- map_IC
    
  }
  rm(p11, p1., p.1, IC_monte, temp)
  gc()
  
  RES <- 
    data.frame(
      n11, n1., n.1,  N, E, 
      "bcpnn_mcmc", "log2",
      EIC, NA, NA,
      alpha, LB, UB,
      stringsAsFactors = FALSE
    )
  colnames(RES) <- 
    c(
      "n11", "drug margin", "event margin", "n..", "E_n11", 
      "est type", "est scale",
      "est", "V_est", "sd_est",
      "alpha", "ci_lo", "ci_hi"
    )
  
  attr(RES, "INPUT.PARAM") <- data.frame(RR0, alpha)
  
  RES
}


### check new function has results equiv to old function
set.seed(1234)
bcpnn_mc0 <- bcpnn_mcmc(28, 942, 17, 31435, NB.MC = 1e6)
bcpnn_mc0


set.seed(1234)
bcpnn_mc <- 
  PhViD::BCPNN(phvid_dat, MIN.n11 = 3, DECISION = 3, RANKSTAT = 2, MC = TRUE, NB.MC = 1e6)$SIGNALS  %>% 
  dplyr::filter(`event effect` == outcome_interest, `drug code` == drug_interest)

expect_equal(
  log(2 ^ bcpnn_mc0$`ci_lo`), #note in the source code phvid uses natural logs instead of log2 ?
  bcpnn_mc$`Q_0.025(log(IC))`,
  tol = 1e-3
)

### confirm new function has sensible results

# Table 1 of soren 2006: headache, should approx be equiv to bcpnn_norm
bcpnn_mcmc(1614, 85304 - 1614, 71209 - 1614, 4640498 - 1614 - (85304 - 1614) - (71209 - 1614), NB.MC = 1e6)
bcpnn_norm(1614, 85304 - 1614, 71209 - 1614, 4640498 - 1614 - (85304 - 1614) - (71209 - 1614))

# Table 1 of soren 2006: akathisia, should approx be equiv to bcpnn_norm
bcpnn_mcmc(328, 85304 - 328, 3001 - 328, 5019555 - 328 - (85304 - 328) - (3001 - 328), NB.MC = 1e6)
bcpnn_norm(328, 85304 - 328, 3001 - 328, 5019555 - 328 - (85304 - 328) - (3001 - 328))


# Singhal et al. p409 table 15. Int J Pharm Pharm Sci, Vol 7, Issue 6, 405-411
bcpnn_mcmc(28, 942, 17, 31435)
bcpnn_norm(28, 942, 17, 31435)
bcpnn_mcmc(122, 1320, 381, 31341)
bcpnn_norm(122, 1320, 381, 31341)

# see Statist. Med. 2006; 25:3740–3757
bcpnn_mcmc(25, 1126 - 25, 87 - 25, 572573 - 25 - (1126 - 25) - (87 - 25))
bcpnn_norm(25, 1126 - 25, 87 - 25, 572573 - 25 - (1126 - 25) - (87 - 25))

# ---- ror ----

is_mult_int_overflow <- function(int1, int2) {
  # note if a * b > .Machine$integer.max where a,b integers ==> overflow
  # then....
  # a > .Machine$integer.max / b ==> overflow
  
  # however, need to test int2 not == 0, 
  # then we won't have int overflow (and cannot divide with)
  if (all(int2 != 0)) {
    return(any(int1 > (.Machine$integer.max / int2)))
  } else {
    return(FALSE)
  }
}

ror <- function(a, b, c, d, RR0 = 1, alpha = 0.05) {

  N <- a + b + c + d
  n11 <- a
  n1. <- a + b
  n.1 <- a + c
  n10 <- n1. - n11
  n01 <- n.1 - n11
  n00 <- N - (n11 + n10 + n01)
  E <- (n1. / N) * n.1
  
  # have to check for integer overflow in R for large values (R has 32 bit (small) integers)
  int_overflow <- is_mult_int_overflow(n10, n01)
  
  logROR <- 0
  if (int_overflow) {
    # turn integers into doubles is the best way to handle
    # may lose precision but likely only REALLY LARGE ints
    logROR <- log((n11 / (as.numeric(n10) * as.numeric(n01))) * n00)
  } else {
    logROR <- log((n11 / (n10 * n01)) * n00)
  }
  var.logROR <- 1/n11 + 1/n10 + 1/n01 + 1/n00
  
  log_LB <- qnorm(alpha / 2, logROR, sqrt(var.logROR))
  log_UB <- qnorm(1 - alpha / 2, logROR, sqrt(var.logROR))
 
  RES <- 
    data.frame(
      n11, n1., n.1,  N, E, 
      "ror", "orig scale",
      exp(logROR), var.logROR, sqrt(var.logROR),
      alpha, exp(log_LB), exp(log_UB),
      stringsAsFactors = FALSE
    )
  colnames(RES) <- 
    c(
      "n11", "drug margin", "event margin", "n..", "E_n11", 
      "est type", "est scale",
      "est", "V_est", "sd_est",
      "alpha", "ci_lo", "ci_hi"
    )
  
  attr(RES, "INPUT.PARAM") <- data.frame(RR0, alpha)
  
  RES
}


# Singhal et al. p409 table 16. Int J Pharm Pharm Sci, Vol 7, Issue 6, 405-411
ror(28, 942, 17, 31435)
ror(122, 1320, 381, 31341)
ror(c(28, 122), c(942, 1320), c(17, 381), c(31435, 31341))

expect_equal(
  rbind(
    ror(28, 942, 17, 31435),
    ror(122, 1320, 381, 31341)
  ),
  ror(c(28, 122), c(942, 1320), c(17, 381), c(31435, 31341))
)

ror_est <- (28/942)/(17/31435)
ror_se <- sqrt(1 / 28 + 1 / 942 + 1 / 17 + 1 / 31435)
expect_equal(
  ror(28, 942, 17, 31435)$est,
  ror_est
)
expect_equal(
  ror(28, 942, 17, 31435)$ci_lo,
  exp(log(ror_est) - qnorm(0.975) * ror_se)
)
expect_equal(
  ror(28, 942, 17, 31435)$ci_hi,
  exp(log(ror_est) + qnorm(0.975) * ror_se)
)


external_ror <-
  DescTools::OddsRatio(
    matrix(c(28, 942, 17, 31435), byrow = TRUE, ncol = 2), 
    conf.level = 0.95, 
    method = "wald"
  )

expect_equal(
  ror(28, 942, 17, 31435)$est,
  external_ror[["odds ratio"]]
)
expect_equal(
  ror(28, 942, 17, 31435)$ci_lo,
  external_ror[["lwr.ci"]]
)
expect_equal(
  ror(28, 942, 17, 31435)$ci_hi,
  external_ror[["upr.ci"]]
)

# ---- prr ----


prr <- function(a, b, c, d, RR0 = 1, alpha = 0.05) {

  N <- a + b + c + d
  n11 <- a
  n1. <- a + b
  n.1 <- a + c
  n10 <- n1. - n11
  n01 <- n.1 - n11
  n00 <- N - (n11 + n10 + n01)
  E <- (n1. / N) * n.1

  logPRR <- log((n11/(n11 + n10))/(n01/(n01 + n00)))
  var.logPRR <- 1/n11 - 1/(n11 + n10) + 1/n01 - 1/(n01 + n00)
  
  log_LB <- qnorm(alpha / 2, logPRR, sqrt(var.logPRR))
  log_UB <- qnorm(1 - alpha / 2, logPRR, sqrt(var.logPRR))

  RES <- 
    data.frame(
      n11, n1., n.1,  N, E, 
      "prr", "orig scale",
      exp(logPRR), var.logPRR, sqrt(var.logPRR),
      alpha, exp(log_LB), exp(log_UB),
      stringsAsFactors = FALSE
    )
  colnames(RES) <- 
    c(
      "n11", "drug margin", "event margin", "n..", "E_n11", 
      "est type", "est scale",
      "est", "V_est", "sd_est",
      "alpha", "ci_lo", "ci_hi"
    )
  
  attr(RES, "INPUT.PARAM") <- data.frame(RR0, alpha)
  
  RES
}

# Singhal et al. p409 table 16. Int J Pharm Pharm Sci, Vol 7, Issue 6, 405-411
prr(28, 942, 17, 31435)
prr(122, 1320, 381, 31341)
prr(c(28, 122), c(942, 1320), c(17, 381), c(31435, 31341))

expect_equal(
  rbind(
    prr(28, 942, 17, 31435),
    prr(122, 1320, 381, 31341)
  ),
  prr(c(28, 122), c(942, 1320), c(17, 381), c(31435, 31341))
)

prr_est <- (28/(28 + 942))/(17/(17 + 31435))
prr_se <- sqrt(1 / 28 - 1 / (28 + 942) + 1 / 17 - 1 / (17 + 31435))
expect_equal(
  prr(28, 942, 17, 31435)$est,
  prr_est
)
expect_equal(
  prr(28, 942, 17, 31435)$ci_lo,
  exp(log(prr_est) - qnorm(0.975) * prr_se)
)
expect_equal(
  prr(28, 942, 17, 31435)$ci_hi,
   exp(log(prr_est) + qnorm(0.975) * prr_se)
)

# ---- battery_signal ----


battery_signal <- function(a, b, c, d, label = "exposure v outcome", orig_scale = FALSE, RR0 = 1, alpha = 0.05, MC.seed = 1234, NB.MC = 1e6) {
  
  set.seed(MC.seed)
  res <- rbind(
    ror(a, b, c, d, RR0 = RR0, alpha = alpha),
    prr(a, b, c, d, RR0 = RR0, alpha = alpha),
    bcpnn_norm(a, b, c, d, RR0 = RR0, alpha = alpha),
    bcpnn_mcmc(a, b, c, d, RR0 = RR0, alpha = alpha, NB.MC = NB.MC)
  )
  res <- as_tibble(res)
  res$lab <- label
  if (orig_scale) {
    res <-
      res %>%
      mutate(
        ci_lo = ifelse(`est scale` == "log2", 2 ^ ci_lo, ci_lo),
        ci_hi = ifelse(`est scale` == "log2", 2 ^ ci_hi, ci_hi),
        est   = ifelse(`est scale` == "log2", 2 ^ est  , est  ),
        `est scale` = ifelse(`est scale` == "log2", "orig scale", `est scale`)
      )
  }
  res
  
}

set.seed(1234)
exp_df <-
  rbind(
    ror(28, 942, 17, 31435),
    prr(28, 942, 17, 31435),
    bcpnn_norm(28, 942, 17, 31435),
    bcpnn_mcmc(28, 942, 17, 31435)
  )
exp_df <- as_tibble(exp_df)
exp_df$lab <- "Cisplatin v Ototoxicity"
expect_equal(
  battery_signal(28, 942, 17, 31435, label = "Cisplatin v Ototoxicity")$est,
  exp_df$est,
  tol = 1e-2
)
expect_equal(
  battery_signal(28, 942, 17, 31435, label = "Cisplatin v Ototoxicity")$ci_lo,
  exp_df$ci_lo,
  tol = 1e-2
)
expect_equal(
  battery_signal(28, 942, 17, 31435, label = "Cisplatin v Ototoxicity")$ci_hi,
  exp_df$ci_hi,
  tol = 1e-2
)

battery_signal(28, 942, 17, 31435, label = "Cisplatin v Ototoxicity", orig_scale = FALSE)
battery_signal(28, 942, 17, 31435, label = "Cisplatin v Ototoxicity", orig_scale = TRUE)




# ---- nx2_to_signal ----

# n rows assumed to be n drugs, 
# 2 cols meant to be 2 outcomes
nx2_to_signal <- function(nx2_tab, orig_scale = FALSE, RR0 = 1, alpha = 0.05, MC.seed = 1234, NB.MC = 1e6) {
  
  tab <- as.data.frame(nx2_tab, stringsAsFactors = FALSE)
  colnames(tab) <- c("drg", "out", "n")
  nr <- nrow(nx2_tab)
  rnms <- rownames(nx2_tab)
  cnms <- colnames(nx2_tab)
  onevall <-
    foreach(i = 1:nr, .combine = bind_rows) %do% {
      this_drg <- rnms[i] 
      not_this_drg <- paste0("* (excl ", this_drg, ")")
      tmp_df <- 
        tab %>%
        dplyr::mutate(drg = if_else(drg == this_drg, drg, not_this_drg)) %>%
        group_by(drg, out) %>% 
        summarise(n_ = sum(n))
      eoi <- tmp_df %>% dplyr::filter(drg == this_drg, out == cnms[2]) %>% pull(n_)
      if (eoi > 2) {
        battery_signal(
          eoi,
          tmp_df %>% dplyr::filter(drg == this_drg, out == cnms[1]) %>% pull(n_),
          tmp_df %>% dplyr::filter(drg == not_this_drg, out == cnms[2]) %>% pull(n_),
          tmp_df %>% dplyr::filter(drg == not_this_drg, out == cnms[1]) %>% pull(n_),
          label = paste(this_drg, "/", not_this_drg, " vs ", cnms[2], "/", cnms[1], sep = ""), 
          orig_scale = orig_scale, RR0 = RR0, alpha = alpha, MC.seed = MC.seed, NB.MC = NB.MC
        )
      } else {
        NULL
      }
    }
  onevall
}





