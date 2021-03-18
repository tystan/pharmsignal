


# ---- bcpnn_norm ----



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



# ---- prr ----



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







