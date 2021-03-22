


# `pharmsignal`

Safety signal detection in pharmacological data (pharmacovigilance).

Implementations of Reporting Odds Ratio (ROR),
Proportional Reporting Ratio (PRR),
Bayesian Confidence Propagation Neural Network Information Component
(BCPNN IC, both MCMC and normal approximation) and
Sequence Symmetry Analysis (SSA).


## Background

Signal detection of adverse drug events often requires analysis of tabulated count data as seen in the below table.

|         | `Event(s) X`| `Event(s) Y` |
|:--------|---------:|----------:|
|**`Drug(s) A`**|       *a*|        *b*|
|**`Drug(s) B`**|       *c*|        *d*|

This package has functions to calculate the following signal detection methods:

* `ror_signal()`: ROR with `100(1-alpha/2)%` confidence intervals 
* `ror_signal()`: PRR with `100(1-alpha/2)%` confidence intervals 
* `bcpnn_norm_signal()`: BCPNN IC using the 'exact' expectation and variance of the the IC posterior distribution of Gould (2003) using a normal approximation to construct the `100(1-alpha/2)%` confidence interval on the log_2 scale (see equations (A3) and (A4) of the Appendix)
* `bcpnn_mcmc_signal()`: BCPNN IC using the maximum a posteriori (m.a.p.) central estimate of the IC with MCMC simulation of the exact empirical distribution for `100(1-alpha/2)%` confidence (credible) regions of Noren (2006) 



## Example usage



```R

library(devtools) # see https://www.r-project.org/nosvn/pandoc/devtools.html
devtools::install_github('tystan/pharmsignal')
library(pharmsignal)

# The example values are taken from Table III (page 3749) of:
# Noren et al. (2006). Stats in Med, 25:3740–3757 (DOI: 10.1002/sim.2473).

n11 <- 25
n1. <- 1126
n.1 <- 87
n.. <- 572573

# convert to the corresponding a, b, c, d values for a 2x2 table
a <- n11
b <- n1. - n11
c <- n.1 - n11
d <- n.. - a - b - c

# perform calculations for values
set.seed(123456)
results_tab <- bcpnn_mcmc_signal(a, b, c, d)

# print in two parts as is a wide results table
n_cols <- ncol(results_tab)
half_cols <- ceiling(n_cols / 2)
results_tab[, 1:half_cols]

# analysis n11 drug margin event margin    n..     E_n11   est_name est_scale
#        1  25        1126           87 572573 0.1710908 bcpnn_mcmc      log2

results_tab[, (half_cols + 1):n_cols]

#      est var_scale var_est sd_est alpha   ci_lo   ci_hi
# 5.247846      <NA>      NA     NA  0.05 4.64927 5.73115


```


