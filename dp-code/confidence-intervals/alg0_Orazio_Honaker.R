######
# This file includes the code for the algorithm proposed by D'Orazio, Honaker and 
# King, as discussed in Section 3.2 in the paper. 
######

library(here)
library(rmutil)
library(tidyverse)
source(here("algorithms/alg5_EXPQ.R"))

#pnormlap <- function(y, mu, sigma, Lap_scale) {
#  alpha = 1 / Lap_scale
#  beta = 1 / Lap_scale

#  term1 <- pnorm((y - mu) / sigma)
#  term2 <- dnorm((y - mu) / sigma)
#  term3 <- beta * Mills_ratio(alpha * sigma - (y - mu) / sigma) - alpha * Mills_ratio(beta * sigma + (y - mu) / sigma)
#  term4 <- alpha + beta

#  term1 - ((term2 * term3) / term4)

#}

#dnormlap <- function(y, mu, sigma, Lap_scale) {
#  alpha = 1 / Lap_scale
#  beta = 1 / Lap_scale

#  term1 <- (alpha * beta) / (alpha + beta)
#  term2 <- dnorm((y - mu) / sigma)
#  term3 <- Mills_ratio(alpha * sigma - (y - mu) / sigma) + Mills_ratio(beta * sigma + (y - mu) / sigma)
#  term1 * term2 * term3
#}

#Mills_ratio <- function(z) {
#  (1 - pnorm(z, 0, 1)) / dnorm(z, 0, 1)
#} 


## fundamental
rnormlap <- function(n, mu, sigma, Lap_scale) {
  rnorm(n, mu, sigma) + rlaplace(n, 0, Lap_scale)
}

## fundamental
priv_se_orazio <- function(X, epsilon, SD_max) {
  n <- length(X)
  M <- round(n / 2)
  max <- SD_max / sqrt(n) + 3 * (SD_max / sqrt(2 * n/M)) * (1 / sqrt(n))
  X <- sample(X, replace = FALSE)
  S <- rep(NA, M)
  grp_size <- n / M
  for (i in 1:M) {
    S[i] <- sd(X[(grp_size * (i - 1) + 1):(grp_size * (i + 1))], na.rm = T) / sqrt(n)
  }
  a <- priv_median_c(S, 1/4, epsilon / 4, 0, max)
  b <- priv_median_c(S, 3/4, epsilon / 4, 0, max)
  mu <- (a + b) / 2
  iqr <- abs(a - b)
  u <- mu + 2 * iqr
  l <- mu - 2 * iqr
  S[S < l] <- l
  S[S > u] <- u
  w <- mean(S)
  Y <- w + rlaplace(1, 0, (u - l)/(.5 * epsilon * M))
  return(Y)
}

## calls priv_se_orazio, rnormlap
priv_interval_orazio <- function(db, a, epsilon, xmin, xmax, SD_max) {
  n <- length(db)
  Mean_est <- mean(db) + rlaplace(1, 0, (xmax - xmin) / (.5 * epsilon * n))
  SE <- priv_se_orazio(db, epsilon / 2, SD_max)
  sims <- rnormlap(10000, 0, max(SE, 0), (xmax - xmin) / (epsilon * .5 * n))
  b <- quantile(sims, 1 - a / 2)
  return(c(Mean_est - b, Mean_est + b))
}

