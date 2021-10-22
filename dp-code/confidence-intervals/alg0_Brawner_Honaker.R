######
# This file includes the code for the algorithm proposed by Brawner and Honaker, as
# discussed in Section 3.3 in the paper. 
######

library(tidyverse)

compute_epsilon <- function(rho, delta) {
  rho + 2 * sqrt(rho * log(sqrt(pi * rho) / delta))
}

inverse_function = function (f, delta, lower = 5e-05, upper = 250) {
  function(y){
    uniroot((function(x) f(x, delta) - y), lower = lower, upper = upper, tol = .Machine$double.eps ^ 0.5)[1][[1]]
  } 
}

#sapply(c(0.01, 0.5, 0.1, 0.5, 1, 2, 5, 10, 15, 20, 100), compute_rho)


#define_rho_i <- function(i, rho, N) {
#  (i * rho) / define_p_i(i, N) 
#}

## calls construct_boot_ci
#avg_cover_boot <- function(reps, n, epsilon, alpha, range) {
#  cov_vec <- rep(NA, reps)
#  moe_vec <- rep(NA, reps)
#  for (i in 1:reps) {
#    interval <- construct_boot_ci(rnorm(n), 50, epsilon, range,  alpha, alpha_prime = .05)
#    cov_vec[i] <- coverage(interval, 0)
#    moe_vec[i] <- diff(interval)
#  }
#  return(c(mean(cov_vec, na.rm = T), mean(moe_vec, na.rm = T)))

#}

#epsilon_rho_equiv = data.frame(epsilons = c(0.01, 0.5, 0.1, 0.5, 1, 2, 5, 10, 15, 20, 100), 
#                               delta = rep(1e-3, 11),
#                               rhos = c())

## fundamental
define_p_i <- function(i, N) {
  choose(N, i) * ((1 / N) ^ i) * (1 - 1 / N) ^ (N - i)
}

## fundamental
estimate_var <- function(k_boot_means, alpha_prime, N, rho, range) {
  sensitivity <- range / N
  k <- length(k_boot_means)
  c_a_prime <- qchisq(alpha_prime, k - 1)
  var(k_boot_means) - (sensitivity ^ 2) / (2 * rho) * ((k * c_a_prime) / (k - 1) - 1)
}

## calls define_p_i
define_sigma_i <- function(sensitivity, i, rho, N) {
  p_i <- define_p_i(i, N)
  i * p_i * (sensitivity ^ 2 / (2 * rho))
}

## calls define_sigma_i
calculate_partition_mean <- function(X_i, i, rho, range, N) {
  sensitivity <- range / N
  sigma_i <- define_sigma_i(sensitivity, i, rho, N)
  
  m_i <- i * sum(X_i) / N
  m_i + rnorm(1, 0, sqrt(sigma_i))
}

## fundamental
partition_bootstrap <- function(data) {
  N <- length(data)
  draws <- rmultinom(1, N, rep(1 / N, N))
  partitions <- list()
  for (i in 1:8) {
    partitions[[i]] <- rep(NA, round(N / 10))
  }
  for (i in 1:8) {
    partitions[[i]] <- data[draws == (i - 1)]
  }
  return(partitions)
}

## calls partition_bootstrap and caluclate_partition_mean
bootstrap_priv_mean <- function(data, rho, range) {
  N <- length(data)
  partitions <- partition_bootstrap(data)
  M_i_vec <- rep(0, N)
  for (i in 1:length(M_i_vec)) {
    if(i > 0 & i <= length(partitions)){
      M_i_vec[i] <- calculate_partition_mean(partitions[[i]], (i - 1), rho, range, N)
    } else{
      sensitivity <- range / N
      sigma_i <- define_sigma_i(sensitivity, (i - 1), rho, N)
      if(sigma_i == 'NaN'){
        M_i_vec[i] = 0
      } else{
        M_i_vec[i] = rnorm(1, 0, sqrt(sigma_i))
      }
    }
  }
  sum(M_i_vec)
}

## calls bootstrap_priv_mean and estimate_var
construct_boot_ci <- function(rep, data, k, epsilon, delta, range, alpha, alpha_prime) {
  compute_rho = inverse_function(compute_epsilon, delta = delta)
  rho = compute_rho(epsilon)
  #rho <- epsilon_rho_equiv$rhos[epsilon_rho_equiv$epsilons == epsilon]
  N <- length(data)
  boot_vec <- rep(NA, k) 
  
  for (i in 1:k) {
    boot_vec[i] <- bootstrap_priv_mean(data, rho/k, range)
  }
  
  mean_est <- mean(boot_vec)
  var_est <- max(0, estimate_var(boot_vec, alpha_prime, N, rho, range))
  se_est <- sqrt(var_est)
  z <- qnorm(1 - alpha / 2)
  c(mean_est - z * se_est, mean_est + z * se_est)
}


