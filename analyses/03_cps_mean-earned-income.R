# load CRAN packages
library(tidyverse)
library(ipumsr)
library(tictoc)
library(furrr)
library(reticulate)

# load functions to implement the formally private mean methods
# we use methods here because the data are symmetrical
source(here::here('dp-code', 'confidence-intervals', 'alg0_Brawner_Honaker.R'))
source(here::here('dp-code', 'confidence-intervals', 'alg2_NOISYVAR.R'))
source(here::here('dp-code', 'confidence-intervals', 'alg4_NOISYMAD.R'))
source(here::here('dp-code', 'confidence-intervals', 'alg6_CENQ_python.R'))
source(here::here('dp-code', 'confidence-intervals', 'alg7_SYMQ_python.R'))
source(here::here('dp-code', 'confidence-intervals', 'alg8_MOD_python.R'))

# load metrics functions
source(here::here("R", "compute_confidence_interval_overlap.R"))

# load function to prep the data
source(here::here("R", "prep_cps03.R"))

# use 14 cores with library(furrr)
# skip this line to run locally
plan(multisession, workers = 14)

# turn earned income into a vector
earned_income <- prep_cps03() %>%
  pull(earned_income_2014)

# parameters of interest
replicate <- 1:100
epsilon <- c(0.01, 0.5, 0.1, 0.5, 1, 5, 10, 15, 20)
delta <- c(10e-3, 10e-7, 10e-10)

# create grids of specifications
testing_grid <- expand_grid(
  epsilon,
  replicate
)

bhm_testing_grid <- expand_grid(
  epsilon,
  delta,
  replicate
) %>%
  filter(epsilon > 0.01) # 0.01 is too small for bhm

# no noise ----------------------------------------------------------------
true_mean <- mean(x = earned_income)

true_se <- sd(x = earned_income) / sqrt(length(earned_income))

true_ci <- c(true_mean - qnorm(0.975) * true_se, true_mean + qnorm(0.975) * true_se)

# SYMQ --------------------------------------------------------------------
# priv_double_ci(rep = 1, db = earned_income, a = 0.05, e = 1, xmin = min(earned_income), xmax = max(earned_income))

tic()
set.seed(797128)
symq_test <- map(
  .x = testing_grid$epsilon,
  .f = ~priv_double_ci(rep = 1, db = earned_income, a = 0.05, e = .x, xmin = min(earned_income), xmax = max(earned_income))
)
toc()

# CENQ --------------------------------------------------------------------
# priv_exp_ci(rep = 1, db = earned_income, a = 0.05, e = 1, xmin = min(earned_income), xmax = max(earned_income))

tic()
set.seed(781033)
cenq_test <- map(
  .x = testing_grid$epsilon,
  .f = ~priv_exp_ci(rep = 1, db = earned_income, a = 0.05, e = .x, xmin = min(earned_income), xmax = max(earned_income))
)
toc()

# MOD ---------------------------------------------------------------------
# priv_exp_ci_abs(rep = 1, db = earned_income, a = 0.05, e = 1, xmin = min(earned_income), xmax = max(earned_income))

tic()
set.seed(892580)
mod_test <- map(
  .x = testing_grid$epsilon,
  .f = ~priv_exp_ci_abs(rep = 1, db = earned_income, a = 0.05, e = .x, xmin = min(earned_income), xmax = max(earned_income))
)
toc()

# noisyvar ----------------------------------------------------------------
# priv_ci(rep = 1, db = earned_income, a = 0.05, e = 1, xmin = min(earned_income), xmax = max(earned_income), p = 0.8)

tic()
noisyvar_test <- future_map(
  .x = testing_grid$epsilon,
  .f = ~priv_ci(rep = 1, db = earned_income, a = 0.05, e = .x, xmin = min(earned_income), xmax = max(earned_income), p = 0.8),
  .options = furrr::furrr_options(seed = 535931),
  .progress = TRUE
)
toc()

# noisymad ----------------------------------------------------------------
# priv_abs_sd_ci(rep = 1, db = earned_income, a = 0.05, e = 1, xmin = min(earned_income), xmax = max(earned_income))

tic()
noisymad_test <- future_map(
  .x = testing_grid$epsilon,
  .f = ~priv_abs_sd_ci(rep = 1, db = earned_income, a = 0.05, e = .x, xmin = min(earned_income), xmax = max(earned_income)),
  .options = furrr::furrr_options(seed = 927010),
  .progress = TRUE
)
toc()

# bhm ---------------------------------------------------------------------
construct_boot_ci(rep = 1,
                  data = earned_income,
                  k = 10,
                  epsilon = 0.1,
                  delta = 10e-10,
                  range = max(earned_income),
                  alpha = 0.05,
                  alpha_prime = 0.05)

tic()
bhm_test <- future_map2(
  .x = bhm_testing_grid$epsilon,
  .y = bhm_testing_grid$delta,
  .f = ~construct_boot_ci(rep = 1, 
                          data = earned_income, 
                          k = 10, 
                          epsilon = .x, 
                          delta = .y, 
                          range = max(earned_income), 
                          alpha = 0.05, 
                          alpha_prime = 0.05),
  .options = furrr::furrr_options(seed = 363574),
  .progress = TRUE
)
toc()

# compile results ---------------------------------------------------------
calc_bias <- function(estimate, noisy_estimate) {
  
  mean(noisy_estimate) - mean(estimate) 
  
}

results <- bind_rows(
  SYMQ = bind_cols(
    testing_grid, 
    lower_bound = map_dbl(symq_test, 1),
    upper_bound = map_dbl(symq_test, 2),
    bias = map_dbl(symq_test, ~calc_bias(estimate = true_ci, noisy_estimate = .x)),
    ci_overlap = map_dbl(symq_test, ~compute_confidence_interval_overlap(
      private_lower = true_ci[1], 
      private_upper = true_ci[2], 
      nonprivate_lower = .x[1], 
      nonprivate_upper = .x[2]
      )
    )
  ),
  CENQ = bind_cols(
    testing_grid, 
    lower_bound = map_dbl(cenq_test, 1),
    upper_bound = map_dbl(cenq_test, 2),
    bias = map_dbl(cenq_test, ~calc_bias(estimate = true_ci, noisy_estimate = .x)),
    ci_overlap = map_dbl(cenq_test, ~compute_confidence_interval_overlap(
      private_lower = true_ci[1], 
      private_upper = true_ci[2], 
      nonprivate_lower = .x[1], 
      nonprivate_upper = .x[2]
      )
    )
  ),
  MOD = bind_cols(
    testing_grid, 
    lower_bound = map_dbl(mod_test, 1),
    upper_bound = map_dbl(mod_test, 2),
    bias = map_dbl(mod_test, ~calc_bias(estimate = true_ci, noisy_estimate = .x)),
    ci_overlap = map_dbl(mod_test, ~compute_confidence_interval_overlap(
      private_lower = true_ci[1], 
      private_upper = true_ci[2], 
      nonprivate_lower = .x[1], 
      nonprivate_upper = .x[2]
      )
    )
  ),
  NOISYVAR = bind_cols(
    testing_grid, 
    lower_bound = map_dbl(noisyvar_test, 1),
    upper_bound = map_dbl(noisyvar_test, 2),
    bias = map_dbl(noisyvar_test, ~calc_bias(estimate = true_ci, noisy_estimate = .x)),
    ci_overlap = map_dbl(noisyvar_test, ~compute_confidence_interval_overlap(
      private_lower = true_ci[1], 
      private_upper = true_ci[2], 
      nonprivate_lower = .x[1], 
      nonprivate_upper = .x[2]
      )
    )
  ),
  NOISYMAD = bind_cols(
    testing_grid, 
    lower_bound = map_dbl(noisymad_test, 1),
    upper_bound = map_dbl(noisymad_test, 2),
    bias = map_dbl(noisymad_test, ~calc_bias(estimate = true_ci, noisy_estimate = .x)),
    ci_overlap = map_dbl(noisymad_test, ~compute_confidence_interval_overlap(
      private_lower = true_ci[1], 
      private_upper = true_ci[2], 
      nonprivate_lower = .x[1], 
      nonprivate_upper = .x[2]
      )
    )
  ),
  bhm = bind_cols(
    bhm_testing_grid,
    lower_bound = map_dbl(bhm_test, 1),
    upper_bound = map_dbl(bhm_test, 2),
    bias = map_dbl(bhm_test, ~calc_bias(estimate = true_ci, noisy_estimate = .x)),
    ci_overlap = map_dbl(bhm_test, ~compute_confidence_interval_overlap(
      private_lower = true_ci[1], 
      private_upper = true_ci[2], 
      nonprivate_lower = .x[1], 
      nonprivate_upper = .x[2]
      )
    )
  ),
  .id = "method"
)

write_csv(results, here::here("results", "03_cps_mean-earned-income.csv"))
