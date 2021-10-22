#load CRAN packages
library(tidyverse)
library(ipumsr)
library(tictoc)
library(furrr)

# load functions to implement the formally private mean methods
source(here::here("dp-code", "confidence-intervals", "alg0_Brawner_Honaker.R"))
source(here::here("dp-code", "confidence-intervals", "alg2_NOISYVAR.R"))
source(here::here("dp-code", "confidence-intervals", "alg4_NOISYMAD.R"))

# load metrics functions
source(here::here("R", "compute_confidence_interval_overlap.R"))

# use 14 cores with library(furrr)
# skip this line to run locally
plan(multisession, workers = 14)

# load data
ddi <- read_ipums_ddi(here::here("data", "cps_00008.xml"))
data <- read_ipums_micro(ddi)

# preprocess data
asec <- data %>%
  filter(ASECFLAG == 1) %>%
  filter(RELATE %in% c("101", "201"))

asec <- asec %>%
  mutate(INCWAGE = zap_labels(INCWAGE))

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
true_mean <- mean(x = asec$INCWAGE)

true_se <- sd(x = asec$INCWAGE) / sqrt(length(asec$INCWAGE))

true_ci <- c(true_mean - 1.96 * true_se, true_mean + 1.96 * true_se)

# noisyvar ----------------------------------------------------------------
#priv_ci(rep = 1, db = asec$INCWAGE, a = 0.05, e = 1, xmin = min(asec$INCWAGE), xmax = max(asec$INCWAGE), p = 0.8)

tic()
noisyvar_test <- future_map(
  .x = testing_grid$epsilon,
  .f = ~priv_ci(rep = 1, db = asec$INCWAGE, a = 0.05, e = .x, xmin = min(asec$INCWAGE), xmax = max(asec$INCWAGE), p = 0.8),
  .options = furrr::furrr_options(seed = 535931),
  .progress = TRUE
)
toc()

# noisymad ----------------------------------------------------------------
#priv_abs_sd_ci(rep = 1, db = asec$INCWAGE, a = 0.05, e = 1, xmin = min(asec$INCWAGE), xmax = max(asec$INCWAGE))

tic()
noisymad_test <- future_map(
  .x = testing_grid$epsilon,
  .f = ~priv_abs_sd_ci(rep = 1, db = asec$INCWAGE, a = 0.05, e = .x, xmin = min(asec$INCWAGE), xmax = max(asec$INCWAGE)),
  .options = furrr::furrr_options(seed = 927010),
  .progress = TRUE
)
toc()

# bhm ---------------------------------------------------------------------
# construct_boot_ci(rep = 1, 
#                   data = asec$INCWAGE, 
#                   k = 10, 
#                   epsilon = 1, 
#                   delta = 10e-7, 
#                   range = max(asec$INCWAGE), 
#                   alpha = 0.05, 
#                   alpha_prime = 0.05)

tic()
bhm_test <- future_map2(
  .x = bhm_testing_grid$epsilon,
  .y = bhm_testing_grid$delta,
  .f = ~construct_boot_ci(rep = 1, 
                          data = asec$INCWAGE, 
                          k = 10, 
                          epsilon = .x, 
                          delta = .y, 
                          range = max(asec$INCWAGE), 
                          alpha = 0.05, 
                          alpha_prime = 0.05),
  .options = furrr::furrr_options(seed = 90210),
  .progress = TRUE
)
toc()

# compile results ---------------------------------------------------------
calc_bias <- function(estimate, noisy_estimate) {
  
  mean(noisy_estimate) - mean(estimate) 
  
}

results <- bind_rows(
  NOISYVAR = bind_cols(
    testing_grid, 
    lower_bound = map_dbl(noisyvar_test, 1),
    upper_bound = map_dbl(noisyvar_test, 2),
    bias = map_dbl(noisyvar_test, ~calc_bias(estimate = true_ci, noisy_estimate = .x)),
    #rmse = sqrt(mean((true_mean - map_dbl(noisyvar_test, mean)) ^ 2)),
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

write_csv(results, here::here("results", "03_cps_mean-income.csv"))
