# load CRAN packages
library(tidyverse)
library(ipumsr)
library(reticulate)
library(tictoc)

# load functions for the formally private quantile methods
source(here::here("R", "quantile_ind_exp.R"))
source(here::here("R", "quantile_joint_exp.R"))
source(here::here("R", "quantile_smooth.R"))

# load data
ddi <- read_ipums_ddi(here::here("data", "cps_00008.xml"))
data <- read_ipums_micro(ddi)

# preprocess data
asec <- data %>%
  filter(ASECFLAG == 1) %>%
  filter(RELATE %in% c("101", "201"))

asec <- asec %>%
  mutate(INCWAGE = zap_labels(INCWAGE))

# calculate the 0.1, 0.2, ..., 0.9 quantiles
probs <- seq(0.1, 0.9, 0.1)

# the data must be sorted for the methods to work
incwage <- sort(pull(asec, INCWAGE))

# parameters of interest
replicate <- 1:100
epsilon <- c(0.01, 0.5, 0.1, 0.5, 1, 5, 10, 15, 20)
delta <- c(10e-3, 10e-7, 10e-10)

# create grids of specifications
testing_grid <- expand_grid(
  epsilon,
  delta,
  replicate
)

# no noise ----------------------------------------------------------------
true_quantiles <- quantile(x = incwage, probs = probs)


# quantile_ind_exp --------------------------------------------------------

# quantile_ind_exp(
#   x = incwage,
#   probs = probs,
#   epsilon = 1,
#   delta = 1e-7,
#   data_low = min(incwage),
#   data_high = max(incwage)
# )

tic()
set.seed(894777)
quantile_ind_exp_test <- map2(
  .x = testing_grid$epsilon,
  .y = testing_grid$delta, 
  .f = ~quantile_ind_exp(
    x = incwage, 
    probs = probs, 
    epsilon = .x, 
    delta = .y, 
    data_low = min(incwage), 
    data_high = max(incwage)
  )
)
toc()

# quantile_joint_exp ------------------------------------------------------

# quantile_joint_exp(
#   x = incwage,
#   probs = probs,
#   epsilon = 1,
#   delta = 1e-7,
#   data_low = min(incwage),
#   data_high = max(incwage)
# )

tic()
set.seed(33633)
quantile_joint_exp_test <- map2(
  .x = testing_grid$epsilon,
  .y = testing_grid$delta,
  .f = ~quantile_joint_exp(
    x = incwage,
    probs = probs,
    epsilon = .x,
    delta = .y,
    data_low = min(incwage),
    data_high = max(incwage)
  )
)
toc()

# quantile_smooth ---------------------------------------------------------

# quantile_smooth(
#   x = incwage,
#   probs = probs,
#   epsilon = 1,
#   delta = 1e-7,
#   data_low = min(incwage),
#   data_high = max(incwage)
# )

tic()
set.seed(206994)
quantile_smooth_test <- map2(
  .x = testing_grid$epsilon,
  .y = testing_grid$delta, 
  .f = ~quantile_smooth(
    x = incwage, 
    probs = probs, 
    epsilon = .x, 
    delta = .y, 
    data_low = min(incwage), 
    data_high = max(incwage)
  )
)
toc()

# compile results ---------------------------------------------------------
calc_bias <- function(estimate, noisy_estimate) {
  
  mean(noisy_estimate - estimate)
  
}

calc_rmse <- function(estimate, noisy_estimate) {
  
  sqrt(mean((noisy_estimate - estimate) ^ 2))
  
}

# this stores results with one row per quantile
detailed_results <- bind_rows(
  quantile_ind_exp = tibble(
    testing_grid, 
    probs = list(probs),
    value = quantile_ind_exp_test
  ) %>%
    unnest(c(probs, value)),
  quantile_joint_exp = tibble(
    testing_grid, 
    probs = list(probs),
    value = quantile_joint_exp_test
  ) %>%
    unnest(c(probs, value)),
  quantile_smooth = tibble(
    testing_grid, 
    probs = list(probs),
    value = quantile_smooth_test
  ) %>%
    unnest(c(probs, value)),
  .id = "method"
) %>%
  mutate(value = as.numeric(value))

write_csv(detailed_results, here::here("results", "03_cps_quantiles-detailed.csv"))

# this stores results with one row per specification
results <- bind_rows(
  quantile_ind_exp = bind_cols(
    testing_grid,
    bias = map_dbl(quantile_ind_exp_test, ~calc_bias(estimate = true_quantiles, noisy_estimate = .x)),
    rmse = map_dbl(quantile_ind_exp_test, ~calc_rmse(estimate = true_quantiles, noisy_estimate = .x))
  ),
  quantile_joint_exp = bind_cols(
    testing_grid,
    bias = map_dbl(quantile_joint_exp_test, ~calc_bias(estimate = true_quantiles, noisy_estimate = .x)),
    rmse = map_dbl(quantile_joint_exp_test, ~calc_rmse(estimate = true_quantiles, noisy_estimate = .x))
  ),
  quantile_smooth = bind_cols(
    testing_grid,
    bias = map_dbl(quantile_smooth_test, ~calc_bias(estimate = true_quantiles, noisy_estimate = .x)),
    rmse = map_dbl(quantile_smooth_test, ~calc_rmse(estimate = true_quantiles, noisy_estimate = .x))
  ),
  .id = "method"
)

write_csv(results, here::here("results", "03_cps_quantiles.csv"))
