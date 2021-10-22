library(tidyverse)
library(LaplacesDemon)
library(broom)
library(tictoc)
library(furrr)

# use 34 cores with library(furrr)
# skip this line to run locally
plan(multisession, workers = 34)

# load the function to prep the data
source(here::here("R", "prep_cps04.R"))

# load functions for formally private regression
source(here::here("dp-code", "DPregression", "mechanisms", "DPBootsLMFunctions.R"))
source(here::here("dp-code", "DPregression", "mechanisms", "AnalyticGaussianMechanism.R"))
source(here::here("dp-code", "DPregression", "mechanisms", "DPBootsAGMFunctions.R"))
source(here::here("dp-code", "DPregression", "mechanisms", "DPBootsWMFunctions_new.R"))
source(here::here("dp-code", "DPregression", "mechanisms", "DPBootsBHMFunctions.R"))
source(here::here("dp-code", "DPregression", "mechanisms", "DPBootsWKLMFunctions.R"))

# load wrapper functions for formally private regression
source(here::here("R", "lm_laplace.R"))
source(here::here("R", "lm_analytic_gauss.R"))
source(here::here("R", "lm_wishart.R"))
source(here::here("R", "lm_bhm.R"))
source(here::here("R", "lm_wklm.R"))

# metrics functions
source(here::here("R", "summarize_lm.R"))
source(here::here("R", "calc_lm_rmse.R"))
source(here::here("R", "compute_confidence_interval_overlap.R"))

# load the data
asec_male <- prep_cps04()$male %>%
  select(
    years_of_educ, 
    potential_experience, 
    potential_experience_squared, 
    potential_experience_cubed, 
    non_white, 
    log_INCWAGE
  )

# create grids of specifications ------------------------------------------

# values for specifications
replicate <- 1:100
epsilon <- c(0.5, 1, 5, 10, 15, 20, 1e6)
epsilon_wishart <- c(0.5, 0.99999999999)
delta <- c(1, 10e-3, 10e-7, 10e-10)
bootstrap_size <- c(10, 25)

# laplace doesn't use delta or n
epsilon_only_testing_grid <- expand_grid(
  epsilon,
  replicate
)

# analytic doesn't us n
analytic_gauss_testing_grid <- expand_grid(
  replicate,
  epsilon,
  delta
)

# wishart only holds when epsilon < 1
wishart_testing_grid <- expand_grid(
  replicate,
  epsilon = epsilon_wishart,
  delta
)

bhm_testing_grid <- expand_grid(
  replicate,
  epsilon,
  delta,
  n = bootstrap_size
)

# wklm doesn't use n
wklm_testing_grid <- expand_grid(
  replicate,
  epsilon,
  delta
)

# inputs for regression methods -------------------------------------------
model_formula <- log_INCWAGE ~ 
  years_of_educ + 
  potential_experience + 
  potential_experience_squared + 
  potential_experience_cubed + 
  non_white

bounds <- list(
  years_of_educ = c(min(asec_male$years_of_educ), max(asec_male$years_of_educ)), 
  potential_experience = c(min(asec_male$potential_experience), max(asec_male$potential_experience)),
  potential_experience_squared = c(min(asec_male$potential_experience_squared), max(asec_male$potential_experience_squared)),
  potential_experience_cubed = c(min(asec_male$potential_experience_cubed), max(asec_male$potential_experience_cubed)),
  non_white = c("0", "1"),
  log_INCWAGE = c(min(asec_male$log_INCWAGE), max(asec_male$log_INCWAGE))
)

reference_classes <- list(non_white = c("0"))

# List of variables (last one is the response)
var_names <- c("years_of_educ", 
               "potential_experience", 
               "potential_experience_squared", 
               "potential_experience_cubed",
               "non_white",
               "log_INCWAGE") # predictors name

# Type of variables
var_types <- c(years_of_educ = "numeric", 
               potential_experience = "numeric", 
               potential_experience_squared  = "numeric", 
               potential_experience_cubed = "numeric", 
               non_white = "factor",
               log_INCWAGE = "numeric") # predictors type

# no noise ----------------------------------------------------------------
true_lm <- lm(log_INCWAGE ~ 
                years_of_educ + 
                potential_experience + 
                potential_experience_squared + 
                potential_experience_cubed + 
                non_white, 
              data = asec_male) %>%
  tidy(conf.int = TRUE)

# laplace mechanism -------------------------------------------------------

# lm_laplace(
#   formula = model_formula,
#   data = asec_male,
#   epsilon = 1e6,
#   bounds = bounds,
#   var_names = var_names,
#   var_types = var_types,
#   reference_classes = reference_classes
# )

tic()
lm_laplace_test <- future_map(
  .x = epsilon_only_testing_grid$epsilon,
  .f = ~lm_laplace(
    formula = model_formula, 
    data = asec_male, 
    epsilon = .x, 
    bounds = bounds,
    var_names = var_names,
    var_types = var_types,
    reference_classes = reference_classes
  ),
  .options = furrr::furrr_options(seed = 559506),
  .progress = TRUE
)
toc()

# analytic gaussian mechanism ---------------------------------------------

# lm_analytic_gauss(
#   formula = model_formula,
#   data = asec_male,
#   epsilon = 1e6,
#   delta = 1,
#   bounds = bounds,
#   var_names = var_names,
#   var_types = var_types,
#   reference_classes = reference_classes
# )

tic()
lm_analytic_gauss_test <- future_map2(
  .x = analytic_gauss_testing_grid$epsilon,
  .y = analytic_gauss_testing_grid$delta, 
  .f = ~lm_analytic_gauss(
    formula = model_formula, 
    data = asec_male, 
    epsilon = .x, 
    delta = .y,
    bounds = bounds,
    var_names = var_names,
    var_types = var_types,
    reference_classes = reference_classes
  ),
  .options = furrr::furrr_options(seed = 330486),
  .progress = TRUE
)
toc()

# wishart mechanism -------------------------------------------------------

# tic()
# lm_wishart(
#   formula = model_formula,
#   data = as.data.frame(asec_male),
#   epsilon = 1,
#   delta = 0.1,
#   bounds = bounds,
#   var_names = var_names,
#   var_types = var_types,
#   reference_classes = reference_classes
# )
# toc()

tic()
lm_wishart_test <- future_map2(
  .x = wishart_testing_grid$epsilon,
  .y = wishart_testing_grid$delta,
  .f = ~lm_wishart(
    formula = model_formula,
    data = as.data.frame(asec_male),
    epsilon = .x,
    delta = .y,
    bounds = bounds,
    var_names = var_names,
    var_types = var_types,
    reference_classes = reference_classes
  ),
  .options = furrr::furrr_options(seed = 148174),
  .progress = TRUE
)
toc()

# brawner and hanaker method ----------------------------------------------

# lm_bhm(
#   formula = model_formula,
#   data = asec_male,
#   epsilon = 1e6,
#   delta = 1,
#   bounds = bounds,
#   var_names = var_names,
#   var_types = var_types,
#   reference_classes = reference_classes,
#   n = 10
# )

tic()
lm_bhm_test <- future_pmap(
  .l = select(bhm_testing_grid, epsilon, delta, n),
  .f = function(epsilon, delta, n) lm_bhm(
    formula = model_formula, 
    data = asec_male, 
    epsilon = epsilon, 
    delta = delta,
    bounds = bounds,
    var_names = var_names,
    var_types = var_types,
    reference_classes = reference_classes,
    n = n
  ),
  .options = furrr::furrr_options(seed = 991191),
  .progress = TRUE
)
toc()

# Wang,  Kifer,  and Lee method -------------------------------------------

# tic()
# lm_wklm(
#   formula = model_formula,
#   data = as.data.frame(asec_male),
#   epsilon = 10,
#   delta = 10e-7,
#   bounds = bounds,
#   var_names = var_names,
#   var_types = var_types,
#   reference_classes = reference_classes,
#   dp_type = "epsilon-DP",
# )
# toc()

tic()
lm_wklm_pure_test <- future_map(
  .x = epsilon_only_testing_grid$epsilon,
  .f = ~lm_wklm(
    formula = model_formula, 
    data = as.data.frame(asec_male), 
    epsilon = .x, 
    delta = NULL,
    bounds = bounds,
    var_names = var_names,
    var_types = var_types,
    dp_type = "epsilon-DP",
    reference_classes = reference_classes
  ),
  .options = furrr::furrr_options(seed = 336863),
  .progress = TRUE
)
toc()

# tic()
# lm_wklm(
#   formula = model_formula,
#   data = as.data.frame(asec_male),
#   epsilon = 10,
#   delta = 10e-7,
#   bounds = bounds,
#   var_names = var_names,
#   var_types = var_types,
#   reference_classes = reference_classes,
#   dp_type = "rho-zCDP",
# )
# toc()

tic()
lm_wklm_ed_test <- future_map2(
  .x = wklm_testing_grid$epsilon,
  .y = wklm_testing_grid$delta, 
  .f = ~lm_wklm(
    formula = model_formula, 
    data = as.data.frame(asec_male), 
    epsilon = .x, 
    delta = .y,
    bounds = bounds,
    var_names = var_names,
    var_types = var_types,
    dp_type = "rho-zCDP",
    reference_classes = reference_classes
  ),
  .options = furrr::furrr_options(seed = 290615),
  .progress = TRUE
)
toc()

# compile results ---------------------------------------------------------

# asymptotic results at the coefficient level
coefficient_results_asymptotic <- bind_rows(
  `Laplace Mechanism` = bind_cols(
    expand_grid(epsilon_only_testing_grid, term = lm_laplace_test[[1]]$asymptotic$term) %>% select(-term),
    map_df(
      lm_laplace_test, .f = ~summarize_lm(results = true_lm, noisy_results = .x$asymptotic)
    )
  ),
  `Analytic Gaussian Mechanism` = bind_cols(
    expand_grid(analytic_gauss_testing_grid, term = lm_analytic_gauss_test[[1]]$asymptotic$term) %>% select(-term),
    map_df(
      lm_analytic_gauss_test, .f = ~summarize_lm(results = true_lm, noisy_results = .x$asymptotic)
    )          
  ),
  `Wishart Mechanism` = bind_cols(
    expand_grid(wishart_testing_grid, term = lm_wishart_test[[1]]$asymptotic$term) %>% select(-term),
    map_df(
      lm_wishart_test, .f = ~summarize_lm(results = true_lm, noisy_results = .x$asymptotic)
    )  
  ),
  `Brawner and Hanaker Method` = bind_cols(
    expand_grid(bhm_testing_grid, term = lm_bhm_test[[1]]$asymptotic$term) %>% select(-term),
    map_df(
      lm_bhm_test, .f = ~summarize_lm(results = true_lm, noisy_results = .x$asymptotic)
    )  
  ),
  `WKL Mechanism (Pure)` = bind_cols(
    expand_grid(epsilon_only_testing_grid, term = lm_wklm_pure_test[[1]]$asymptotic$term) %>% select(-term),
    map_df(
      lm_wklm_pure_test, .f = ~summarize_lm(results = true_lm, noisy_results = .x$asymptotic)
    )  
  ),
  `WKL Mechanism (ED)` = bind_cols(
    expand_grid(wklm_testing_grid, term = lm_wklm_ed_test[[1]]$asymptotic$term) %>% select(-term),
    map_df(
      lm_wklm_ed_test, .f = ~summarize_lm(results = true_lm, noisy_results = .x$asymptotic)
    )  
  ),  
  .id = "method"
)

write_csv(coefficient_results_asymptotic, file = here::here("results", "04_cps_male-regression_coefficient_asymptotic.csv"))

# bootstrap results at the coefficient level
coefficient_results_bootstrap <- bind_rows(
  `Laplace Mechanism` = bind_cols(
    expand_grid(epsilon_only_testing_grid, term = lm_laplace_test[[1]]$bootstrap$term) %>% select(-term),
    map_df(
      lm_laplace_test, .f = ~summarize_lm(results = true_lm, noisy_results = .x$bootstrap)
    )
  ),
  `Analytic Gaussian Mechanism` = bind_cols(
    expand_grid(analytic_gauss_testing_grid, term = lm_analytic_gauss_test[[1]]$bootstrap$term) %>% select(-term),
    map_df(
      lm_analytic_gauss_test, .f = ~summarize_lm(results = true_lm, noisy_results = .x$bootstrap)
    )          
  ),
  `Wishart Mechanism` = bind_cols(
    expand_grid(wishart_testing_grid, term = lm_wishart_test[[1]]$bootstrap$term) %>% select(-term),
    map_df(
      lm_wishart_test, .f = ~summarize_lm(results = true_lm, noisy_results = .x$bootstrap)
    )  
  ),
  `Brawner and Hanaker Method` = bind_cols(
    expand_grid(bhm_testing_grid, term = lm_bhm_test[[1]]$bootstrap$term) %>% select(-term),
    map_df(
      lm_bhm_test, .f = ~summarize_lm(results = true_lm, noisy_results = .x$bootstrap)
    )  
  ),
  `WKL Mechanism (Pure)` = bind_cols(
    expand_grid(epsilon_only_testing_grid, term = lm_wklm_pure_test[[1]]$bootstrap$term) %>% select(-term),
    map_df(
      lm_wklm_pure_test, .f = ~summarize_lm(results = true_lm, noisy_results = .x$bootstrap)
    )  
  ),
  `WKL Mechanism (ED)` = bind_cols(
    expand_grid(wklm_testing_grid, term = lm_wklm_ed_test[[1]]$bootstrap$term) %>% select(-term),
    map_df(
      lm_wklm_ed_test, .f = ~summarize_lm(results = true_lm, noisy_results = .x$bootstrap)
    )  
  ),
  .id = "method"
)

write_csv(coefficient_results_bootstrap, file = here::here("results", "04_cps_male-regression_coefficient_bootstrap.csv"))

# asymptotic results at the model level
rmse_results_asymptotic <- bind_rows(
  `Laplace Mechanism` = bind_cols(
    epsilon_only_testing_grid, 
    rmse = map_dbl(
      lm_laplace_test, .f = ~calc_lm_rmse(
        data = asec_male, 
        formula = model_formula,
        truth = asec_male$log_INCWAGE,
        noisy_results = .x$asymptotic
      )
    )
  ),
  `Analytic Gaussian Mechanism` = bind_cols(
    analytic_gauss_testing_grid,
    rmse = map_dbl(
      lm_analytic_gauss_test, .f = ~calc_lm_rmse(
        data = asec_male, 
        formula = model_formula,
        truth = asec_male$log_INCWAGE,
        noisy_results = .x$asymptotic
      )
    )          
  ),
  `Wishart Mechanism` = bind_cols(
    wishart_testing_grid, 
    rmse = map_dbl(
      lm_wishart_test, .f = ~calc_lm_rmse(
        data = asec_male, 
        formula = model_formula,
        truth = asec_male$log_INCWAGE,
        noisy_results = .x$asymptotic
      )
    )  
  ),
  `Brawner and Hanaker Method` = bind_cols(
    bhm_testing_grid,
    rmse = map_dbl(
      lm_bhm_test, .f = ~calc_lm_rmse(
        data = asec_male, 
        formula = model_formula,
        truth = asec_male$log_INCWAGE,
        noisy_results = .x$asymptotic
      )
    )  
  ),
  `WKL Mechanism (Pure)` = bind_cols(
    epsilon_only_testing_grid,
    rmse = map_dbl(
      lm_wklm_pure_test, .f = ~calc_lm_rmse(
        data = asec_male, 
        formula = model_formula,
        truth = asec_male$log_INCWAGE,
        noisy_results = .x$asymptotic
      )
    )  
  ),
  `WKL Mechanism (ED)` = bind_cols(
    wklm_testing_grid,
    rmse = map_dbl(
      lm_wklm_ed_test, .f = ~calc_lm_rmse(
        data = asec_male, 
        formula = model_formula,
        truth = asec_male$log_INCWAGE,
        noisy_results = .x$asymptotic
      )
    )  
  ),  
  .id = "method"
)

write_csv(rmse_results_asymptotic, file = here::here("results", "04_cps_male-regression_rmse_asymptotic.csv"))

# bootstrap results at the model level
rmse_results_bootstrap <- bind_rows(
  `Laplace Mechanism` = bind_cols(
    epsilon_only_testing_grid, 
    rmse = map_dbl(
      lm_laplace_test, .f = ~calc_lm_rmse(
        data = asec_male, 
        formula = model_formula,
        truth = asec_male$log_INCWAGE,
        noisy_results = .x$bootstrap
      )
    )
  ),
  `Analytic Gaussian Mechanism` = bind_cols(
    analytic_gauss_testing_grid,
    rmse = map_dbl(
      lm_analytic_gauss_test, .f = ~calc_lm_rmse(
        data = asec_male, 
        formula = model_formula,
        truth = asec_male$log_INCWAGE,
        noisy_results = .x$bootstrap
      )
    )          
  ),
  `Wishart Mechanism` = bind_cols(
    wishart_testing_grid, 
    rmse = map_dbl(
      lm_wishart_test, .f = ~calc_lm_rmse(
        data = asec_male, 
        formula = model_formula,
        truth = asec_male$log_INCWAGE,
        noisy_results = .x$bootstrap
      )
    )  
  ),
  `Brawner and Hanaker Method` = bind_cols(
    bhm_testing_grid,
    rmse = map_dbl(
      lm_bhm_test, .f = ~calc_lm_rmse(
        data = asec_male, 
        formula = model_formula,
        truth = asec_male$log_INCWAGE,
        noisy_results = .x$bootstrap
      )
    )  
  ),
  `WKL Mechanism (Pure)` = bind_cols(
    epsilon_only_testing_grid,
    rmse = map_dbl(
      lm_wklm_pure_test, .f = ~calc_lm_rmse(
        data = asec_male, 
        formula = model_formula,
        truth = asec_male$log_INCWAGE,
        noisy_results = .x$bootstrap
      )
    )  
  ),
  `WKL Mechanism (ED)` = bind_cols(
    wklm_testing_grid,
    rmse = map_dbl(
      lm_wklm_ed_test, .f = ~calc_lm_rmse(
        data = asec_male, 
        formula = model_formula,
        truth = asec_male$log_INCWAGE,
        noisy_results = .x$bootstrap
      )
    )  
  ),  
  .id = "method"
)

write_csv(rmse_results_bootstrap, file = here::here("results", "04_cps_male-regression_rmse_bootstrap.csv"))
