# load CRAN packages
library(tidyverse)

# load scripts that prep the data
source(here::here("R", "prep_cps03.R"))

# load scripts with the functions to implement the formally private methods
source(here::here("R", "lap_mech.R"))
source(here::here("R", "stat_inference.R"))
source(here::here("R", "count_inference.R"))
source(here::here("R", "laplace_sanitizer.R"))
source(here::here("R", "smart_round.R"))
source(here::here("R", "gauss_mech.R"))
source(here::here("R", "gauss_sanitizer.R"))
source(here::here("R", "gauss_renyi.R"))
source(here::here("R", "hist_laplace.R"))
source(here::here("R", "hist_laplace_multi.R"))
source(here::here("R", "hist_gaussian.R"))
source(here::here("R", "hist_gaussian_multi.R"))

# prepare the data
mortenson_cps <- prep_cps03()

# values for specifications
histogram_breaks <- seq(from = 0, to = 30000, by = 1000)

replicate <- 1:1000
epsilon <- c(0.01, 0.05, 0.1, 0.5, 1, 5, 10, 15, 20)
delta <- c(10e-3, 10e-7, 10e-10)

# no noise ----------------------------------------------------------------
true_histogram <- hist(mortenson_cps$earned_income_2014, breaks = histogram_breaks)$counts

# laplace mechanism -------------------------------------------------------

# hist_laplace(mortenson_cps$earned_income_2014, breaks = histogram_breaks, epsilon = 1)$counts

# hist_laplace_multi(mortenson_cps$earned_income_2014, breaks = histogram_breaks, epsilon = 1)$counts

# create a grid of specifications (laplace does not use delta)
laplace_testing_grid <- expand_grid(
  epsilon,
  replicate
)

# iterate tests
set.seed(1650456)
hist_laplace_test <- map(
  .x = laplace_testing_grid$epsilon,
  .f = ~hist_laplace(mortenson_cps$earned_income_2014, breaks = histogram_breaks, epsilon = .x)$counts
)

set.seed(909444)
hist_laplace_multi_test <- map(
  .x = laplace_testing_grid$epsilon,
  .f = ~hist_laplace_multi(mortenson_cps$earned_income_2014, breaks = histogram_breaks, epsilon = .x)$counts
)

# compile results
laplace_results <- 
  bind_rows(
    `Laplace Sanitizer` = bind_cols(
      laplace_testing_grid,
      bias = map_dbl(
        hist_laplace_test, .f = ~mean(.x - true_histogram)
      ),
      rmse = map_dbl(
        hist_laplace_test, .f = ~sqrt(mean((.x - true_histogram) ^ 2))
      )
    ),
    `Laplace Multiple Queries` = bind_cols(  
      laplace_testing_grid,
      bias = map_dbl(
        hist_laplace_multi_test, .f = ~mean(.x - true_histogram)
      ),
      rmse = map_dbl(
        hist_laplace_multi_test, .f = ~sqrt(mean((.x - true_histogram) ^ 2))
      )
    ),
    .id = "method"
  )

# gaussian mechanism ------------------------------------------------------

# run method once
# hist_gaussian(mortenson_cps$earned_income_2014, breaks = histogram_breaks, epsilon = 1, delta = 1e-5)$counts

# hist_gaussian_multi(mortenson_cps$earned_income_2014, breaks = histogram_breaks, epsilon = 1, delta = 1e-5)$counts

# create a grid of specifications
gaussian_testing_grid <- expand_grid(
  epsilon,
  replicate,
  delta
)

# iterate tests
set.seed(3278686)
hist_gaussian_test <- map2(
  .x = gaussian_testing_grid$epsilon,
  .y = gaussian_testing_grid$delta,
  .f = ~hist_gaussian(mortenson_cps$earned_income_2014, breaks = histogram_breaks, epsilon = .x, delta = .y)$counts
)

set.seed(993876)
hist_gaussian_multi_test <- map2(
  .x = gaussian_testing_grid$epsilon,
  .y = gaussian_testing_grid$delta,
  .f = ~hist_gaussian_multi(mortenson_cps$earned_income_2014, breaks = histogram_breaks, epsilon = .x, delta = .y)$counts
)

# compile results
gaussian_results <- 
  bind_rows(
    `Gaussian Mechanism` = bind_cols(
      gaussian_testing_grid,
      bias = map_dbl(
        hist_gaussian_test, .f = ~mean(.x - true_histogram)
      ),
      rmse = map_dbl(
        hist_gaussian_test, .f = ~sqrt(mean((.x - true_histogram) ^ 2)),
      )
    ),
    `Gaussian Multiple Queries` = bind_cols(
      gaussian_testing_grid,
      bias = map_dbl(
        hist_gaussian_multi_test, .f = ~mean(.x - true_histogram)
      ),
      rmse = map_dbl(
        hist_gaussian_multi_test, .f = ~sqrt(mean((.x - true_histogram) ^ 2))
      )
    ),
    .id = "method"
  )
  
# combine laplacian and gaussian results and save the summarized results
results <- bind_rows(
  laplace_results,
  gaussian_results
)

write_csv(results, here::here("results", "03_cps_histograms.csv"))

# combine laplacian and gaussian results and save the detailed results
detailed_results <- bind_rows(
  `Laplace Sanitizer` = tibble(
    laplace_testing_grid,
    bin_id = list(1:length(true_histogram)),
    noisy_count = hist_laplace_test
  ),
  `Laplace Multiple Queries` = tibble(
    laplace_testing_grid,
    bin_id = list(1:length(true_histogram)),
    noisy_count = hist_laplace_multi_test
  ),
  `Gaussian Mechanism` = tibble(
    gaussian_testing_grid,
    bin_id = list(1:length(true_histogram)),
    noisy_count = hist_gaussian_test
  ),
  `Gaussian Multiple Queries`  = tibble(
    gaussian_testing_grid,
    bin_id = list(1:length(true_histogram)),
    noisy_count = hist_gaussian_multi_test
  ),
  .id = "method"
) %>%
  unnest(cols = c(bin_id, noisy_count)) %>%
  select(method, replicate, epsilon, delta, bin_id, noisy_count)

write_csv(detailed_results, here::here("results", "03_cps_histograms_detailed.csv"))
