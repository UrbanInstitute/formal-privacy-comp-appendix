library(tidyverse)
library(ggridges)
library(ggpol)

theme_set(theme_minimal())
theme_update(plot.title.position = "plot")

#####
### Histograms
#####
histogram_results_summary = readRDS('report_figures/cps_histogram_summaries.rds')

histogram_results_plot_1 = histogram_results_summary %>%
  filter(method %in% c('Gaussian Mechanism', 'Laplace Sanitizer'),
         measure %in% c('max_error', 'max_error_rev')) %>%
  ## take max over relative both directions
  group_by(method, epsilon, delta, replicate) %>%
  mutate(max_error_both = pmax(value[measure == 'max_error'], value[measure == 'max_error_rev'])) %>%
  ## drop duplicate
  filter(measure == 'max_error') %>%
  ungroup()


png('report_figures/cps/income_hist_max_relative.png', width = 1600, height = 900)
ggplot(histogram_results_plot_1 %>%
         filter(epsilon >= 0.1, 
                epsilon <= 5) %>%
         mutate(epsilon_char = recode(epsilon, 
                                      `0.1` = 'Epsilon = 0.1',
                                      `0.5` = 'Epsilon = 0.5',
                                      `1` = 'Epsilon = 1.0',
                                      `5` = 'Epsilon = 5.0')), 
       aes(x = max_error_both, y = method, fill = as.factor(delta))) +
  geom_violin(draw_quantiles = c(0.5)) +
  facet_wrap(~epsilon_char, dir = 'v', scales = 'free_x') +
  scale_fill_manual(name = 'Delta', 
                    values = c('#e41a1c', '#377eb8', '#4daf4a'), 
                    na.value = 'grey50',
                    labels = c('10e-10', '10e-07', '10e-03')) +
  theme_minimal(base_size = 24) +
  scale_x_continuous(labels = scales::percent) +
  ylab('') +
  xlab('Max Error (%) of Cumulative Sums Relative to True Total')
dev.off()  

## get RMSE
#histogram_results = read_csv(here::here("results", "01_cps_histograms.csv"), guess_max = 36000)
histogram_results_plot_2 = histogram_results_summary %>%
  filter(method %in% c('Gaussian Mechanism', 'Laplace Sanitizer'),
         measure %in% c('mean_relative_error', 'mean_relative_error_rev')) %>%
  ## take max over relative both directions
  group_by(method, epsilon, delta, replicate) %>%
  mutate(mean_relative_error_both = mean(value)) %>%
  ## drop duplicate
  filter(measure == 'mean_relative_error') %>%
  ungroup()

png('report_figures/cps/income_hist_mean_cumulative.png', width = 1600, height = 900)
ggplot(histogram_results_plot_2 %>%
         filter(method %in% c('Gaussian Mechanism', 'Laplace Sanitizer'),
                epsilon >= 0.1, 
                epsilon <= 5) %>%
         mutate(epsilon_char = recode(epsilon, 
                                      `0.1` = 'Epsilon = 0.1',
                                      `0.5` = 'Epsilon = 0.5',
                                      `1` = 'Epsilon = 1.0',
                                      `5` = 'Epsilon = 5.0')), 
       aes(x = mean_relative_error_both, y = method, fill = as.factor(delta))) +
  geom_violin(draw_quantiles = c(0.5)) +
  facet_wrap(~epsilon_char, dir = 'v', scales = 'free_x') +
  scale_fill_manual(name = 'Delta', 
                    values = c('#e41a1c', '#377eb8', '#4daf4a'), 
                    na.value = 'grey50',
                    labels = c('10e-10', '10e-07', '10e-03')) +
  scale_x_continuous(labels = scales::percent) +
  theme_minimal(base_size = 24) +
  ylab('') +
  xlab('Mean Error (%) of Cumulative Totals')
dev.off()  



#####
### Means
#####
mean_income_results = read_csv(here::here("results", "03_cps_mean-income.csv"), guess_max = 4500)

mean_earned_income_results = read_csv(here::here("results", "03_cps_mean-earned-income.csv"), guess_max = 4500)

## add some values
mean_income_results = mean_income_results %>%
  mutate(dp_mean = (upper_bound - lower_bound) / 2 + lower_bound,
         true_mean = dp_mean - bias,
         true_upper = 32826.76,
         true_lower = 32168.22,
         covered = case_when(true_mean >= lower_bound & true_mean <= upper_bound ~ TRUE,
                             TRUE ~ FALSE),
         ci_ratio = (upper_bound - lower_bound) / (true_upper - true_lower))

mean_earned_income_results = mean_earned_income_results %>%
  mutate(dp_mean = (upper_bound - lower_bound) / 2 + lower_bound,
         true_mean = dp_mean - bias,
         true_upper = 15981.18,
         true_lower = 15659.83,
         covered = case_when(true_mean >= lower_bound & true_mean <= upper_bound ~ TRUE,
                             TRUE ~ FALSE),
         ci_ratio = (upper_bound - lower_bound) / (true_upper - true_lower))


## quadrant plot
png('report_figures/cps/mean_income_quadrant_size.png', width = 1600, height = 1800)
ggplot(mean_income_results %>%
         filter(epsilon >= 0.1, 
                epsilon <= 5, 
                #method != 'bhm',
                delta == 0.01 | is.na(delta), 
                ci_overlap != -Inf) %>%
         mutate(epsilon_char = recode(epsilon, 
                                      `0.1` = 'Epsilon = 0.1',
                                      `0.5` = 'Epsilon = 0.5',
                                      `1` = 'Epsilon = 1.0',
                                      `5` = 'Epsilon = 5.0'),
                method = recode(method, 
                                `bhm` = 'BHM'),
                relative_bias = bias / true_mean * 100),
       #group_by(method, epsilon, delta) %>%
       #summarize_at(vars(ci_ratio, ci_overlap), mean, na.rm = TRUE), 
       aes(x = ci_ratio, y = ci_overlap)
           #size = 6)
           #size = abs(relative_bias))
) +
  geom_point(alpha = 0.4,
             size = 6) +
  #geom_density_2d(bins = 5) +
  #xlim(c(0, 2.5)) +
  #ylim(c(0, 1)) +
  geom_vline(xintercept = 1, color = 'red') +
  geom_hline(yintercept = 0.5, color = 'red') +
  facet_wrap(~epsilon_char + method, ncol = 3) +
  ylab('Confidence Interval Overlap') +
  #ylab('Relative Bias (%)') +
  xlab('Confidence Interval Ratio') +
  #scale_size_binned(name = 'Relative Bias (%)', range = c(0, 12), n.breaks = 5) +
  #scale_size(name = 'CIO', breaks = seq(0, 1, 0.1)) +
  theme_minimal(base_size = 24) +
  theme(legend.position = 'top')
#theme_bw(base_size = 12)
dev.off()


png('report_figures/cps/mean_earned_income_quadrant_size.png', width = 1600, height = 1800)
ggplot(mean_earned_income_results %>%
         filter(epsilon >= 0.1, 
                epsilon <= 5, 
                #method != 'bhm',
                delta == 0.01 | is.na(delta), 
                ci_overlap != -Inf,
                !method %in% c('CENQ', 'MOD')) %>%
         mutate(epsilon_char = recode(epsilon, 
                                      `0.1` = 'Epsilon = 0.1',
                                      `0.5` = 'Epsilon = 0.5',
                                      `1` = 'Epsilon = 1.0',
                                      `5` = 'Epsilon = 5.0'),
                method = recode(method, 
                                `bhm` = 'BHM'),
                relative_bias = bias / true_mean * 100),
       #group_by(method, epsilon, delta) %>%
       #summarize_at(vars(ci_ratio, ci_overlap), mean, na.rm = TRUE), 
       aes(x = ci_ratio, y = ci_overlap)
       #size = 6)
       #size = abs(relative_bias))
) +
  geom_point(alpha = 0.4,
             size = 6) +
  #geom_density_2d(bins = 5) +
  #xlim(c(0, 2.5)) +
  #ylim(c(0, 1)) +
  geom_vline(xintercept = 1, color = 'red') +
  geom_hline(yintercept = 0.5, color = 'red') +
  facet_wrap(~epsilon_char + method, ncol = 4) +
  ylab('Confidence Interval Overlap') +
  #ylab('Relative Bias (%)') +
  xlab('Confidence Interval Ratio') +
  #scale_size_binned(name = 'Relative Bias (%)', range = c(0, 12), n.breaks = 5) +
  #scale_size(name = 'CIO', breaks = seq(0, 1, 0.1)) +
  theme_minimal(base_size = 24) +
  theme(legend.position = 'top')
#theme_bw(base_size = 12)
dev.off()



## relative bias
png('report_figures/cps/mean_income_relative_bias.png', width = 1600, height = 900)
ggplot(mean_income_results %>% 
         filter(epsilon >= 0.1, 
                epsilon <= 5, 
                delta == 0.01 | is.na(delta)) %>%
         mutate(epsilon_char = recode(epsilon, 
                                      `0.1` = 'Epsilon = 0.1',
                                      `0.5` = 'Epsilon = 0.5',
                                      `1` = 'Epsilon = 1.0',
                                      `5` = 'Epsilon = 5.0')), 
       aes(x = bias / true_mean, y = method)) +
  geom_violin(draw_quantiles = c(0.5), fill = 'grey50') +
  facet_wrap(~epsilon_char, dir = 'v', nrow = 2, scales = 'free_x') +
  #scale_fill_manual(name = 'delta', 
  #                  values = c('#e41a1c', '#377eb8', '#4daf4a'),
  #                  na.value = 'grey50',
  #                  labels = c('10e-10', '10e-07', '10e-03')) +
  theme_minimal(base_size = 24) +
  scale_y_discrete(labels = c('BHM', 'NOISYMAD', 'NOISYVAR')) +
  scale_x_continuous(labels = scales::percent) +
  #ggtitle('Bias for Mean Income') +
  ylab('') +
  xlab('Relative Bias (%)')
dev.off()

png('report_figures/cps/mean_earned_income_relative_bias.png', width = 1600, height = 900)
ggplot(mean_earned_income_results %>% 
         filter(epsilon >= 0.1, 
                epsilon <= 5, 
                delta == 0.01 | is.na(delta),
                !method %in% c('CENQ', 'MOD')) %>%
         mutate(epsilon_char = recode(epsilon, 
                                      `0.1` = 'Epsilon = 0.1',
                                      `0.5` = 'Epsilon = 0.5',
                                      `1` = 'Epsilon = 1.0',
                                      `5` = 'Epsilon = 5.0')), 
       aes(x = bias / true_mean, y = method)) +
  geom_violin(draw_quantiles = c(0.5), fill = 'grey50') +
  facet_wrap(~epsilon_char, dir = 'v', nrow = 2, scales = 'free_x') +
  #scale_fill_manual(name = 'delta', 
  #                  values = c('#e41a1c', '#377eb8', '#4daf4a'),
  #                  na.value = 'grey50',
  #                  labels = c('10e-10', '10e-07', '10e-03')) +
  theme_minimal(base_size = 24) +
  scale_y_discrete(labels = c('BHM', 'NOISYMAD', 'NOISYVAR', 'SYMQ')) +
  scale_x_continuous(labels = scales::percent) +
  #ggtitle('Bias for Mean Income') +
  ylab('') +
  xlab('Relative Bias (%)')
dev.off()


#####
### Quantiles
#####
quantile_results = read_csv(here::here("results", "03_cps_quantiles.csv"))

ggplot(quantile_results %>% 
         filter(method != 'quantile_smooth', epsilon >= 0.1, epsilon <= 5), 
       aes(x = bias, y = method, fill = as.factor(delta))) +
  geom_violin(draw_quantiles = c(0.5)) +
  facet_wrap(~epsilon, dir = 'v', nrow = 2) +
  scale_fill_manual(name = 'delta', 
                    values = c('#e41a1c', '#377eb8', '#4daf4a'), 
                    na.value = 'grey50',
                    labels = c('10e-10', '10e-07', '10e-03')) +
  theme_minimal(base_size = 24) +
  ggtitle('Bias for Quantiles') +
  ylab('')


## each individual quantile
cps_truth = tibble(
  probs = seq(0.1, 0.9, 0.1),
  truth = c(0, 0, 0, 2500, 16000, 27000, 39000, 52000, 78000)
)

cps_quant = read_csv('results/03_cps_quantiles-detailed.csv')

cps_quant = cps_quant %>%
  left_join(cps_truth, by = 'probs') %>%
  mutate(bias = value - truth,
         relative_bias = case_when(truth != 0 ~ bias / truth,
                                   truth == 0 ~ NA_real_),
         delta = case_when(method == 'quantile_joint_exp' ~ NA_real_,
                           TRUE ~ delta))

## absolute error for zero quantiles
png('report_figures/cps/quantile_absolute_bias.png', width = 1600, height = 900)
ggplot(cps_quant %>%
         filter(method != 'quantile_smooth', epsilon >= 0.1, epsilon <= 5, truth == 0) %>%
         group_by(method, epsilon, delta, replicate) %>%
         summarize(bias = mean(bias)) %>%
         mutate(epsilon_char = recode(epsilon, 
                                      `0.1` = 'Epsilon = 0.1',
                                      `0.5` = 'Epsilon = 0.5',
                                      `1` = 'Epsilon = 1.0',
                                      `5` = 'Epsilon = 5.0')), 
       aes(x = bias, y = method, fill = as.factor(delta))) +
  geom_violin(draw_quantiles = c(0.5)) +
  facet_wrap(~epsilon_char, dir = 'v', nrow = 2, scales = 'free_x') +
  scale_fill_manual(name = 'Delta', 
                    values = c('#e41a1c', '#377eb8', '#4daf4a'), 
                    na.value = 'grey50',
                    labels = c('10e-10', '10e-07', '10e-03')) +
  scale_y_discrete(labels = c('AppIndExp', 'JointExp')) +
  theme_minimal(base_size = 24) +
  #ggtitle('Bias for Quantiles') +
  ylab('') +
  xlab('Absolute Bias')
dev.off()


## relative error for non-zero quantiles
png('report_figures/cps/quantile_relative_bias.png', width = 1600, height = 900)
ggplot(cps_quant %>%
         filter(method != 'quantile_smooth', epsilon >= 0.1, epsilon <= 5, truth != 0) %>%
         group_by(method, epsilon, delta, replicate) %>%
         summarize(relative_bias = median(relative_bias, na.rm = TRUE)) %>%
         mutate(epsilon_char = recode(epsilon, 
                                      `0.1` = 'Epsilon = 0.1',
                                      `0.5` = 'Epsilon = 0.5',
                                      `1` = 'Epsilon = 1.0',
                                      `5` = 'Epsilon = 5.0')), 
       aes(x = relative_bias, y = method, fill = as.factor(delta))) +
  geom_violin(draw_quantiles = c(0.5)) +
  facet_wrap(~epsilon_char, dir = 'v', nrow = 2) +
  scale_fill_manual(name = 'Delta', 
                    values = c('#e41a1c', '#377eb8', '#4daf4a'), 
                    na.value = 'grey50',
                    labels = c('10e-10', '10e-07', '10e-03')) +
  scale_y_discrete(labels = c('AppIndExp', 'JointExp')) +
  scale_x_continuous(labels = scales::percent) +
  theme_minimal(base_size = 24) +
  #ggtitle('Bias for Quantiles') +
  ylab('') +
  xlab('Relative Bias (%)')
dev.off()





#####
## Regression
#####
true_coefficients = tibble(term = c('(Intercept)', 
                                    'non_white1', 
                                    'potential_experience', 
                                    'potential_experience_squared',
                                    'potential_experience_cubed',
                                    'years_of_educ'),
                           true_estimate = c(6.186, 
                                             -1.530e-02, 
                                             1.540e-01, 
                                             -5.422e-03, 
                                             5.814e-05, 
                                             1.613e-01),
                           true_se = c(2.052e-02, 
                                       9.007e-03, 
                                       1.929e-03, 
                                       1.027e-04, 
                                       1.511e-06,
                                       1.393e-03)) %>%
  mutate(ci_lower = true_estimate - 1.96 * true_se,
         ci_upper = true_estimate + 1.96 * true_se,
         ci_length = ci_upper - ci_lower) %>%
  select(term, ci_length)

coefficients_asymptotic <- read_csv(here::here("results", "04_cps_female-regression_coefficient_asymptotic.csv"),
                                    col_types = cols(
                                      method = col_character(),
                                      epsilon = col_double(),
                                      replicate = col_double(),
                                      term = col_character(),
                                      estimate_results = col_double(),
                                      estimate_noisy_results = col_double(),
                                      bias = col_double(),
                                      ci_overlap = col_double(),
                                      ci_length = col_double(),
                                      sign_match = col_logical(),
                                      ci_contains_zero = col_logical(),
                                      delta = col_double(),
                                      n = col_double()
                                    )) %>%
  mutate(method = if_else(method == "Brawner and Hanaker Method", "BHM", method)) %>%
  mutate(method = if_else(!is.na(n), paste0(method, " (n = ", n, ")"), method)) %>%
  left_join(true_coefficients,
            by = "term", suffix = c("_noisy", "_original")) %>%
  mutate(relative_ci_length = ci_length_noisy / ci_length_original,
         relative_bias = bias / estimate_results * 100,
         term = recode(term,
                       `(Intercept)` = 'Intercept',
                       `non_white1` = 'Non-White',
                       `potential_experience` = 'Potential Experience',
                       `potential_experience_squared` = 'Potential Experience Squared',
                       `potential_experience_cubed` = 'Potential Experience Cubed',
                       `years_of_educ` = 'Years of Education'),
         method = recode(method,
                         `WKL Mechanism (ED)` = 'Regularized Normal Mechanism'))


coefficients_bootstrap <- read_csv(here::here("results", "04_cps_female-regression_coefficient_bootstrap.csv"),
                                   col_types = cols(
                                     method = col_character(),
                                     epsilon = col_double(),
                                     replicate = col_double(),
                                     term = col_character(),
                                     estimate_results = col_double(),
                                     estimate_noisy_results = col_double(),
                                     bias = col_double(),
                                     ci_overlap = col_double(),
                                     ci_length = col_double(),
                                     sign_match = col_logical(),
                                     ci_contains_zero = col_logical(),
                                     delta = col_double(),
                                     n = col_double()
                                   )) %>%
  mutate(method = if_else(method == "Brawner and Hanaker Method", "BHM", method)) %>%
  mutate(method = if_else(!is.na(n), paste0(method, " (n = ", n, ")"), method)) %>%
  left_join(true_coefficients,
            by = "term", suffix = c("_noisy", "_original")) %>%
  mutate(relative_ci_length = ci_length_noisy / ci_length_original,
         relative_bias = bias / estimate_results * 100,
         term = recode(term,
                       `(Intercept)` = 'Intercept',
                       `non_white1` = 'Non-White',
                       `potential_experience` = 'Potential Experience',
                       `potential_experience_squared` = 'Potential Experience Squared',
                       `potential_experience_cubed` = 'Potential Experience Cubed',
                       `years_of_educ` = 'Years of Education'),
         method = recode(method,
                         `WKL Mechanism (ED)` = 'Regularized Normal Mechanism'))



## bias plot
png('report_figures/cps/female_regression_absolute_bias_all.png', width = 1200, height = 1600)
coefficients_bootstrap %>% 
  filter(epsilon == 5,
         delta %in% c(NA, 1e-06),
         term != 'Intercept',
         !method %in% c("BHM (n = 10)", "BHM (n = 25)", "Wishart Mechanism", "WKL Mechanism (Pure)")) %>%
  mutate(epsilon_char = recode(epsilon, 
                               `5` = 'Epsilon = 5.0'),
         term = factor(term,
                       levels = c('Non-White', 
                                  'Potential Experience', 
                                  'Potential Experience Squared',
                                  'Potential Experience Cubed',
                                  'Years of Education')),
         relative_bias = bias / estimate_results) %>%
  group_by(method, epsilon, delta, term) %>%
  ungroup() %>%
  ggplot(aes(x = estimate_noisy_results, y = method)) +
  geom_vline(aes(xintercept = estimate_results), col = 'red') +
  geom_density_ridges(alpha = 0.1, fill = "#00AFBB", quantile_lines = TRUE, quantiles = c(0.1, 0.9)) +
  facet_wrap(~epsilon_char + term, ncol = 1, dir = 'v', scales = 'free') +
  theme_minimal(base_size = 24) +
  ylab('') +
  xlab('DP Estimates')
dev.off()

## relative bias plot
png('report_figures/cps/female_regression_relative_bias_all.png', width = 1200, height = 1600)
coefficients_bootstrap %>% 
  filter(epsilon == 5,
         delta %in% c(NA, 1e-06),
         term != 'Intercept',
         !method %in% c("BHM (n = 10)", "BHM (n = 25)", "Wishart Mechanism", "WKL Mechanism (Pure)")) %>%
  mutate(epsilon_char = recode(epsilon, 
                               `5` = 'Epsilon = 5.0'),
         term = factor(term,
                       levels = c('Non-White', 
                                  'Potential Experience', 
                                  'Potential Experience Squared',
                                  'Potential Experience Cubed',
                                  'Years of Education')),
         relative_bias = bias / estimate_results) %>%
  group_by(method, epsilon, delta, term) %>%
  #filter(relative_bias >= quantile(relative_bias, 0.1) & relative_bias <= quantile(relative_bias, 0.9)) %>%
  ungroup() %>%
  ggplot(aes(x = relative_bias, y = method)) +
  geom_vline(aes(xintercept = 0), col = 'red') +
  geom_density_ridges(alpha = 0.1, fill = "#00AFBB", quantile_lines = TRUE, quantiles = c(0.1, 0.9)) +
  facet_wrap(~epsilon_char + term, ncol = 1, dir = 'v', scales = 'free') +
  theme_minimal(base_size = 24) +
  scale_x_continuous(labels = scales::percent) +
  ylab('') +
  xlab('Relative Bias (%)')
dev.off()

## only one method 
png('report_figures/cps/female_regression_absolute_bias_AGM.png', width = 1600, height = 900)
coefficients_bootstrap %>% 
  filter(epsilon <= 20,
         delta %in% c(NA, 1e-06),
         term != 'Intercept',
         method %in% c("Analytic Gaussian Mechanism")) %>%
  mutate(epsilon_char = factor(recode(epsilon, 
                                      `0.5` = 'Epsilon = 0.5',
                                      `1` = 'Epsilon = 1.0',
                                      `5` = 'Epsilon = 5.0',
                                      `10` = 'Epsilon = 10.0',
                                      `15` = 'Epsilon = 15.0',
                                      `20` = 'Epsilon = 20.0'),
                               levels = c('Epsilon = 0.5', 'Epsilon = 1.0', 'Epsilon = 5.0', 'Epsilon = 10.0', 'Epsilon = 15.0', 'Epsilon = 20.0')),
         relative_bias = bias / estimate_results,
         term = factor(term,
                       levels = c('Non-White', 
                                  'Potential Experience', 
                                  'Potential Experience Squared',
                                  'Potential Experience Cubed',
                                  'Years of Education')),
         method = recode(method, `Analytic Gaussian Mechanism` = 'AGM')) %>%
  group_by(method, epsilon, delta, term) %>%
  #filter(relative_bias >= quantile(relative_bias, 0.1) & relative_bias <= quantile(relative_bias, 0.9)) %>%
  ungroup() %>%
  ggplot(aes(x = estimate_noisy_results, y = method)) +
  geom_vline(aes(xintercept = estimate_results), col = 'red') +
  geom_density_ridges(alpha = 0.1, fill = "#00AFBB", quantile_lines = TRUE, quantiles = c(0.1, 0.9)) +
  facet_wrap(~epsilon_char + term, nrow = 5, dir = 'v', scales = 'free') +
  theme_minimal(base_size = 20) +
  ylab('') +
  xlab('DP Estimates')
dev.off()

## only one method - relative bias
png('report_figures/cps/female_regression_relative_bias_AGM.png', width = 1600, height = 900)
coefficients_bootstrap %>% 
  filter(epsilon <= 20,
         delta %in% c(NA, 1e-06),
         term != 'Intercept',
         method %in% c("Analytic Gaussian Mechanism")) %>%
  mutate(epsilon_char = factor(recode(epsilon, 
                                      `0.5` = 'Epsilon = 0.5',
                                      `1` = 'Epsilon = 1.0',
                                      `5` = 'Epsilon = 5.0',
                                      `10` = 'Epsilon = 10.0',
                                      `15` = 'Epsilon = 15.0',
                                      `20` = 'Epsilon = 20.0'),
                               levels = c('Epsilon = 0.5', 'Epsilon = 1.0', 'Epsilon = 5.0', 'Epsilon = 10.0', 'Epsilon = 15.0', 'Epsilon = 20.0')),
         relative_bias = bias / estimate_results,
         term = factor(term,
                       levels = c('Non-White', 
                                  'Potential Experience', 
                                  'Potential Experience Squared',
                                  'Potential Experience Cubed',
                                  'Years of Education')),
         method = recode(method, `Analytic Gaussian Mechanism` = 'AGM')) %>%
  group_by(method, epsilon, delta, term) %>%
  #filter(relative_bias >= quantile(relative_bias, 0.1) & relative_bias <= quantile(relative_bias, 0.9)) %>%
  ungroup() %>%
  ggplot(aes(x = relative_bias, y = method)) +
  geom_vline(aes(xintercept = 0), col = 'red') +
  geom_density_ridges(alpha = 0.1, fill = "#00AFBB", quantile_lines = TRUE, quantiles = c(0.1, 0.9)) +
  scale_x_continuous(labels = scales::percent) +
  facet_wrap(~epsilon_char + term, nrow = 5, dir = 'v', scales = 'free') +
  theme_minimal(base_size = 20) +
  ylab('') +
  xlab('Relative Bias (%)')
dev.off()


## cofidence interval quadrant plots
png('report_figures/cps/female_regression_quadrant_bootstrap_all.png', width = 1600, height = 900)
coefficients_bootstrap %>%
  filter(epsilon == 5,
         delta %in% c(NA, 1e-06),
         term != 'Intercept',
         #term == 'mcg_first_dollar',
         !method %in% c("BHM (n = 10)", "BHM (n = 25)", "Wishart Mechanism", "WKL Mechanism (Pure)")) %>%
  mutate(epsilon_char = recode(epsilon, 
                               `5` = 'Epsilon = 5.0'),
         term = factor(term,
                       levels = c('Non-White', 
                                  'Potential Experience', 
                                  'Potential Experience Squared',
                                  'Potential Experience Cubed',
                                  'Years of Education'))) %>%
  group_by(method, epsilon, delta, term) %>%
  filter(ci_overlap >= quantile(ci_overlap, 0.1) & ci_overlap <= quantile(ci_overlap, 0.9),
         relative_ci_length <= quantile(relative_ci_length, 0.9) & relative_ci_length >= quantile(relative_ci_length, 0.1)) %>%
  ungroup() %>%
  ggplot(aes(x = log(relative_ci_length), y = ci_overlap)) +
  geom_point(alpha = 0.4, size = 4) +
  geom_vline(xintercept = 0, color = 'red') +
  geom_hline(yintercept = 0.5, color = 'red') +
  facet_wrap(~method + term, ncol = 3, dir = 'v') +
  ylab('Confidence Interval Overlap') +
  xlab('Log Confidence Interval Ratio') +
  theme_minimal(base_size = 24) +
  theme(legend.position = 'top')
dev.off()

png('report_figures/cps/female_regression_quadrant_asymptotic_all.png', width = 1600, height = 900)
coefficients_asymptotic %>%
  filter(epsilon == 5,
         delta %in% c(NA, 1e-06),
         term != 'Intercept',
         #term == 'mcg_first_dollar',
         !method %in% c("BHM (n = 10)", "BHM (n = 25)", "Wishart Mechanism", "WKL Mechanism (Pure)")) %>%
  mutate(epsilon_char = recode(epsilon, 
                               `5` = 'Epsilon = 5.0'),
         term = factor(term,
                       levels = c('Non-White', 
                                  'Potential Experience', 
                                  'Potential Experience Squared',
                                  'Potential Experience Cubed',
                                  'Years of Education'))) %>%
  group_by(method, epsilon, delta, term) %>%
  filter(ci_overlap >= quantile(ci_overlap, 0.1) & ci_overlap <= quantile(ci_overlap, 0.9),
         relative_ci_length <= quantile(relative_ci_length, 0.9) & relative_ci_length >= quantile(relative_ci_length, 0.1)) %>%
  ungroup() %>%
  ggplot(aes(x = relative_ci_length, y = ci_overlap)) +
  geom_point(alpha = 0.4, size = 4) +
  geom_vline(xintercept = 1, color = 'red') +
  geom_hline(yintercept = 0.5, color = 'red') +
  facet_wrap(~method + term, ncol = 3, dir = 'v', scales = 'free_y') +
  ylab('Confidence Interval Overlap') +
  xlab('Confidence Interval Ratio') +
  theme_minimal(base_size = 24) +
  theme(legend.position = 'top')
dev.off()


## just Analytic Gaussian Mechanism - or other
png('report_figures/cps/female_regression_quadrant_bootstrap_AGM.png', width = 1600, height = 900)
coefficients_bootstrap %>%
  filter(epsilon <= 20,
         delta %in% c(NA, 1e-06),
         term != 'Intercept',
         #term == 'mcg_first_dollar',
         method %in% c("Analytic Gaussian Mechanism")) %>%
  mutate(epsilon_char = factor(recode(epsilon, 
                                      `0.5` = 'Epsilon = 0.5',
                                      `1` = 'Epsilon = 1.0',
                                      `5` = 'Epsilon = 5.0',
                                      `10` = 'Epsilon = 10.0',
                                      `15` = 'Epsilon = 15.0',
                                      `20` = 'Epsilon = 20.0'),
                               levels = c('Epsilon = 0.5', 'Epsilon = 1.0', 'Epsilon = 5.0', 'Epsilon = 10.0', 'Epsilon = 15.0', 'Epsilon = 20.0')),
         term = factor(term,
                       levels = c('Non-White', 
                                  'Potential Experience', 
                                  'Potential Experience Squared',
                                  'Potential Experience Cubed',
                                  'Years of Education'))) %>%
  group_by(method, epsilon, delta, term) %>%
  filter(ci_overlap >= quantile(ci_overlap, 0.1) & ci_overlap <= quantile(ci_overlap, 0.9),
         relative_ci_length <= quantile(relative_ci_length, 0.9) & relative_ci_length >= quantile(relative_ci_length, 0.1)) %>%
  ungroup() %>%
  ggplot(aes(x = log(relative_ci_length), y = ci_overlap)) +
  geom_point(alpha = 0.4, size = 4) +
  geom_vline(xintercept = 0, color = 'red') +
  geom_hline(yintercept = 0.5, color = 'red') +
  facet_wrap(~epsilon_char + term, nrow = 5, dir = 'v') +
  ylab('Confidence Interval Overlap') +
  xlab('Log Confidence Interval Ratio') +
  theme_minimal(base_size = 24) +
  theme(legend.position = 'top')
dev.off()


png('report_figures/cps/female_regression_quadrant_asymptotic_AGM.png', width = 1600, height = 900)
coefficients_asymptotic %>%
  filter(epsilon <= 20,
         delta %in% c(NA, 1e-06),
         term != 'Intercept',
         #term == 'mcg_first_dollar',
         method %in% c("Analytic Gaussian Mechanism")) %>%
  mutate(epsilon_char = factor(recode(epsilon, 
                                      `0.5` = 'Epsilon = 0.5',
                                      `1` = 'Epsilon = 1.0',
                                      `5` = 'Epsilon = 5.0',
                                      `10` = 'Epsilon = 10.0',
                                      `15` = 'Epsilon = 15.0',
                                      `20` = 'Epsilon = 20.0'),
                               levels = c('Epsilon = 0.5', 'Epsilon = 1.0', 'Epsilon = 5.0', 'Epsilon = 10.0', 'Epsilon = 15.0', 'Epsilon = 20.0')),
         term = factor(term,
                       levels = c('Non-White', 
                                  'Potential Experience', 
                                  'Potential Experience Squared',
                                  'Potential Experience Cubed',
                                  'Years of Education'))) %>%
  group_by(method, epsilon, delta, term) %>%
  filter(ci_overlap >= quantile(ci_overlap, 0.1) & ci_overlap <= quantile(ci_overlap, 0.9),
         relative_ci_length <= quantile(relative_ci_length, 0.9) & relative_ci_length >= quantile(relative_ci_length, 0.1)) %>%
  ungroup() %>%
  ggplot(aes(x = relative_ci_length, y = ci_overlap)) +
  geom_point(alpha = 0.4, size = 4) +
  geom_vline(xintercept = 1, color = 'red') +
  geom_hline(yintercept = 0.5, color = 'red') +
  facet_wrap(~epsilon_char + term, nrow = 5, dir = 'v', scales = 'free_y') +
  ylab('Confidence Interval Overlap') +
  xlab('Confidence Interval Ratio') +
  theme_minimal(base_size = 24) +
  theme(legend.position = 'top')
dev.off()


###
###
## sign and significance match
png('report_figures/cps/female_regression_sign_signif_bootstrap_all.png', width = 1600, height = 900)
coefficients_bootstrap %>%
  filter(epsilon == 5,
         delta %in% c(NA, 1e-06),
         term != 'Intercept',
         #term == 'mcg_first_dollar',
         !method %in% c("BHM (n = 10)", "BHM (n = 25)", "Wishart Mechanism", "WKL Mechanism (Pure)")) %>%
  mutate(epsilon_char = recode(epsilon, 
                               `5` = 'Epsilon = 5.0'),
         term = factor(term,
                       levels = c('Non-White', 
                                  'Potential Experience', 
                                  'Potential Experience Squared',
                                  'Potential Experience Cubed',
                                  'Years of Education')),
         sign_match = factor(sign_match, levels = c('TRUE', 'FALSE'))) %>%
  ggplot(aes(x = sign_match, y = !ci_contains_zero)) +
  geom_confmat(text.colour = 'black', text.size = 8) +
  scale_fill_distiller(name = 'Percent of Results', type = 'seq', palette = 'Reds', direction = 1) +
  facet_wrap(~method + term, dir = 'v', nrow = 5) +
  ylab('Significance Match') +
  xlab('Sign Match') +
  theme_minimal(base_size = 24) +
  theme(legend.position = 'none')
dev.off()

png('report_figures/cps/female_regression_sign_signif_asymptotic_all.png', width = 1600, height = 900)
coefficients_asymptotic %>%
  filter(epsilon == 5,
         delta %in% c(NA, 1e-06),
         term != 'Intercept',
         #term == 'mcg_first_dollar',
         !method %in% c("BHM (n = 10)", "BHM (n = 25)", "Wishart Mechanism", "WKL Mechanism (Pure)")) %>%
  mutate(epsilon_char = recode(epsilon, 
                               `5` = 'Epsilon = 5.0'),
         term = factor(term,
                       levels = c('Non-White', 
                                  'Potential Experience', 
                                  'Potential Experience Squared',
                                  'Potential Experience Cubed',
                                  'Years of Education')),
         sign_match = factor(sign_match, levels = c('TRUE', 'FALSE'))) %>%
  ggplot(aes(x = sign_match, y = !ci_contains_zero)) +
  geom_confmat(text.colour = 'black', text.size = 8) +
  scale_fill_distiller(name = 'Percent of Results', type = 'seq', palette = 'Reds', direction = 1) +
  facet_wrap(~method + term, dir = 'v', nrow = 5) +
  ylab('Significance Match') +
  xlab('Sign Match') +
  theme_minimal(base_size = 24) +
  theme(legend.position = 'none')
dev.off()

## just WKL (ED) - or other methods
png('report_figures/cps/female_regression_sign_signif_bootstrap_AGM.png', width = 1600, height = 900)
coefficients_bootstrap %>%
  filter(epsilon <= 20,
         delta %in% c(NA, 1e-06),
         term != 'Intercept',
         #term == 'mcg_first_dollar',
         method %in% c("Analytic Gaussian Mechanism")) %>%
  mutate(epsilon_char = factor(recode(epsilon, 
                                      `0.5` = 'Epsilon = 0.5',
                                      `1` = 'Epsilon = 1.0',
                                      `5` = 'Epsilon = 5.0',
                                      `10` = 'Epsilon = 10.0',
                                      `15` = 'Epsilon = 15.0',
                                      `20` = 'Epsilon = 20.0'),
                               levels = c('Epsilon = 0.5', 'Epsilon = 1.0', 'Epsilon = 5.0', 'Epsilon = 10.0', 'Epsilon = 15.0', 'Epsilon = 20.0')),
         term = factor(term,
                       levels = c('Non-White', 
                                  'Potential Experience', 
                                  'Potential Experience Squared',
                                  'Potential Experience Cubed',
                                  'Years of Education')),
         sign_match = factor(sign_match, levels = c('TRUE', 'FALSE'))) %>%
  ggplot(aes(x = sign_match, y = !ci_contains_zero)) +
  geom_confmat(text.colour = 'black', text.size = 8) +
  scale_fill_distiller(name = 'Percent of Results', type = 'seq', palette = 'Reds', direction = 1) +
  facet_wrap(~epsilon_char + term, dir = 'v', nrow = 5) +
  ylab('Significance Match') +
  xlab('Sign Match') +
  theme_minimal(base_size = 24) +
  theme(legend.position = 'none')
dev.off()

png('report_figures/cps/female_regression_sign_signif_asymptotic_AGM.png', width = 1600, height = 900)
coefficients_asymptotic %>%
  filter(epsilon <= 20,
         delta %in% c(NA, 1e-06),
         term != 'Intercept',
         #term == 'mcg_first_dollar',
         method %in% c("Analytic Gaussian Mechanism")) %>%
  mutate(epsilon_char = factor(recode(epsilon, 
                                      `0.5` = 'Epsilon = 0.5',
                                      `1` = 'Epsilon = 1.0',
                                      `5` = 'Epsilon = 5.0',
                                      `10` = 'Epsilon = 10.0',
                                      `15` = 'Epsilon = 15.0',
                                      `20` = 'Epsilon = 20.0'),
                               levels = c('Epsilon = 0.5', 'Epsilon = 1.0', 'Epsilon = 5.0', 'Epsilon = 10.0', 'Epsilon = 15.0', 'Epsilon = 20.0')),
         term = factor(term,
                       levels = c('Non-White', 
                                  'Potential Experience', 
                                  'Potential Experience Squared',
                                  'Potential Experience Cubed',
                                  'Years of Education')),
         sign_match = factor(sign_match, levels = c('TRUE', 'FALSE'))) %>%
  ggplot(aes(x = sign_match, y = !ci_contains_zero)) +
  geom_confmat(text.colour = 'black', text.size = 8) +
  scale_fill_distiller(name = 'Percent of Results', type = 'seq', palette = 'Reds', direction = 1) +
  facet_wrap(~epsilon_char + term, dir = 'v', nrow = 5) +
  ylab('Significance Match') +
  xlab('Sign Match') +
  theme_minimal(base_size = 24) +
  theme(legend.position = 'none')
dev.off()




