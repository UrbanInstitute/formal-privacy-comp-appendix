library(tidyverse)
library(ggridges)

source(here::here("R", "prep_cps03.R"))

mortenson_cps = prep_cps03()

# values for specifications
breaks = seq(from = 0, to = 30000, by = 1000)

# no noise ----------------------------------------------------------------
true_histogram = hist(mortenson_cps$earned_income_2014, breaks = breaks, plot = FALSE)$counts
#true_histogram = tibble()

## load in detailed results
noisy_results = read_csv('../../03_cps_histograms_detailed.csv')
noisy_results

max_error_function = function(noisy_histogram, true_histogram){
  max_low = max(abs(cumsum(true_histogram) - cumsum(noisy_histogram))) / sum(true_histogram)
  max_high = max(abs(cumsum(rev(true_histogram)) - cumsum(rev(noisy_histogram)))) / sum(rev(true_histogram))
  
  max_out = max(max_low, max_high)
  return(max_out)
}

mean_cum_error_function = function(noisy_histogram, true_histogram){
  mean_low = mean(cumsum(abs(true_histogram - noisy_histogram)) / cumsum(true_histogram))
  mean_high = mean(cumsum(abs(rev(true_histogram) - rev(noisy_histogram))) / cumsum(rev(true_histogram)))
  
  mean_out = mean(mean_low, mean_high)
  return(mean_out)
}

## iterate over each 
temp_output_tibble = vector('list', 72)
count = 1
for(temp_method in unique(noisy_results$method)){
  for(temp_eps in unique(noisy_results$epsilon)){
    if(temp_method %in% c('Gaussian Mechanism', 'Gaussian Multiple Queries')){
      for(temp_delta in setdiff(unique(noisy_results$delta), NA)){
        temp_test = noisy_results %>%
          filter(method == temp_method, epsilon == temp_eps, delta == temp_delta)
        
        new_out = new_out_alt = new_out_rev = new_out_rev_alt = canberra_out = l1_out = rep(NA, max(noisy_results$replicate))
        for(temp_rep in 1:max(noisy_results$replicate)){
          histogram_results = filter(temp_test, replicate == temp_rep) %>% pull(noisy_count)
          
          #ks_out[temp_rep] = ks.test(true_histogram, histogram_results)$statistic
          
          new_out[temp_rep] = max(abs(cumsum(true_histogram) - cumsum(histogram_results))) / sum(true_histogram)
          new_out_alt[temp_rep] = mean(cumsum(abs(true_histogram - histogram_results)) / cumsum(true_histogram))
          
          new_out_rev[temp_rep] = max(abs(cumsum(rev(true_histogram)) - cumsum(rev(histogram_results)))) / sum(rev(true_histogram))
          new_out_rev_alt[temp_rep] = mean(cumsum(abs(rev(true_histogram) - rev(histogram_results))) / cumsum(rev(true_histogram)))
          
          canberra_out[temp_rep] = dist(rbind(true_histogram, histogram_results), method = 'canberra')
          
          l1_out[temp_rep] = dist(rbind(true_histogram, histogram_results), method = 'manhattan')
        }
        temp_output_tibble[[count]] = tibble(method = rep(temp_method, max(noisy_results$replicate) * 6), 
                                             epsilon = rep(temp_eps, max(noisy_results$replicate) * 6), 
                                             delta = rep(temp_delta, max(noisy_results$replicate) * 6), 
                                             measure = rep(c('max_error', 
                                                             'max_error_rev', 
                                                             'mean_relative_error', 
                                                             'mean_relative_error_rev', 
                                                             'canberra', 
                                                             'l1'), 
                                                           each = max(noisy_results$replicate)),
                                             replicate = rep(1:max(noisy_results$replicate), 6),
                                             value = c(new_out, new_out_rev, new_out_alt, new_out_rev_alt, canberra_out, l1_out))
        count = count + 1
      }
    } else{
      temp_test = noisy_results %>%
        filter(method == temp_method, epsilon == temp_eps)
      
      new_out = new_out_alt = new_out_rev = new_out_rev_alt = canberra_out = l1_out = rep(NA, max(noisy_results$replicate))
      for(temp_rep in 1:max(noisy_results$replicate)){
        histogram_results = filter(temp_test, replicate == temp_rep) %>% pull(noisy_count)
        
        #ks_out[temp_rep] = ks.test(true_histogram, histogram_results)$statistic
        
        new_out[temp_rep] = max(abs(cumsum(true_histogram) - cumsum(histogram_results))) / sum(true_histogram)
        new_out_alt[temp_rep] = mean(cumsum(abs(true_histogram - histogram_results)) / cumsum(true_histogram))
        
        new_out_rev[temp_rep] = max(abs(cumsum(rev(true_histogram)) - cumsum(rev(histogram_results)))) / sum(rev(true_histogram))
        new_out_rev_alt[temp_rep] = mean(cumsum(abs(rev(true_histogram) - rev(histogram_results))) / cumsum(rev(true_histogram)))
        
        canberra_out[temp_rep] = dist(rbind(true_histogram, histogram_results), method = 'canberra')
        
        l1_out[temp_rep] = dist(rbind(true_histogram, histogram_results), method = 'manhattan')
      }
      temp_output_tibble[[count]] = tibble(method = rep(temp_method, max(noisy_results$replicate) * 6), 
                                           epsilon = rep(temp_eps, max(noisy_results$replicate) * 6), 
                                           delta = rep(NA, max(noisy_results$replicate) * 6), 
                                           measure = rep(c('max_error', 
                                                           'max_error_rev', 
                                                           'mean_relative_error', 
                                                           'mean_relative_error_rev', 
                                                           'canberra', 
                                                           'l1'), 
                                                         each = max(noisy_results$replicate)),
                                           replicate = rep(1:max(noisy_results$replicate), 6),
                                           value = c(new_out, new_out_rev, new_out_alt, new_out_rev_alt, canberra_out, l1_out))
      count = count + 1
    }
    cat(temp_eps, '\n')
  }
  cat(temp_method, '\n')
}

output_tibble = bind_rows(temp_output_tibble)
output_tibble

## save output
saveRDS(output_tibble, 'report_figures/cps_histogram_summaries.rds')






