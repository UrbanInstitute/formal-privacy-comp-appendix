#' Summarize formally private linear models
#'
#' @param results Results from tidy.lm from library(broom)
#' @param noisy_results Results from a formally private regression
#' with output formatted like tidy.lm from library(broom)
#'
#' @return A tibble with coefficient estimates bias, confidence interval 
#' overlap, confidence interval length, sign match, and if the confidence 
#' interval contains zero.
#' 
summarize_lm <- function(results, noisy_results) {
  
  combined_data <- dplyr::left_join(
    x = results, 
    y = noisy_results, 
    suffix = c("_results", "_noisy_results"),
    by = "term"
  )
  
  results <- combined_data %>%
    dplyr::mutate(
      bias = estimate_noisy_results - estimate_results,
      ci_overlap = compute_confidence_interval_overlap(
        private_lower = conf.low_results,
        private_upper = conf.high_results,
        nonprivate_lower = conf.low_noisy_results,
        nonprivate_upper = conf.high_noisy_results
      ),
      ci_length = conf.high_noisy_results - conf.low_noisy_results,
      sign_match = estimate_results * estimate_noisy_results > 0,
      ci_contains_zero = conf.low_noisy_results < 0 & conf.high_noisy_results > 0
    ) %>%
    dplyr::select(
      term, 
      estimate_results, 
      estimate_noisy_results,
      bias, 
      ci_overlap, 
      ci_length, 
      sign_match,
      ci_contains_zero
    )
  
  return(results)
  
}
