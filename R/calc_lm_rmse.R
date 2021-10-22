#' Calculate rmse with noisy coefficients
#'
#' @param data A tibble with the original data
#' @param formula A formula for the estimate model
#' @param truth A numeric vector with the true outcome variable
#' @param noisy_results A noisy table of model estimates
#'
#' @return A numeric with the calculated rmse
#'
calc_lm_rmse <- function(data, formula, truth, noisy_results) {
  
  # create a design matrix
  pred_data <- modelr::model_matrix(formula, data = data) %>%
    as.matrix()
  
  # test names
  model_matrix_names <- colnames(pred_data)
  estimate_names <- noisy_results$term
  stopifnot(all(model_matrix_names == estimate_names))
  
  # calculate rmse
  rmse <- bind_cols(
    truth = truth, 
    estimate = as.vector(pred_data %*% matrix(noisy_results$estimate, ncol = 1))
  ) %>%
    yardstick::rmse(truth = .data$truth, estimate = .data$estimate)
  
  return(rmse$.estimate)
  
}