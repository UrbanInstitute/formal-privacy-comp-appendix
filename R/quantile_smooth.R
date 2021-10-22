#' Calculate formally private quantiles with smooth
#'
#' @param x a non-decreasing numeric vector whose sample quantiles are wanted
#' @param probs numeric vector of probabilities with values in [0,1].
#' @param epsilon numeric value for the desired value of epsilon 
#' @param delta numeric value for the desired value of delta 
#' @param data_low numeric value for the minimum value in the data 
#' @param data_high numeric value for the maximum value in the data 
#'
#' @return a numeric vector of formally private quantiles
#' 
quantile_smooth <- function(x, probs, epsilon, delta, data_low, data_high) {
  
  # test if x is sorted
  if(!(all(diff(x) >= 0))) {
    
    stop("x vector decreases!")
    
  }
  
  # python setup
  reticulate::use_condaenv(condaenv = "r-reticulate", conda = "auto", required = TRUE)
  
  sys <- reticulate::import("sys", convert = FALSE)
  sys$path$insert(0L, 'dp-code/quantile-functions/')
  
  smooth <- reticulate::source_python(file = here::here("dp-code", "quantile-functions", "smooth.py"))
  gcpy <- reticulate::import(module = "gc")

  divided_epsilon <- epsilon / length(probs)
  divided_delta <- delta / length(delta)
  
  x <- reticulate::r_to_py(x) # needed on older versions of reticulate
  
  quantiles <- py$smooth(
    sorted_data = x, 
    data_low = data_low, 
    data_high = data_high, 
    qs = probs, 
    divided_eps = divided_epsilon,
    divided_delta = divided_delta
  )

  gcpy$collect()
  rm(smooth, gcpy)
  gc()

  return(quantiles)
}
