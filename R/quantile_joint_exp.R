#' Calculate formally private quantiles with joint_exp
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
quantile_joint_exp <- function(x, probs, epsilon, delta, data_low, data_high) {
  
  # test if x is sorted
  if(!(all(diff(x) >= 0))) {
    
    stop("x vector decreases!")
    
  }
  
  # python setup
  reticulate::use_condaenv(condaenv = "r-reticulate", conda = "auto", required = TRUE)
  
  sys <- reticulate::import("sys", convert = FALSE)
  sys$path$insert(0L, 'dp-code/quantile-functions/')
  
  joint_exp <- reticulate::source_python(file = here::here("dp-code", "quantile-functions", "joint_exp.py"))
  gcpy <- reticulate::import(module = "gc")

  x <- reticulate::r_to_py(x) # needed on older versions of reticulate
  
  quantiles <- py$joint_exp(
    sorted_data = x, 
    data_low = data_low, 
    data_high = data_high, 
    qs = probs, 
    eps = epsilon, 
    swap = TRUE
  )

  gcpy$collect()
  rm(joint_exp, gcpy)
  gc()

  return(quantiles)
}
