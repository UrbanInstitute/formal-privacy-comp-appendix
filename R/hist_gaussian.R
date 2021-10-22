#' Construct a histogram with a Gaussian mechanism
#'
#' @param x A numeric vector
#' @param breaks A numeric vector of breaks for the histogram 
#' @param epsilon A numeric value for epsilon
#' @param delta A numeric value for delta
#'
#' @return A list with the inputted epsilon (numeric), delta (numeric), 
#' breaks (numeric vector), and counts (numeric vector) of the noisy counts 
#' within each break.
#' 
hist_gaussian <- function(x, breaks, epsilon, delta) {

  # true histogram
  true_histogram <- hist(x, breaks = breaks, plot = FALSE)
  
  # generate object with counts
  tab <- x %>% 
    cut(breaks = true_histogram$breaks) %>%
    janitor::tabyl()
  
  breaks <- tab$.
  
  counts <- tab$n

  # apply Guassian mechanism
  noisy_counts <- gauss_san(x = counts, eps = epsilon, delta = delta)
  
  output <- list(
    epsilon = epsilon,
    delta = delta,
    breaks = breaks,
    counts = noisy_counts
  )
  
  return(output)
  
}
