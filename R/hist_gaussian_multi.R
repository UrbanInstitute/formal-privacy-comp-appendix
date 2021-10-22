#' Construct a histogram with a Gaussian mechanism for multiple queries
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
hist_gaussian_multi <- function(x, breaks, epsilon, delta) {

  # true histogram
  true_histogram <- hist(x, breaks = breaks, plot = FALSE)
  
  # generate object with counts
  tab <- x %>% 
    cut(breaks = true_histogram$breaks) %>%
    janitor::tabyl()
  
  breaks <- tab$.
  
  counts <- tab$n

  m <- length(counts)
  
  # apply Guassian mechanism
  noisy_counts <- map_dbl(
    .x = counts,
    .f = ~max(.x + gauss_renyi(m = m, eps = epsilon, delta = delta), 0)
  )
  
  output <- list(
    epsilon = epsilon,
    delta = delta,
    breaks = breaks,
    counts = noisy_counts
  )
  
  return(output)
  
}
