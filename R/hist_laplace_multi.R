#' Construct a histogram with a Laplace sanitizer
#'
#' @param x A numeric vector
#' @param breaks A numeric vector of breaks for the histogram 
#' @param epsilon A numeric value for epsilon
#'
#' @return A list with the inputted epsilon (numeric), 
#' breaks (numeric vector), and counts (numeric vector) of the noisy counts 
#' within each break.
#' 
hist_laplace_multi <- function(x, breaks, epsilon) {
  
  # global sensitivity
  gs <- 1
  
  # true histogram
  true_histogram <- hist(x, breaks = breaks, plot = FALSE)
  
  # generate object with counts
  tab <- x %>% 
    cut(breaks = true_histogram$breaks) %>%
    janitor::tabyl()
  
  breaks <- tab$.
    
  counts <- tab$n
  
  m <- length(counts)
  
  # Apply Laplace mechanism
  noisy_counts <- map_dbl(
    .x = counts, 
    .f = ~max(.x + lap_mech(eps = epsilon / m, gs = gs), 0)
  )
  
  output <- list(
    epsilon = epsilon,
    breaks = breaks,
    counts = noisy_counts
  )
    
  return(output)
    
}
