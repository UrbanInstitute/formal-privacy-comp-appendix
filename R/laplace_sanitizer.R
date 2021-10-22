#################################################################################
# 	Laplace sanitizer
#################################################################################
# Function implementing the Laplace sanitizer from Abowd and Vilhuber (2008)
#
# INPUTS:
#        x: vector of counts (num)
#      eps: epsilon, privacy loss (num)
#       gs: the global sensitivity, default is 1 for counting queries (num)
#       lb: the lower bound, default is 0 for counting queries (num)
# OUTPUTS:
#    x_san: vector of sanitized counts by the Laplace mechanism

lap_san <- function(x, eps, gs = 1, lb = 0) {
  
  # Checking whether x is a vector
  if (is.null(dim(x)) == FALSE) {
    stop("The x must be vector.")
  }
  
  # Releasing sanitized queries using the Laplace mechanism
  x_san <- sapply(x, function(x) max(x + lap_mech(eps, gs), lb)) %>%
    `/`(sum(.)) %>%
    `*`(sum(x)) %>%
    smart_round(.)
  
  return(x_san)
}