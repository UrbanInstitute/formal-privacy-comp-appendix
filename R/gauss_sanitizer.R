#################################################################################
# 	Gaussian Sanitizer
#################################################################################
# Function implementing Gaussian mechanism instead of Laplace
#
# INPUTS:
#        x: vector of counts (num)
#      eps: epsilon, privacy loss (num)
#    delta: delta tuning parameter that is a non-zero probability (num)
# OUTPUTS:
#    x_san: vector of sanitized counts by the Gaussian mechanism

gauss_san <- function(x, eps, delta) {
  # Checking whether x is a vector
  if (is.null(dim(x)) == FALSE) {
    stop("The x must be vector.")
  }
  
  # Releasing sanitized queries using the Gaussian mechanism.
  x_san <- sapply(x, function(x) max(x + gauss_mech(eps, delta, 1), 0)) %>%
    `/`(sum(.)) %>%
    `*`(sum(x)) %>%
    smart_round(.)
  
  return(x_san)
}