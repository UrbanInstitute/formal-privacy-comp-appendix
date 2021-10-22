#################################################################################
#	Gaussian Mechanism - Renyi Composition
#################################################################################

# Function for adding approximate or (epsilon, delta)-DP noise to a statistic or
# m queries with Renyi Composition where the sensivitiy is 1
#
# INPUTS:
#   m: number of "folds" (int)
#   eps: epsilon, privacy loss (num)
#   delta: delta tuning parameter that is a non-zero probability (num)
# OUTPUTS:
#   x: single draw from a normal distribution that satisifies approx-DP (num)
gauss_renyi <- function (m, eps, delta) {
  
  # Checking for proper values
  if (any(eps <= 0)) {
    stop("The eps must be positive.")
  }
  if (any(delta < 0 || delta > 1)) {
    stop("The delta must be between 0 and 1")
  }
  if (any(m <= 0)) {
    stop("The m must be positive.")
  }
  
  # Calculate the updated epsilon
  p <- log(1 / delta) / eps
  sig <- m / (2 * eps) * (2 * sqrt(3 * p^2 + p) + (2 * p + 1))
  eps0 <- sqrt(1 / sig)
  
  # Standard Deviation for the Gaussian Mechanism under Renyi Composition
  gauss_sd <- 1 / eps0
  
  # Draw from a normal distribution
  x <- rnorm(1, mean = 0, sd = gauss_sd)
  
  return(x)
}