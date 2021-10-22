#################################################################################
#	Gaussian Mechanism - Multiple "single" Draws
#################################################################################

# Function for adding approximate or (epsilon, delta)-DP noise to a statistic or
# query.
#
# INPUTS:
#   n: number of draws (int)
#   eps: epsilon, privacy loss (num)
#   delta: delta tuning parameter that is a non-zero probability (num)
#   gs: L2 global sensitivity of the specific query (num)
# OUTPUTS:
#   x: n number of draws from a normal distribution that satisifies approx-DP (num)
gauss_mech_multi <- function (n, eps, delta, gs) {
  
  # Checking for proper values
  if (any(eps <= 0)) {
    stop("The eps must be positive.")
  }
  if (any(delta < 0 || delta > 1)) {
    stop("The delta must be between 0 and 1")
  }
  if (any(gs <= 0)) {
    stop("The gs must be positive.")
  }
  
  # Standard Deviation for the Gaussian Mechanism
  gauss_sd <- (gs / eps) * sqrt(2 * log(1.25 / delta))
  
  # Draw from a normal distribution
  x <- rnorm(n, mean = 0, sd = gauss_sd)
  
  return(x)
}