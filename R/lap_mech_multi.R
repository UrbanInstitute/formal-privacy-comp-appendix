#################################################################################
#	Laplace Mechanism - Multiple Draws
#################################################################################

# Function for drawing random Laplace distribution values, code originally from
# Laplace's Demon, but now in the R pkg rmutil
#
# INPUTS:
#          n: number of draws (int)
#        eps: epsilon, privacy loss (num)
#         gs: L1 global sensitivity of the specific query (num)
# OUTPUTS:
#         x: n number of draws from the lap. dist. that satisifies approx-DP (num)
lap_mech_multi <- function(n, eps, gs) {
  
  # Checking for proper values
  if (any(eps <= 0)) {
    stop("The eps must be positive.")
  }
  if (any(n <= 0)) {
    stop("The n must be positive.")
  }
  if (any(gs <= 0)) {
    stop("The gs must be positive.")
  }
  
  # Calculating the scale
  scale <- gs / eps

  r <- r2 <- runif(n)
  temp <- which(r > 0.5)
  r2[temp] <- 1 - r[temp]
  x <- 0 - sign(r - 0.5) * scale * log(2 * r2)

  return(x)
}