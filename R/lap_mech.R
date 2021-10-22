#################################################################################
#	Laplace Mechanism
#################################################################################

# Function for drawing random Laplace distribution values, code originally from
# Laplace's Demon, but now in the R pkg rmutil
#
# INPUTS:
#        eps: epsilon, privacy loss (num)
#         gs: L1 global sensitivity of the specific query (num)
# OUTPUTS:
#         x: a single draw from the lap. dist. that satisifies approx-DP (num)
lap_mech <- function(eps, gs) {
  
  # Checking for proper values
  if (any(eps <= 0)) {
    stop("The eps must be positive.")
  }
  if (any(gs <= 0)) {
    stop("The GS must be positive.")
  }
  
  # Calculating the scale
  scale <- gs / eps

  r <- runif(1)

  if(r > 0.5) {
    r2 <- 1 - r
    x <- 0 - sign(r - 0.5) * scale * log(2 * r2)
  } else {
    x <- 0 - sign(r - 0.5) * scale * log(2 * r)
  }
  
  return(x)
}