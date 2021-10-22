#################################################################################
#	DP Count Statistical Inference
#################################################################################

# Function to summarize results from replicates
# INPUTS:
#     x: true value (num)
#     y: estimate (num)
#     n: number of records (int)
#     z: confidence interval overlap value
# OUTPUTS:
#   estimate: estimate from the replicates (num)
#   variance: variance of estimate (num)
#      ci_up: confidence interval upper bound (num)
#      ci_lo: confidence interval lower bound (num)
#         ci: if the confidence interval captured x (int)
#   ci_width: confidence interval width (num)

count_inference <- function(x, y, n, z) {
  
  alpha <- (1 - z) / 2
  
  est <- y
  p <- est / n
  variance <- n * p * (1 - p)
  ci_lo <- (p + sqrt(p * (1 - p) / n) * qnorm(alpha)) * n
  ci_up <- (p + sqrt(p * (1 - p) / n) * qnorm(1 - alpha)) * n
  ci_width <- ci_up - ci_lo
  ci <- ifelse(x <= ci_up & x >= ci_lo, 1, 0)
  
  return(list(
    "estimate" = est,
    "variance" = variance,
    "ci_up" = ci_up,
    "ci_lo" = ci_lo,
    "ci" = ci,
    "ci_width" = ci_width
  ))
}
