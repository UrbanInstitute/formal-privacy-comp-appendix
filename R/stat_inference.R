#################################################################################
#	Statistical Inference Function
#################################################################################

# Function to summarize results from replicates
# INPUTS:
#     x: true value (num)
#     y: vector of values (num)
#     N: number of replicates (int)
# OUTPUTS:
#   estimate: estimate from the replicates (num)
#       bias: bias (num)
#       rmse: Root Mean Square Error (num)

stat_inference <- function(x, y, N) {
  est <- mean(y)
  bias <- mean(est - x)
  rmse <- sqrt(sum((y - est))^2 / (N))
  
  return(list(
    "estimate" = est,
    "bias" = bias,
    "rmse" = rmse
  ))
}