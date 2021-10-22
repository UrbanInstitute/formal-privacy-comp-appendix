#################################################################################
#	Smart Round
#################################################################################
#
# Smartly round numbers to integers. This functions works
# accurately on grouped dataframes (ie applies the rounding within each group)
#
# INPUTS:
#     x: vector of values (num)
# OUTPUTS:
#     y: a vector of integers (int)

smart_round <- function(x) {
  y <- floor(x)
  indices <- tail(order(x - y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  
  return(y)
}