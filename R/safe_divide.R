#################################################################################
#	Safe Divide
#################################################################################

# Safely divide values
# INPUTS:
#     num: numerator (num)
#   denom: denominator (num)
# OUTPUTS:
#     0 or num / denom (num)

safe_divide <- function(num, denom) {
  if (num == 0) {
    return(0)
  } else {
    return(num / denom)
  }
}