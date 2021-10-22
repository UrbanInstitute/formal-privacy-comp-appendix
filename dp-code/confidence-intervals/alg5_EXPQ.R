######
# This file includes the code to compute exponential quantiles as described in Algorithm
# 5, EXPQ. To reduce computation time, the main steps are done in c language. (See file
# "exp_mech.c")
######

#install.packages("rmutil")
#install.packages("here")
library(rmutil)
library(here)

system(paste0("R CMD SHLIB ", here("dp-code", "confidence-intervals", "exp_mech.c")))
dyn.load(here("dp-code", "confidence-intervals", "exp_mech.so"))


############
# The function priv_median_c() takes up to 6 parameters and returns a differentially 
# private quantile estimate: 
#   db:         Data as a vector of numbers.
#   q:          The quantile needed; q = 0.5 for median.
#   e:          Privacy parameter. 
#   xmin:       Given lower bound of the data range.
#   xmax:       Given upper bound of the data range.
#   is_sorted:  TRUE if the given data has been sorted; FALSE if not.
############
priv_median_c <- function(db, q, e, xmin, xmax, is_sorted = FALSE) {
  n <- length(db)
  qi <- floor((n - 1) * q + 1.5)
  
  db[db > xmax] <- xmax
  db[db < xmin] <- xmin
  db <- c(xmin, db, xmax)
  if(!is_sorted) {
    db <- sort(db)
  }
  
  result <- .C("priv_median_exp", db = as.double(db),
               n = as.integer(n), e = as.double(e), qi = as.integer(qi), 
               probs = as.double(1:(n + 1)), r = as.double(runif(1)), priv_qi = as.integer(n + 1))
  
  priv_qi <- result$priv_qi
  
  return(runif(1, db[priv_qi], db[priv_qi + 1]))
}