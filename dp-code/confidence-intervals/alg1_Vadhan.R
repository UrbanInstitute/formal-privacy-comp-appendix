######
# This file includes the code for the algorithm proposed by Karwa and Vadhan, as
# discussed in Section 3.1 in the paper. 
######
library('rmutil')

## not used
#pub_interval <- function(db, a) {
#  m <- mean(db)
#  radius <- sqrt(variance(db, m) / length(db)) * qnorm(1 - a / 2)
#  return(c(m - radius, m + radius))
#}

## not used
#pub_range <- function(db) {
#  return(c(min(db), max(db)))
#}

## fundamental: why is this used instead of which.max??
#maxi <- function(v) {
#  l <- 1
#  for(i in 2:length(v)) {
#    if(v[l] < v[i]) {
#      l <- i
#    }
#  }
#  return(l)
#}

## fundamental
variance <- function(x, m) {
  return((1 / (length(x) - 1)) * sum((x - m) ^ 2))
}

## fundamental
pub_histogram_learner <- function(db, bins) {
  #outputs a normalized histogram of db separated by the intervals in bins
  db <- sort(db)
  probs <- rep(0, length(bins) - 1)
  db_i <- 1
  while(db[db_i] < bins[1]) {
    db_i <- db_i + 1
    if(db_i > length(db)) {
      return(probs / sum(probs))
    }
  }
  for(i in 1:length(probs)) {
    while(db[db_i] < bins[i + 1]) {
      probs[i] <- probs[i] + 1
      db_i <- db_i + 1
      if(db_i > length(db)) {
        return(probs / sum(probs))
      }
    }
  }
  return(probs / sum(probs))
}

## calls pub_histogram_learner
priv_histogram_learner <- function(db, bins, e) {
  probs <- pub_histogram_learner(db, bins)
  return(probs + rlaplace(length(probs), 0, 2 / e / length(db)))
}

## calls priv_histogram_learner
priv_std <- function(db, a, e, stdmin, stdmax) {
  bins_base <- floor(log2(stdmin) - 2)
  bins <- 2 ^ (bins_base:ceiling(log2(stdmax) + 2))
  y <- 1:floor(length(db) / 2)
  for(i in 1:length(y)) {
    y[i] <- abs(db[2 * i] - db[2 * i - 1])
  }
  
  l <- which.max(priv_histogram_learner(y, bins, e))
  return(2 ^ ((l + bins_base - 1) + 2))
}

## calls priv_histogram_learner
priv_mean <- function(db, e, std, r) {
  rn <- ceiling(r / std)
  bins_base <- -rn
  bins <- ((bins_base:(rn + 1)) - .5) * std
  
  l <- which.max(priv_histogram_learner(db, bins, e))
  return((l + bins_base - 1) * std)
}

## calls priv_std, priv_mean
priv_range <- function(db, a1, a2, e1, e2, stdmin, stdmax, r) {
  priv_std <- priv_std(db, a1, e1, stdmin, stdmax)
  priv_mean <- priv_mean(db, e2, priv_std, r)
  
  radius <- 4 * priv_std * sqrt(log(length(db) / a2))
  return(c(priv_mean - radius, priv_mean + radius))
}

## calls priv_range, variance
priv_vadhan <- function(db, a0, a1, a2, a3, e1, e2, e3, stdmin, stdmax, r) { #stdmin, stdmax, r) {
  n <- length(db)
  #priv_range is the most significant source of error
  xrange <- priv_range(db, a3 / 2, a3 / 2, e3 / 2, e3 / 2, stdmin, stdmax, r)
  xmin <- xrange[1]
  xmax <- xrange[2]
  xdist <- xmax - xmin
  
  #clamp
  db[db < xmin] <- xmin
  db[db > xmax] <- xmax
  
  mean_var <- xdist / (e1 * n)
  priv_mean <- mean(db) + rlaplace(1, 0, mean_var)
  if(priv_mean < xmin) {
    priv_mean = xmin
  } else if(priv_mean > xmax) {
    priv_mean = xmax
  }
  
  var_var <- xdist ^ 2 / (e2 * (n - 1))
  #priv_var <- public_var + extra_var + lap_noise
  priv_var <- variance(db, priv_mean) + var_var * log(1 / a2) + rlaplace(1, 0, var_var)
  if(priv_var < 0 || priv_var > stdmax ^ 2) {
    priv_var = stdmax ^ 2
  } 
  
  #mean_var*log(1/a1) is the second most significant source of error
  priv_radius <- sqrt(priv_var / n) * qt(1 - a0 / 2, n - 1) + mean_var * log(1 / a1)
  
  return(c(priv_mean - priv_radius, priv_mean + priv_radius))
}

## calls priv_vadhan
priv_vadhan_ci <- function(rep, db, a, e, stdmin, stdmax, xmin, xmax) {
  return(priv_vadhan(db, a / 4, a / 4, a / 4, a / 4, e / 3, e / 3, e / 3, stdmin, stdmax, max(abs(xmax), abs(xmin))))
}



