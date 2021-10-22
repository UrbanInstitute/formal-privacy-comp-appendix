##################################################################
##################################################################
### This function outputs augmented matrix associated with A
##################################################################
##################################################################
# INPUTS:
### Names:    list with variables names
### Type:     list with variables types (numeric or factor) 
### Bounds:   list with upper and bounds for numerical variables
###           and levels for categorical variable
### Ref:      list with reference levels for categorical predictors
### A:        dataframe with predictors and response
### epsilon:  privacy parameter
### delta:  privacy parameter
# OUTPUTS:
### augmented matrix associated with A
# CB: == Claire Bowen notes during code review

# CB: A_WKL() code is the same as in `DPBootsWMFunctions_new.R`
A_WKL = function(A, Names, Type, Bounds, Ref)
{
  ### Applying hard-thresholding to numerical variables
  for(j in 1:ncol(A))
  {
    if(Type[Names[j]] == "numeric")
    {
      # CB: Forcing the variable to be a numeric
      A[, Names[j]] = as.numeric(A[,Names[j]])
      A[, Names[j]] = ifelse(A[,Names[j]] < Bounds[[Names[j]]][1], Bounds[[Names[j]]][1], A[,Names[j]])
      A[, Names[j]] = ifelse(A[,Names[j]] > Bounds[[Names[j]]][2], Bounds[[Names[j]]][2], A[,Names[j]])
    }
  }
  
  ### Removing and adding levels to categorical predictors
  addLevel = function(x, newlevel = NULL) {
    if(is.factor(x)) {
      if (is.na(match(newlevel, levels(x))))
        return(factor(x, levels = c(levels(x), newlevel)))
    }
    return(x)
  }
  
  # CB: Checking if any of the variables are categorical predictors and
  # correcting the levels if some values are dropped
  for(j in 1:ncol(A))
  {
    if(Type[Names[j]] == "factor")
    {
      # CB: Forcing the variable to be a factor
      A[,Names[j]] = as.factor(A[, Names[j]])
      
      # CB: Figuring out which values are within the bounds for that variable
      # by setting them to be true
      Ind = A[, Names[j]] %in% Bounds[[Names[j]]]
      
      # CB: If data are set up correctly, then all data values should be included.
      A = A[Ind, ]
      
      # CB: using addLevels() and droplevels() to adjust the factor levels.
      A[,Names[j]] = droplevels(A[, Names[j]])
      
      # CB: Checking if there are new levels due to the bound adjustments.
      Newlevels = Bounds[[Names[j]]][!(Bounds[[Names[j]]] %in% levels(A[, Names[j]]))]
      if(length(Newlevels) != 0)
        A[, Names[j]] = addLevel(A[, Names[j]], Newlevels)
      A[, Names[j]] = relevel(A[, Names[j]], ref = as.character(Ref[Names[j]]))
    }
  }
  
  ### If the final A (after thresholding and removing and adding levels)
  ### is empty, the following code creates a fake dataset with single row
  # CB: Tested with A <- A[0,], code below generates random values for that row
  if(nrow(A) == 0)
  {
    for(j in 1:length(Names))
    {
      if(Type[Names[j]] == "numeric")
        A[1,Names[j]] = runif(1,Bounds[[Names[j]]][1], Bounds[[Names[j]]][2])
      if(Type[Names[j]] == "factor")
        A[1,Names[j]] = sample(Bounds[[Names[j]]], 1)
    }
  }
  
  ### Rescale Data
  # CB: Normalize data between 0 to 1
  # CB: Tested with apply(A, 2, range), lower bound 0, upper bounds vary
  for(j in 1:ncol(A))
  {
    if(Type[Names[j]] == "numeric")
    {
      A[, Names[j]] = (A[, Names[j]] - Bounds[[Names[j]]][1]) / diff(Bounds[[Names[j]]])
    }
  }
  
  A
}


##################################################################
##################################################################
### This function samples from spherical laplace
##################################################################
##################################################################
# INPUTS:
### n:  sample size
### p:  dimension
### input.par.gama: parameter
# OUTPUTS:
### draws

# CB: Spherical Laplace mechanism draws
# Equation 3 from https://arxiv.org/pdf/1804.03794.pdf
rsphelaplace = function(n, p, input.par.gamma)
{
  
  return(rbind(t(sapply(1:n, function(i){z = rnorm(p); z = z / sqrt(sum(z^2)); z * rgamma(1, p, input.par.gamma)}))))
  
  if(FALSE)
  {
    nsave = 5000 + n * 100
    B = matrix(NA, nsave, p)
    Z = c()
    b = rep(0, p)
    lp = c()
    
    for(J in 1:nsave)
    {
      # logpost
      logpost = function(b) - input.par.gamma*sqrt(sum(b^2))
      
      x1 = x0 = b
      px0 = length(x0)
      if(J < 1000) w = rep(input.par.gamma,px0)
      
      if(J > 1000 & J < 2000) w = apply(apply(B[1:(J - 1), ], 2, range), 2, diff)
      if(J > 2000 & J < 3000) w = apply(apply(B[1000:(J - 1),], 2, range), 2, diff)
      
      ## Slice sampler
      # Step a)
      fx0 = log(runif(1))+logpost(x0)
      
      # Step b)
      u = runif(px0)
      
      LL = x0 - w * u
      RR = LL + w
      
      # Step c)
      x1 = sapply(1:length(x0), function(i1) runif(1, LL[i1], RR[i1]))
      
      z = 0
      while( fx0 >= logpost(x1))
      {
        LL = ifelse(x1 < x0, x1, LL)
        RR = ifelse(x1 >= x0, x1, RR)
        x1 = sapply(1:length(x0), function(i1) runif(1, LL[i1], RR[i1]))
        z = z + 1
      }
      b = x1
      
      B[J,] = b 
      Z[J] = z
      lp[J] = logpost(b)
    }
    B[5000 + (1:n) * 100, ]
  }
}

##################################################################
##################################################################
### This function outputs a noisy version of positive-definite
### matrix (Algorithm 2)
##################################################################
##################################################################
# INPUTS:
### M: positive-definite matrix to be privatized
### Sens: Sensitivity of M under bounded DP
### phi:  privacy parameter
### input.par.c: regularization 
### typeDP: "epsilon-DP" or "rho-zCDP"
# OUTPUTS:
### Noisy version of M

#CB: Lines 1 - 8 of Algorithm 2
Algorithm2_WKL_a = function(M, Sens, phi, typeDP, ndraws = 1)
{
  input.par.c = 0.001 # Default parameter used in the paper
  
  # CB: Conversion from zCDP to eps-DP
  if(typeDP == "epsilon-DP")
  {
    epsilon = phi
    input.par.gamma = phi / Sens
    b = rbind(rsphelaplace(ndraws, prod(dim(M)), input.par.gamma))
  }
  if(typeDP == "rho-zCDP")
  {
    rho = phi
    s = sqrt((Sens^2) / (2 * phi)) 
    b = matrix(rnorm(prod(dim(M)) * ndraws, 0, s), nrow = ndraws)
  }
  ftmp = function(i)
  {
    eta = matrix(b[i, ], ncol(M))
    tildeM = M + eta
    tildeM = (tildeM + t(tildeM)) / 2
    
    colnames(tildeM) = rownames(tildeM) = colnames(M)
    
    tildeM
  }
  lapply(1:ndraws, ftmp)
}

#CB: Lines 9 - 14 of Algorithm 2
Algorithm2_WKL_b = function(M)
{
  input.par.c = 0.001 # Default parameter used in the paper
  ftmp = function(tildeM)
  {
    Gamma = diag(eigen(tildeM)$val)
    V = eigen(tildeM)$vector
    
    Gamma = diag(diag(ifelse(Gamma < 2 * input.par.c, 2 * input.par.c,Gamma)))
    
    tildeM = V %*% Gamma %*% t(V)
    
    tildeM = (tildeM + t(tildeM)) / 2
    
    colnames(tildeM) = rownames(tildeM) = colnames(M[[1]])
    
    tildeM
  }
  lapply(M, ftmp)
}

##################################################################
##################################################################
### This function computes differentially private regression 
### coefficients an confidence intervals using Ferrando's approach and
## algorithm 2  in
### Wang, Yue, Daniel Kifer, and Jaewoo Lee. 2019. 
### Differentially Private Confidence Intervals for Empirical Risk Minimization� 
### Journal of Privacy and Confidentiality 9 (1). https://doi.org/10.29012/jpc.660.
##################################################################
##################################################################
# INPUTS:
### Dataset:  confidential data (data.frame)
### Names:    list with variables names (last variable should be
###           be the response)
### Type:     list with variables types (numeric or factor) 
### Bounds:   list with upper and bounds for numerical variables
###           and levels for categorical variable
### Ref:      list with reference levels for categorical predictors
### epsilon,delta:      privacy parameters
### typeDP: "epsilon-DP" or "rho-zCDP"
### alpha:    significance level
# OUTPUTS:
### asymptotic: DP estimates using asymptotic+plug-in argument 
### bootstrap:  DP estimates using bootstrap
### bootsBeta:  draws for the bootstrap

# CB: The main function called to produce the outputs for the results
DpWKL.lm = function(Dataset, Names, Type, Bounds, Ref, epsilon, delta, typeDP, alpha)
{
  if(typeDP == "epsilon-DP") phi = epsilon
  
  # CB: Conversion from zCDP to eps-DP, there is a pi for spherical data
  if(typeDP == "rho-zCDP")
  {
    compute_epsilon = function(rho, delta)
      rho + 2 * sqrt(rho * log(sqrt(pi * rho) / delta))
    
    rho = uniroot(function(rho) 
      (epsilon - compute_epsilon(rho, delta = delta)), 
      c(delta^2 / pi, 1e10))$root
    phi = rho
  }
  
  ##################
  ### Compute augmented matrix A
  A = A_WKL(Dataset, Names, Type, Bounds, Ref)
  # CB: Create the design matrix
  A = model.matrix(~., data = A)
  # CB: Normalizing again
  A = A / sqrt(ncol(A))
  # CB: ATA (transpose)
  AtA = t(A) %*% A
  
  ##### Noisy statistics
  Sens = 1
  # CB: Adding noise to the design matrix via Algorithm 2 in two functions
  NAtA0 = Algorithm2_WKL_a(AtA,Sens, phi = phi, typeDP)
  NAtA = Algorithm2_WKL_b(NAtA0)[[1]]
  colnames(NAtA) = rownames(NAtA) = colnames(AtA)
  # CB: Number of predictor coefficients, including intercept
  p = ncol(NAtA) - 1
  
  tilde.tXX = NAtA[1:p, 1:p]
  hat.n = NAtA[1, 1] #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  hat.Q = tilde.tXX / hat.n
  
  ##### Regularization
  # CB: Not sure if definiting R is necessary since it's full of 0's.
  R = diag(0, p + 1) # We should not need any regularization using this algorithm 
  # CB: Flagging that this is hardcoded assuming the last variable is the output
  # Values used for Algorithm 3 bootstrap method
  RNAtA = NAtA + R
  Rtilde.tXX = RNAtA[1:p, 1:p]
  Rhat.n = RNAtA[1, 1] * ncol(A)
  Rhat.Q = Rtilde.tXX / Rhat.n
  
  ##### Asymptotic Beta and sigma2 (Plug-in)
  hat.Beta = solve(Rtilde.tXX) %*% RNAtA[1:p, p + 1]
  
  # CB: Fitting a model based on the normalized data without the x.intercept
  # Maybe need to go back and make this not hard coded in.
  #fit = lm(logINCWAGE~.-1, data = data.frame(A))
  
  # CB: Calculating the standard error
  hat.sigma2 = ((RNAtA[p + 1,p + 1] - RNAtA[1:p, p + 1] %*%
                   solve(Rtilde.tXX) %*%
                   RNAtA[1:p, p + 1]) / (Rhat.n - ncol(RNAtA) - 1))[1, 1]
  
  # CB: Note that `alpha()` is a function normally.
  hat.l = hat.Beta - qnorm(1 - alpha / 2) * sqrt(hat.sigma2 * diag(solve(Rtilde.tXX)))
  hat.u = hat.Beta + qnorm(1 - alpha / 2) * sqrt(hat.sigma2 * diag(solve(Rtilde.tXX)))
  
  hat.se = sqrt(hat.sigma2*diag(solve(Rtilde.tXX)))
  
  
  ##### Bootstrap
  # We approximate added noisy using this mechanism multiple times
  NNAtA = Algorithm2_WKL_a(NAtA, Sens, phi = phi,typeDP,ndraws = 10000)
  Noise = list()
  for(iter in 1:10000) Noise[[iter]] = NNAtA[[iter]] - NAtA
  
  Rhat.n0 = Rhat.n
  Rhat.Q0 = Rhat.Q
  
  # CB: fboots for Ferrando Bootstrap method, Algorithm 3
  fboots = function(iter, Noise, Reg)
  {
    N = Noise[[iter]]
    
    tilde.inv_n.V = N[1:p,1:p] # Noise for tXX
    tilde.inv_n.V = tilde.inv_n.V / Rhat.n0
    tilde.inv_n.W = N[1:p, p + 1] / Rhat.n0
    tilde.inv_n.tXu = rmvn(1, mu = rep(0, p), Sigma = hat.sigma2 * Rhat.Q0 / Rhat.n0)
    
    if(Reg == 1)
      tilde.Beta = solve(Rhat.Q) %*% (Rhat.Q - tilde.inv_n.V) %*% hat.Beta + # there is with the dimensions problem Algorithm 3
      solve(Rhat.Q) %*% t(tilde.inv_n.tXu+tilde.inv_n.W)
    
    # CB: For when the error is large.
    if(Reg == 2)
      tilde.Beta = solve(Algorithm2_WKL_b(list(Rhat.Q))[[1]]) %*% Algorithm2_WKL_b(list(Rhat.Q - tilde.inv_n.V))[[1]] %*% hat.Beta + # there is with the dimensions problem Algorithm 3
      solve(Rhat.Q) %*% t(tilde.inv_n.tXu+tilde.inv_n.W)
    
    frescale(tilde.Beta, Names, Bounds, Type)
  }
  
  # CB: apply the fboots() on the noise generated
  tilde.BETA = sapply(1:10000, function(iter) try(fboots(iter,Noise, Reg = 1)))
  tilde.Beta = rowMeans(tilde.BETA)
  hat.Beta = frescale(hat.Beta, Names, Bounds, Type)
  Error = (tilde.Beta - hat.Beta) / hat.Beta
  if(max(abs(Error)) > 1.5) 
    tilde.BETA = sapply(1:10000, function(iter) try(fboots(iter, Noise, Reg = 2)))
  
  # CB: Calculating the CI
  tilde.l = apply(tilde.BETA, 1, quantile,prob = (alpha / 2))
  tilde.u = apply(tilde.BETA, 1, quantile,prob = (1 - alpha / 2))
  tilde.se = apply(tilde.BETA, 1, sd)
  
  # Outputs
  # CB: Rescaling the values, hat.Beta already rescaled earlier in code
  hat.l = frescale(hat.l, Names, Bounds, Type)
  hat.u = frescale(hat.u, Names, Bounds, Type)
  hat.se = (hat.u - hat.l) / (2 * qnorm(1 - alpha / 2))
  asymptotic = data.frame(Beta = hat.Beta, CI.l = hat.l, CI.u = hat.u, se = hat.se)
  
  # asymptotic[1,-1] = NA 
  boots = data.frame(Beta = tilde.Beta, CI.l = tilde.l, CI.u = tilde.u, se = tilde.se)
  rownames(boots) = rownames(asymptotic)
  bootsBeta = t(tilde.BETA)
  colnames(bootsBeta) = row.names(asymptotic)
  
  list(asymptotic = asymptotic,
       bootstrap = boots,
       bootsBeta = bootsBeta)
}

#######################
# This function scale back the estimates to the original scale
#######################
# CB: rescale the coefficients after having normalized the data
frescale = function(hat.Beta, Names, Bounds, Type)
{
  d = length(Names)
  ly = Bounds[[d]][1]
  
  uy_ly = diff(Bounds[[d]])
  rs.hat.Beta = uy_ly * hat.Beta
  aux.hat.Beta = rs.hat.Beta - rs.hat.Beta
  
  lx = Bounds
  ux_lx = Bounds
  for(j in 1:(d - 1))
  {
    lx[[j]] = NA
    ux_lx[[j]] = NA
    if(Type[Names[j]] == "numeric")
    {
      lx[Names[j]] = Bounds[[j]][1]
      ux_lx[Names[j]] = diff(Bounds[[j]])
      rs.hat.Beta[Names[j], ] = rs.hat.Beta[Names[j], ] / ux_lx[Names[j]][[1]]
      aux.hat.Beta[Names[j], ] = rs.hat.Beta[Names[j], ] * lx[Names[j]][[1]] / ux_lx[Names[j]][[1]]
    }
  }
  
  if(rownames(rs.hat.Beta)[1] == "(Intercept)")
    rs.hat.Beta[1, ] = rs.hat.Beta[1, ] + ly - sum(aux.hat.Beta)
  
  rs.hat.Beta
}