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

A_WKL = function(A, Names, Type, Bounds, Ref)
{
  ### Applying hard-thresholding to numerical varaibles
  for(j in 1:ncol(A))
  {
    if(Type[Names[j]] == "numeric")
    {
      A[,Names[j]] = as.numeric(A[,Names[j]])
      A[,Names[j]] = ifelse(A[,Names[j]] < Bounds[[Names[j]]][1], Bounds[[Names[j]]][1], A[,Names[j]])
      A[,Names[j]] = ifelse(A[,Names[j]] > Bounds[[Names[j]]][2], Bounds[[Names[j]]][2], A[,Names[j]])
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
  for(j in 1:ncol(A))
  {
    if(Type[Names[j]] == "factor")
    {
      A[,Names[j]] = as.factor(A[,Names[j]])
      Ind = A[,Names[j]] %in% Bounds[[Names[j]]]
      A = A[Ind,]
      A[,Names[j]] = droplevels(A[,Names[j]])
      Newlevels = Bounds[[Names[j]]][!(Bounds[[Names[j]]]%in%levels(A[,Names[j]]))]
      if(length(Newlevels) != 0)
        A[,Names[j]] = addLevel(A[,Names[j]],Newlevels)
      A[,Names[j]] = relevel(A[,Names[j]], ref = as.character(Ref[Names[j]]))
    }
  }
  ### If the final A (after thresholding and removing and adding levels)
  ### is empty, the following code creates a fake dataset with single row
  if(nrow(A) == 0)
  {
    for(j in 1:length(Names))
    {
      if(Type[Names[j]] == "numeric")
        A[1,Names[j]] = runif(1,Bounds[[Names[j]]][1],Bounds[[Names[j]]][2])
      if(Type[Names[j]] == "factor")
        A[1,Names[j]] = sample(Bounds[[Names[j]]],1)
    }
  }
  
  ### Rescale Data
  for(j in 1:ncol(A))
  {
    if(Type[Names[j]] == "numeric")
    {
      A[,Names[j]] = (A[,Names[j]] - Bounds[[Names[j]]][1])/diff(Bounds[[Names[j]]])
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
rwhishartcentered = function(n,p,epsilon,delta)
{
  # Generating Wishart error
  B = sqrt(1) # Bound on the l2-norm of A = (Y,X)
  d = p  # Number of columns in A
  k = floor(d+(14/epsilon^2)*2*log(4/delta))
  mW0 = k*B^2*diag(1,d)
  
  return(rbind(t(sapply(1:n,function(i) as.vector(rWishart(1,k,B*diag(d))[,,1]-mW0)))))
  if(FALSE)
  {
    # Generating error
    ftmp = function(iter)
    {
      v = matrix(rnorm(k*d,0,sqrt(B)),nrow=k,ncol=d) 
      Error = t(v)%*%v
      as.vector(Error-mW0)
    }
    t(sapply(1:n,ftmp))
  }
  
}

##################################################################
##################################################################
### This function outputs a noisy version of positive-definite
### matrix (Algorithm 2) Sheffet + regularization
### Upadhyay, Jalaj, and Sarvagya Upadhyay. 
### Sheffet, Or. "Old techniques in differentially private linear regression." 
### In Algorithmic Learning Theory, pp. 789-827. PMLR, 2019.
##################################################################
##################################################################
# INPUTS:
### AtA: positive-definite matrix to be privatized
### epsilon and delta:  privacy parameters
# OUTPUTS:
### Noisy version of AtA
Algorithm2_Whishart = function(AtA, epsilon,delta,ndraws = 1)
{
  d = ncol(AtA)
  b = rbind(rwhishartcentered(ndraws,d,epsilon,delta))
  ftmp = function(i)
  {
    eta = matrix(b[i,],d)
    NAtA = AtA + eta
    
    d = d  # Number of columns in A
    B = 1 # Bound on the l2-norm of A = (Y,X)
    k = floor(d+(14/epsilon^2)*2*log(4/delta))
    mW0 = k*B^2*diag(1,d)
    mW1 = (B^2)*(sqrt(k)-(sqrt(d)+sqrt(2*log(4/delta))))^2*diag(1,d)
    
    ##### Regularization
    R = diag(0,d)
    
    # regularization proportional equal to quantile of the
    if(!is.positive.definite(NAtA+R))
      R = -mW1+mW0
    
    # regularization proportonal equal to 3 times the eigen value of NAtA
    if(!is.positive.definite(NAtA+R))
      R = diag(-3*min(eigen(NAtA)$val),d)
    
    ##### Noisy statistics Beta
    RNAtA = NAtA+R
  }
  lapply(1:ndraws,ftmp)
}


##################################################################
##################################################################
### This function computes differentially private regression 
### coefficients, confidence intervals, and standard errors using 
### the Wishart mechanism.
### This implementation is an adaptation of the approach proposed
### by Ferrando, Cecilia, Shufan Wang, and Daniel Sheldon. 
### "General-Purpose Differentially-Private Confidence Intervals." 
### arXiv preprint arXiv:2006.07749 (2020).
##################################################################
##################################################################
# INPUTS:
### Names:    list with variables names (last variable should be
###           be the response)
### Type:     list with variables types (numeric or factor) 
### Bounds:   list with upper and bounds for numerical variables
###           and levels for categorical variable
### Ref:      list with reference levels for categorical predictors
### Dataset:  confidential data (data.frame)
### epsilon:  privacy parameter
### delta:    privacy parameter
### alpha:    significance level
# OUTPUTS:
### asymptotic: DP estimates using asymptotic+plug-in argument 
### bootstrap:  DP estimates using bootstrap
### bootsBeta:  draws for the bootstrap
DpBoots_WM.lm = function(Dataset, Names, Type, Bounds, Ref, epsilon, delta, alpha)
{
  ##################
  ### Compute augmented matrix A
  A = A_WKL(Dataset, Names, Type, Bounds, Ref)
  summary(A)
  A = model.matrix(~., data = A)
  summary(A)
  A = A /sqrt(ncol(A))
  range(apply(A, 1, function(x) sqrt(sum(x^2))))
  
  AtA = t(A)%*%A
  
  ##### Noisy statistics
  NAtA = Algorithm2_Whishart(AtA,epsilon,delta,ndraws = 1)[[1]]
  colnames(NAtA) = rownames(NAtA) = colnames(AtA)
  p = ncol(NAtA)-1
  
  round(NAtA/AtA,5)
  
  tilde.tXX = NAtA[1:p,1:p]
  hat.n = NAtA[1,1]*ncol(A) #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  nrow(A);hat.n
  
  hat.Q = tilde.tXX/hat.n
  
  ##### Regularization
  RNAtA = NAtA
  Rtilde.tXX = RNAtA[1:p,1:p]
  Rhat.n = RNAtA[1,1]*ncol(A)
  Rhat.Q = Rtilde.tXX/Rhat.n
  
  ##### Asymptotic Beta and sigma2 (Plug-in)
  hat.Beta = solve(Rtilde.tXX)%*%RNAtA[1:p,p+1]
  
  hat.Beta0 = hat.Beta
  
  hat.Beta; solve(RNAtA[1:p,1:p])%*%RNAtA[1:p,p+1]
  RNAtA[1:p,1:p]/AtA[1:p,1:p]
  solve(RNAtA[1:p,1:p])/solve(AtA[1:p,1:p])
  RNAtA[1:p,p+1]/AtA[1:p,p+1]
  
  
  hat.sigma2 = ((RNAtA[p+1,p+1]-RNAtA[1:p,p+1]%*%
                   solve(Rtilde.tXX)%*%
                   RNAtA[1:p,p+1])/(Rhat.n-ncol(RNAtA)-1))[1,1]
  
  hat.l = hat.Beta-qnorm(1-alpha/2)*sqrt(hat.sigma2*diag(solve(Rtilde.tXX)))
  hat.u = hat.Beta+qnorm(1-alpha/2)*sqrt(hat.sigma2*diag(solve(Rtilde.tXX)))
  
  hat.se = sqrt(hat.sigma2*diag(solve(Rtilde.tXX)))
  
  # fit = summary(fit)
  # Original = cbind(fit[["coefficients"]][,1],confint(lm(logINCWAGE~.-1,data=data.frame(A))),fit[["coefficients"]][,2])
  # colnames(Original) = c("Beta","CI.l","CI.u","se")
  # 
  # Original
  # 
  
  ##### Bootstrap
  # We approximate added noisy using this mechanism multiple times
  Noise = list()
  for(iter in 1:10000) Noise[[iter]] = matrix(rwhishartcentered(1,p+1,epsilon,delta),p+1,p+1)
  
  fboots = function(iter,Noise)
  {
    # N1 = rlaplace(length(diagAtA),0,DeltaAtA/epsilon)
    # N2 = rlaplace(length(unique_offdiagAtA),0,DeltaAtA/epsilon)
    # N = diag(N1)
    # for(j in 1:length(unique_offdiagAtA))
    # N  = (AtA==unique_offdiagAtA[j])*N2[j] + N
    N = Noise[[iter]]
    
    tilde.inv_n.V = N[1:p,1:p] # Noise for tXX
    tilde.inv_n.V = tilde.inv_n.V/Rhat.n
    tilde.inv_n.W = N[1:p,p+1]/Rhat.n
    tilde.inv_n.tXu = rmvn(1,mu = rep(0,p), Sigma = hat.sigma2*Rhat.Q/Rhat.n)
    
    tilde.Beta = solve(Rhat.Q)%*%(Rhat.Q - tilde.inv_n.V)%*%hat.Beta + # there is with the dimensions problem Algorithm 3
      solve(Rhat.Q)%*%t(tilde.inv_n.tXu+tilde.inv_n.W)
    
    frescale(tilde.Beta,Names,Bounds,Type)
  }
  tilde.BETA = sapply(1:10000,function(iter) try(fboots(iter,Noise)))
  
  if(is.list(tilde.BETA))
  {
    Ind = which(simplify2array(lapply(tilde.BETA,is.numeric)))
    tilde.BETA = sapply(Ind,function(iter) tilde.BETA[[iter]])
  }
  
  tilde.Beta = rowMeans(tilde.BETA)
  tilde.l = apply(tilde.BETA,1,quantile,prob = (alpha/2))
  tilde.u = apply(tilde.BETA,1,quantile,prob = (1-alpha/2))
  tilde.se = apply(tilde.BETA,1,sd)
  
  # Outputs
  hat.Beta = frescale(hat.Beta,Names,Bounds,Type)
  hat.l = frescale(hat.l,Names,Bounds,Type)
  hat.u = frescale(hat.u,Names,Bounds,Type)
  hat.se = (hat.u - hat.l)/(2*qnorm(1-alpha/2))
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
frescale = function(hat.Beta,Names,Bounds,Type)
{
  d = length(Names)
  rescale1 = rep(diff(Bounds[[d]]),nrow(hat.Beta))
  names(rescale1) = rownames(hat.Beta)
  rescale2 = rep(0,nrow(hat.Beta))
  names(rescale2) = rownames(hat.Beta)
  for(j in 1:(d-1))
  {
    if(Type[Names[j]] == "numeric")
    {
      rescale1[Names[j]] = rescale1[Names[j]]/diff(Bounds[[Names[j]]])
      rescale2[Names[j]] = -diff(Bounds[[Names[d]]])*Bounds[[Names[j]]][1]/diff(Bounds[[Names[j]]])
    }
  }
  
  aux = hat.Beta*rescale1
  aux[1] = (hat.Beta*rescale1)[1]+sum((hat.Beta*rescale2))+Bounds[[Names[d]]][1]
  
  aux
}

