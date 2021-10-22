######
# This functions adapts the code for the algorithm proposed by Brawner and Honaker, as
# discussed in Section 3.3 in the Du paper. 
######

compute_epsilon <- function(rho, delta) {
  rho + 2 * sqrt(rho * log(sqrt(pi * rho)/delta))
}

define_p_i <- function(i, N) {
  choose(N, i) * ((1 / N) ^ i) * (1 - 1 / N) ^ (N - i)
}

partition_bootstrap <- function(A, p_i) {
  N <- nrow(A)
  weights <- rmultinom(1, N, rep(1 / N, N))
  partitions = list()
  for(a in 1:length(p_i)){
    partitions[[a]] = which(weights == a)
  }
  return(partitions)
}

define_sigma_i <- function(sensitivity, i, rho, N) {
  p_i <- define_p_i(i, N)
  i * p_i * (sensitivity ^ 2 / (2 * rho))
}


##################################################################
##################################################################
### This function outputs the l2 sensitivity rescaling variables
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
### Sensitivy

Sensitivity_BHM = function(A, Names, Type, Bounds, Ref, epsilon, delta)
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
  ### Output
  # AtA matrix
  a = model.matrix(as.formula(paste("~",paste(Names,collapse = "+"))),data = A)
  AtA = matrix(NA, ncol = ncol(a), nrow = ncol(a))
  colnames(AtA) = rownames(AtA) = colnames(a)
  # Sensitivy for AtA entries
  SensitivyAtA = AtA
  # Epsilon used at each entry of AtA
  EpsilonSplitAtA = AtA
  count.epsilon = 0
  # Entries to be sanitized in AtA
  QuerySplitAtA = AtA
  count.queries = 0
  
  ### Computing sensitivy and identifying entries to be sanitized
  # Intercept
  j = 1
  a = model.matrix(as.formula(paste("~",paste(Names[c(j,j)],collapse = "+"))),data = A)
  ata = (t(a)%*%a)[1,1,drop = F]
  AtA[rownames(ata),colnames(ata)] = ata
  SensitivyAtA[rownames(ata),colnames(ata)] = 1
  count.queries = count.queries + 1
  QuerySplitAtA[rownames(ata),colnames(ata)] = count.queries
  count.epsilon = count.epsilon + 1
  EpsilonSplitAtA[rownames(ata),colnames(ata)] = count.epsilon
  
  # Diagonal and first row
  for(j in 1:length(Names))
  {
    if(Type[j] == "numeric")
    {
      a = model.matrix(as.formula(paste("~",paste(Names[c(j,j)],collapse = "+"))),data = A)
      ata = (t(a)%*%a)
      ata = diag(ata)[-1]
      AtA[names(ata),names(ata)] = ata
      
      L = optim(mean(Bounds[[Names[j]]]), function(x) x^2, 
                lower = Bounds[[Names[j]]][1], upper = Bounds[[Names[j]]][2],
                method = "L-BFGS-B")$value
      U = -optim(mean(Bounds[[Names[j]]]), function(x) -x^2, 
                 lower = Bounds[[Names[j]]][1], upper = Bounds[[Names[j]]][2],
                 method = "L-BFGS-B")$value
      
      SensitivyAtA[names(ata),names(ata)] = U-L
      count.queries = count.queries + 1
      QuerySplitAtA[names(ata),names(ata)] = count.queries
      count.epsilon = count.epsilon + 1
      EpsilonSplitAtA[names(ata),names(ata)] = count.epsilon
      
      ata = (t(a)%*%a)[1,2,drop = F]
      AtA[rownames(ata),colnames(ata)] = ata
      SensitivyAtA[rownames(ata),colnames(ata)] = diff(Bounds[[Names[j]]])
      count.queries = count.queries + 1
      QuerySplitAtA[rownames(ata),colnames(ata)] = count.queries
      count.epsilon = count.epsilon + 1
      EpsilonSplitAtA[rownames(ata),colnames(ata)] = count.epsilon
    }
    if(Type[j] == "factor")
    {
      a = model.matrix(as.formula(paste("~",paste(Names[c(j,j)],collapse = "+"))),data = A)
      ata = (t(a)%*%a)
      ata = diag(ata)[-1]
      
      if(length(ata) == 1)
      {
        AtA[names(ata),names(ata)] = ata
        AtA[1,names(ata)] = ata
        
        SensitivyAtA[names(ata),names(ata)] = diag(rep(1,length(ata)))
        SensitivyAtA[1,names(ata)] = 1
        
        count.queries = count.queries + 1
        QuerySplitAtA[names(ata),names(ata)] = count.queries
        QuerySplitAtA[1,names(ata)] = count.queries
        
        count.epsilon = count.epsilon + 1
        EpsilonSplitAtA[names(ata),names(ata)] = count.epsilon
        EpsilonSplitAtA[1,names(ata)] = count.epsilon
      }
      if(length(ata)>1)
      {
        AtA[names(ata),names(ata)] = diag(ata)
        AtA[1,names(ata)] = ata
        
        SensitivyAtA[names(ata),names(ata)] = diag(rep(1,length(ata)))
        SensitivyAtA[1,names(ata)] = 1
        
        count.queries = count.queries + length(ata)
        QuerySplitAtA[names(ata),names(ata)] = diag(count.queries:(count.queries-length(ata)+1))
        QuerySplitAtA[1,names(ata)] = (count.queries:(count.queries-length(ata)+1))
        
        count.epsilon = count.epsilon + 1
        EpsilonSplitAtA[names(ata),names(ata)] = diag(rep(count.epsilon,length(ata)))
        EpsilonSplitAtA[1,names(ata)] = count.epsilon
      }
    }
  }
  
  # Off-diagonal and first row
  J = combn(1:ncol(A),2)
  
  for(j in 1:ncol(J))
  {
    j1 = J[1,j]; j2 = J[2,j]
    a = model.matrix(as.formula(paste("~",paste(Names[c(j1,j2)],collapse = "+"))),data = A)
    ata = (t(a)%*%a)
    rowNames = colnames(model.matrix(as.formula(paste("~",paste(Names[c(j1,j1)],collapse = "+"))),data = A))[-1]
    colNames = colnames(model.matrix(as.formula(paste("~",paste(Names[c(j2,j2)],collapse = "+"))),data = A))[-1]
    AtA[rowNames,colNames] = ata[rowNames,colNames]
    
    if(sum(Type[Names[c(j1,j2)]] == "numeric") == 2)
    {
      L = optim(c(mean(Bounds[[Names[j1]]]),mean(Bounds[[Names[j2]]])), function(x) prod(x), 
                lower = c(Bounds[[Names[j1]]][1],Bounds[[Names[j2]]][1]), 
                upper = c(Bounds[[Names[j1]]][2],Bounds[[Names[j2]]][2]),
                method = "L-BFGS-B")$value
      U = -optim(c(mean(Bounds[[Names[j1]]]),mean(Bounds[[Names[j2]]])), function(x) -prod(x), 
                 lower = c(Bounds[[Names[j1]]][1],Bounds[[Names[j2]]][1]), 
                 upper = c(Bounds[[Names[j1]]][2],Bounds[[Names[j2]]][2]),
                 method = "L-BFGS-B")$value
      SensitivyAtA[rowNames,colNames] = U-L
      count.queries = count.queries + 1
      QuerySplitAtA[rowNames,colNames] = count.queries
      count.epsilon = count.epsilon + 1
      EpsilonSplitAtA[rowNames,colNames] = count.epsilon
    }
    
    if(sum(Type[c(j1,j2)] == "numeric") == 1)
    {
      SensitivyAtA[rowNames,colNames] = diff(Bounds[[Names[c(j1,j2)][Type[Names[c(j1,j2)]] == "numeric"]]])
      count.queries = count.queries + length(ata[rowNames,colNames])
      QuerySplitAtA[rowNames,colNames] = (count.queries:(count.queries-length(ata[rowNames,colNames])+1))
      count.epsilon = count.epsilon + 1
      EpsilonSplitAtA[rowNames,colNames] = count.epsilon
    }
    
    if(sum(Type[c(j1,j2)] == "numeric") == 0)
    {
      SensitivyAtA[rowNames,colNames] = 1
      count.queries = count.queries + length(ata[rowNames,colNames])
      QuerySplitAtA[rowNames,colNames] = (count.queries:(count.queries-length(ata[rowNames,colNames])+1))
      count.epsilon = count.epsilon + 1
      EpsilonSplitAtA[rowNames,colNames] = count.epsilon
    }
  }
  
  AtA = ifelse(is.na(AtA),0,AtA)
  AtA = t(AtA)+AtA*upper.tri(diag(1,ncol(AtA)))
  
  SensitivyAtA = ifelse(is.na(SensitivyAtA),0,SensitivyAtA)
  SensitivyAtA = t(SensitivyAtA)+SensitivyAtA*upper.tri(diag(1,ncol(AtA)))
  
  EpsilonSplitAtA = ifelse(is.na(EpsilonSplitAtA),0,EpsilonSplitAtA)
  EpsilonSplitAtA = t(EpsilonSplitAtA)+EpsilonSplitAtA*upper.tri(diag(1,ncol(AtA)))
  
  QuerySplitAtA = ifelse(is.na(QuerySplitAtA),0,QuerySplitAtA)
  QuerySplitAtA = t(QuerySplitAtA)+QuerySplitAtA*upper.tri(diag(1,ncol(AtA)))
  
  list(sensitivity = sqrt(max(EpsilonSplitAtA)), # Sensitivity - rescale variables
       QuerySplitAtA = QuerySplitAtA,
       SensitivyAtA = SensitivyAtA,
       A = A)
}
##################################################################
##################################################################
### This function outputs each Mi
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
### partititions:  labels for partition
# OUTPUTS:
### AtA:      AtA matrix
### NAtA:     sanitized AtA matrix (Mi)

# Mi_BHM = function(A, Names, Type, Bounds, Ref, epsilon, delta, sigma2, i, partitions,
#                   QuerySplitAtA, SensitivyAtA)
Mi_BHM = function(A, Names, Type, Bounds, Ref, epsilon, delta, sigma2, i, partitions,
                    QuerySplitAtA, SensitivyAtA)
{
  # ### Applying hard-thresholding to numerical varaibles
  # for(j in 1:ncol(A))
  # {
  #   if(Type[Names[j]] == "numeric")
  #   {
  #     A[,Names[j]] = as.numeric(A[,Names[j]])
  #     A[,Names[j]] = ifelse(A[,Names[j]] < Bounds[[Names[j]]][1], Bounds[[Names[j]]][1], A[,Names[j]])
  #     A[,Names[j]] = ifelse(A[,Names[j]] > Bounds[[Names[j]]][2], Bounds[[Names[j]]][2], A[,Names[j]])
  #   }
  # }
  # ### Removing and adding levels to categorical predictors
  # addLevel = function(x, newlevel = NULL) {
  #   if(is.factor(x)) {
  #     if (is.na(match(newlevel, levels(x))))
  #       return(factor(x, levels = c(levels(x), newlevel)))
  #   }
  #   return(x)
  # }
  # for(j in 1:ncol(A))
  # {
  #   if(Type[Names[j]] == "factor")
  #   {
  #     A[,Names[j]] = as.factor(A[,Names[j]])
  #     Ind = A[,Names[j]] %in% Bounds[[Names[j]]]
  #     A = A[Ind,]
  #     A[,Names[j]] = droplevels(A[,Names[j]])
  #     Newlevels = Bounds[[Names[j]]][!(Bounds[[Names[j]]]%in%levels(A[,Names[j]]))]
  #     if(length(Newlevels) != 0)
  #       A[,Names[j]] = suppressWarnings(addLevel(A[,Names[j]],Newlevels))
  #     A[,Names[j]] = relevel(A[,Names[j]], ref = as.character(Ref[Names[j]]))
  #   }
  # }
  # ### If the final A (after thresholding and removing and adding levels)
  # ### is empty, the following code creates a fake dataset with single row
  # if(nrow(A) == 0)
  # {
  #   for(j in 1:length(Names))
  #   {
  #     if(Type[Names[j]] == "numeric")
  #       A[1,Names[j]] = runif(1,Bounds[[Names[j]]][1],Bounds[[Names[j]]][2])
  #     if(Type[Names[j]] == "factor")
  #       A[1,Names[j]] = sample(Bounds[[Names[j]]],1)
  #   }
  # }
  ### Output
  # AtA matrix
  a = model.matrix(as.formula(paste("~",paste(Names,collapse = "+"))),data = A)
  AtA = matrix(NA, ncol = ncol(a), nrow = ncol(a))
  colnames(AtA) = rownames(AtA) = colnames(a)

  AtA = t(model.matrix(~., data = A))%*%model.matrix(~., data = A)
  
  # Define a version of AtA that has entries with sensitivity equal 1 
  AtAsens1 = AtA/SensitivyAtA
  AtAsens1 = ifelse(is.na(AtAsens1),0,AtAsens1)
  
  # Generating error
  Error = AtA-AtA
  error = rnorm(max(QuerySplitAtA),0,sqrt(sigma2))
  for(j in 1:max(QuerySplitAtA))
  {
    Error = ifelse(QuerySplitAtA == j,error[j],Error)
  }
  round(Error,3)
  
  # Noisy version of AtA
  NAtA = (i*AtAsens1+Error)*SensitivyAtA
  
  if(length(partitions[[i]])==0)
    NAtA = (Error)*SensitivyAtA
  
  list(AtA = AtA,
       NAtA = NAtA)
}
##################################################################
##################################################################
### This function computes differentially private regression 
### coefficients, confidence intervals, and standard errors using 
### the Laplace mechanism.
### This implementation is an adaptation of the approach proposed
### by  by Brawner and Honaker, as discussed in Section 3.3 in the Du paper.
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
### NN:        number of bootstrap samples
### alpha:    significance level
# OUTPUTS:
### bootstrap:  DP estimates using bootstrap
### bootstrapAN:  DP estimates using bootstrap + Asymptotic Normality
### bootsBeta:  draws for the bootstrap
DpBoots_BHM.lm = function(Dataset, Names, Type, Bounds, Ref, epsilon, delta, NN, alpha)
{
  ### Data
  A = cbind(Dataset[,Names])
  Epsilon0 = epsilon
  epsilon = Epsilon0/NN
  
  Delta0 = delta
  delta = Delta0/NN
  
  
  Hat.beta = list()
  
  aux = Sensitivity_BHM(A, Names, Type, Bounds, Ref, epsilon=epsilon, delta=delta)
  sensitivity = aux$sensitivity  
  QuerySplitAtA = aux$QuerySplitAtA
  SensitivyAtA = aux$SensitivyAtA
  AA = aux$A
  
  for(J in 1:NN)
  {
    ######
    # Algorithm proposed by Brawner and Honaker, as
    # discussed in Section 3.3 in the Du paper. 
    ######
    
    rho = uniroot(function(rho) (epsilon - compute_epsilon(rho, delta = delta)), c(delta^2/pi,1e10))$root
    p_i = sapply(1:nrow(A),define_p_i, N = nrow(A))
    p_i = p_i[!is.na(p_i)]
    p_i = p_i[p_i != 0]
    
    partitions = partition_bootstrap(A,p_i)
    
    # sensitivity = Sensitivity_BHM(A, Names, Type, Bounds, Ref, epsilon=epsilon, delta=delta)
    
    sigma2_i = c()
    for(i in 1:length(p_i)) 
      sigma2_i[i] = define_sigma_i(sensitivity, i, rho, nrow(A))
    
    compute_M_i = function(i)
    {
      Mi_BHM(AA[partitions[[i]],], 
             Names, Type, Bounds, Ref, epsilon, delta, 
             sigma2=sigma2_i[i], i = i, partitions = partitions, 
             QuerySplitAtA, SensitivyAtA)$NAtA
    }
    M_i = sapply(1:length(partitions),compute_M_i)
    M = matrix(rowSums(M_i),sqrt(nrow(M_i)))
    
    aux = colnames(Mi_BHM(AA[partitions[[1]],], 
                          Names, Type, Bounds, Ref, epsilon, delta, 
                          sigma2=sigma2_i[1], i = 1, partitions = partitions,
                          QuerySplitAtA, SensitivyAtA)$NAtA)
    colnames(M) = rownames(M) = aux
    
    Hat.beta[[J]] = solve(M[-ncol(M),-ncol(M)])%*%M[1:(ncol(M)-1),ncol(M)]
    
  }
  
  hat.BETA = Hat.beta[[1]]
  for(J in 2:length(Hat.beta)) hat.BETA = cbind(hat.BETA,Hat.beta[[J]])
  
  
  tilde.Beta = rowMeans(hat.BETA)
  tilde.l = apply(hat.BETA,1,quantile,prob = (alpha/2))
  tilde.u = apply(hat.BETA,1,quantile,prob = (1-alpha/2))
  tilde.se = apply(hat.BETA,1,sd)
  tilde.l.AN = tilde.Beta - qnorm(1-alpha/2)*tilde.se 
  tilde.u.AN = tilde.Beta +qnorm(1-alpha/2)*tilde.se
  
  # Outputs
  boots = data.frame(Beta = tilde.Beta, CI.l = tilde.l, CI.u = tilde.u, se = tilde.se)
  bootsAN = data.frame(Beta = tilde.Beta, CI.l = tilde.l.AN, CI.u = tilde.u.AN)
  
  list(bootstrap = boots,
       bootstrapAN = bootsAN,
       bootsBeta = hat.BETA)
}

# out = DpBoots_BHM.lm(Dataset, Names, Type, Bounds, Ref, epsilon, delta, NN=10)
# 
# cbind(lm(logincome~.,data=A)$coef,confint(lm(logincome~.,data=A)))
# out$bootstrap
# out$bootstrapAN
