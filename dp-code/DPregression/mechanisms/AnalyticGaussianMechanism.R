# n = 10 # Number of draws
# Delta = 1 # Sensitivity
# epsilon = 0.1 # Privacy parameter
# delta = exp(-5) # Privacy parameter
# 
rAnalyticGaussianMechanism = function(n,Delta,epsilon,delta)
{
  delta0 = pnorm(0) - exp(epsilon + pnorm(-sqrt(2*epsilon),log.p = T))
  # Compute alpha
  if(delta >= delta0)
  {
    Bp = function(v,epsilon)
      pnorm(sqrt(epsilon*v)) - exp(epsilon+pnorm(-sqrt(epsilon*(v+2)),log.p = T))
    aux = 0.001
    while (Bp(aux,epsilon) < delta) aux = aux + 0.001
    v.star = optim(1,function(v) abs(Bp(v,epsilon)-delta), method = "Brent", lower = 0, upper = aux)$par
    Bp(v.star,epsilon) <= delta
    alpha = sqrt(1+(v.star/2))-sqrt(v.star/2)
  }
  if(delta < delta0)
  {
    Bm = function(u,epsilon)
      pnorm(-sqrt(epsilon*u)) - exp(epsilon+pnorm(-sqrt(epsilon*(u+2)),log.p = T))
    aux = 0
    while (Bm(aux,epsilon) >= delta) aux = aux + 0.001
    u.star = optim(1,function(u) abs(Bm(u,epsilon)-delta), method = "Brent", lower = aux, upper =  aux + 0.001)$par
    Bm(u.star,epsilon) <= delta
    alpha = sqrt(1+(u.star/2))+sqrt(u.star/2)
  }
  rnorm(n,0,alpha*Delta/sqrt(2*epsilon))
}
