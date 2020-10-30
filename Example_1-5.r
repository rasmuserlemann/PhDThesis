#Conditional simulations
simcount = function(n,t)
{
  barsvec=randombars(n,t)
  tell=c()
  mm=max(barsvec)
  for (k in 0:mm)
  {
    tell[k+1]=0
    for (i in 1:n)
    {
      tell[k+1] = tell[k+1] + 1*(barsvec[i]==k)
    }
  }
  return(tell)
}

# Draws from the stars and bars
#
alg1 = function(n, t){
  N = 0
  V = 0
  I = t+n-1
  k = c()
  while (N<n-1){
    p = runif(1)
    if (p<(n-1-N)/(t+n-1-N-V)){
      k = c(I,k)
      N = N+1
      I = I-1
    }
    else{
      V = V+1
      I = I-1
    }
  }
  return(k)
}
randombars = function(n,t){
  s1 = alg1(n, t)
  final = c(s1[1]-1)
  for (ind in 2:length(s1)){
    val = s1[ind]-ind-(s1[ind-1]-(ind-1))
    final = append(final, val)
  }
  final = append(final, t-(s1[length(s1)]-(n-1)))
  return (final)
}

# Function that takes 'd' to 'count' 
ct = function(dat){
  out = c()
  M = max(dat)
  for (v in 0:M)
    out = c(out, sum(dat == v))
  out
}

#Cramer-von Mises test
CM <- function(count)
{
  m=length(count)
  n=sum(count) # number of obs
  vec=c(0:(m-1))
  t=sum(vec*count) #sum of xi
  phat = n/(n+t)
  O=cumsum(count)
  E=c()
  ph=c()
  for (k in 1:m) E[k]=n*(1-(1-phat)^k)
  for (k in 1:m) ph[k]=phat*(1-phat)^(k-1)
  Z=O-E
  mm = m
  while (ph[mm]*(1-phat) >= 0.001/n)
  {
    ph[mm+1] = ph[mm]*(1-phat)
    Z[mm+1] = Z[mm]-n*ph[mm+1]
    mm = mm+1
  }
  CM=sum(Z^2*ph)/n
  return(CM)
}

n = 100
d1 = rgeom(n, 0.5)
d2 = rdweibull(n,0.7,0.8,zero=TRUE)
datacount1 = ct(d1)
datacount2 = ct(d2)
#Critical values
cramer1 = CM(datacount1)
cramer2 = CM(datacount2)

#Histograms
hist(d1, xlab="", ylab = "", main = "")
hist(d2, xlab="", ylab = "", main = "")

#Sufficient statistic values
t1 = sum(d1)
t2 = sum(d2)

#MLE estimates under the null hypothesis
mle1 = n/(t1+n)
mle2 = n/(t2+n)

#Cramer-von Mises test distribution via parametric bootstrapping
test1 = c()
test2 = c()
L = 10**(5)
for (k in 1:L){
  dp1 = rgeom(n, mle1)
  dp2 = rgeom(n, mle2)
  c1 = ct(dp1)
  c2 = ct(dp2)
  test1 = c(test1, CM(c1))
  test2 = c(test2, CM(c2))
}
plot(density(test1), ylim = c(0,11.7), xlim = c(0,0.5), xlab="", ylab = "", main = "")
plot(density(test2), ylim = c(0,9.1), xlim = c(0,0.9), xlab="", ylab = "", main = "")

#Conditional test distributions
condtest1 = c()
condtest2 = c()
for (k in 1:L){
  condc1 = simcount(n, t1)
  condc2 = simcount(n, t2)
  condtest1 = c(condtest1, CM(condc1))
  condtest2 = c(condtest2, CM(condc2))
}
par(new=TRUE)
plot(density(condtest1), ylim = c(0,11.7), xlim = c(0,0.5), xlab="", ylab = "", main = "", col = "blue")
plot(density(condtest2), ylim = c(0,9.1), xlim = c(0,0.9), xlab="", ylab = "", main = "", col = "blue")
abline(v=cramer1, col="green")
abline(v=cramer2, col="green")
