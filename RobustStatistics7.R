library(depth)

# crude implementation of median-of-means estimator
depthMOM = function(x, m = 11) # m is the number of blocks
{ 
  N = length(x[,1])
  NN = N
  y = x
  smeans = matrix(NA, m, length(x[1,]))
  for(i in 1:m)
  { 
    yy <- sample(1:NN, floor(N/m))
    subsample <- y[yy,]
    y = y[-yy,]
    NN = NN - floor(N/m)
    smeans[i,] = colMeans(subsample)
  }
  rest = med(smeans, approx = TRUE) #  "approx" allows approximate solution, the exact might be not possible
  return(rest$median)
}

# a)
N = 1100
R = 200
set.seed(12)
results <- replicate(R, {
  x1 = rnorm(N)
  x2 = rexp(N)
  xm = data.frame(x1, x2)
  smean = colMeans(xm)
  MOM = depthMOM(xm, m = 11)
  c(smean, MOM)
})
#  the warnings are about the small number of iterations, you can increase it
par(mfrow = c(3,2))
boxplot(results[1,],results[3,]); abline(h=0) # horizontal line marks the true value of the parameter 
boxplot(results[2,],results[4,]); abline(h=1)

# b)
set.seed(12)
results <- replicate(R, {
  x1 = rnorm(N)
  x2 = rexp(N)
  x1[N] = 100
  x2[N] = 100
  xm = data.frame(x1, x2)
  smean = colMeans(xm)
  MOM = depthMOM(xm, m = 11)
  c(smean, MOM)
})
# with one heavy outlier, MOM works fine, mean is biased
boxplot(results[1,],results[3,]); abline(h=0) # horizontal line marks the true value of the parameter 
boxplot(results[2,],results[4,]); abline(h=1)

# c)
set.seed(1)
results <- replicate(R, {
  x1 = rnorm(N)
  x2 = rexp(N)
  xm = data.frame(x1, x2)
  eps <- rbinom(N, 1, 0.01)
  xm <- (1-eps)*xm + eps*cbind(rep(3, N), rep(5, N))
  smean = colMeans(xm)
  MOM = depthMOM(xm, m = 11)
  c(smean, MOM)
})
# with epsilon = 0.01 the number of outliers is larger than the MOM can tolerate,
# hence the difference between teh estimators is very small.
boxplot(results[1,],results[3,]); abline(h=0) # horizontal line marks the true value of the parameter 
boxplot(results[2,],results[4,]); abline(h=1)
#

#---------------------------------------------------------------------------------------------------
# we can make a few more runs 
set.seed(1)
results <- replicate(R, {
  x1 = rnorm(N)
  x2 = rexp(N)
  xm = data.frame(x1, x2)
  eps <- rbinom(N, 1, 0.005)
  xm <- (1-eps)*xm + eps*cbind(rep(3, N), rep(5, N))
  smean = colMeans(xm)
  MOM = depthMOM(xm, m = 22)
  c(smean, MOM)
})
# MOM works better
par(mfrow = c(3,2))
boxplot(results[1,],results[3,]); abline(h=0) # horizontal line marks the true value of the parameter 
boxplot(results[2,],results[4,]); abline(h=1)

set.seed(1)
results <- replicate(R, {
  x1 = rnorm(N)
  x2 = rexp(N)
  xm = data.frame(x1, x2)
  eps <- rbinom(N, 1, 0.005)
  xm <- (1-eps)*xm + eps*cbind(rep(5, N), rep(15, N))
  smean = colMeans(xm)
  MOM = depthMOM(xm, m = 22)
  c(smean, MOM)
})
# MOM works better
boxplot(results[1,],results[3,]); abline(h=0) # horizontal line marks the true value of the parameter 
boxplot(results[2,],results[4,]); abline(h=1)

set.seed(1)
results <- replicate(R, {
  x1 = rnorm(N)
  x2 = rexp(N)
  xm = data.frame(x1, x2)
  eps <- rbinom(N, 1, 0.005)
  xm <- (1-eps)*xm + eps*cbind(rep(5, N), rep(15, N))
  smean = colMeans(xm)
  MOM = depthMOM(xm, m = 44)
  c(smean, MOM)
})
# MOM works better
boxplot(results[1,],results[3,]); abline(h=0) # horizontal line marks the true value of the parameter 
boxplot(results[2,],results[4,]); abline(h=1)