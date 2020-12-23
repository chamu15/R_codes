## correlation measures
install.packages("robustbase")
library("mvtnorm")
library("robustbase")

# Exercise 1
# control parameters
n <- 20
# n <- 30
# n <- 50
p <- 10
rho <- 0.2
Sigma <- rho^t(sapply(1:p, function(i, j) abs(i-j), 1:p))
R <- 1000

# Fisher-consistency correction
corPearson <- function(...) cor(..., method = "pearson")
corSpearman <- function(...) 2 * sin(pi/6 * cor(..., method = "spearman"))
corKendall <- function(...) sin(pi/2 * cor(..., method = "kendall"))

# run simulation
set.seed(1)
results <- replicate(R, {
  # generate complete data
  x <- rmvnorm(n, sigma = Sigma)
  scales <- apply(x, 2, mad, na.rm = TRUE)
  SS <- corSpearman(x) * tcrossprod(scales)
  SK <- corKendall(x) * tcrossprod(scales)
  c(Spearman = any(eigen(SS)$values < 0), Kendall = any(eigen(SK)$values < 0))
})
rowSums(results)
# The resulting covariance matrices are not always positive (semi)definite. 
# If something like this happens in analysis of real data, then this contradicts to the geometric intuition,
# that the volume of the cloud of points is positive.
# It can also lead to intractable Mahalanobis distances.

#------------------------------------------------------------------------------
# Exercise 2
# control parameters
N <- 100
p <- 2
mu <- rep(0,p)
sgm <- matrix(c(1,0.5,0.5,1), 2, 2)
R <- 500

# subquestion a
set.seed(123)
results <- replicate(R, {
  # generate complete data
  x <- rmvnorm(N, mean = mu, sigma = sgm)
  sm <- colMeans(x)
  sc <- cov(x)
  mcd <- covMcd(x)
  c(sc, mcd$cov, sm, mcd$center)
})
rowMeans(results) # The output is a vector 12*1
# Average sample variances are 1st and 4th, average sample covariance is 2nd or 3rd.
# Average MCD estimated variances are 5th and 8th, average MCD estimated covariance is 6th or 7th.
# Sample means are 9th and 10th
# MCD estimates of expectation are 11th and 12th
# On average we estimate the true parameters.

# subquestion b
set.seed(123) 
results <- replicate(R, {
  # generate complete data
  x <- rmvnorm(N, mean = mu, sigma = sgm)
  eps <- rbinom(N, 1, 0.05) # contamination with probability 0.05
  xc <- (1-eps)*x + eps*cbind(rep(2, N), rep(-2,N)) # either from original data or outlier
  sm <- colMeans(xc)
  sc <- cov(xc)
  mcd <- covMcd(xc)
  c(sc, mcd$cov, sm, mcd$center)
})
rowMeans(results)
# The sample covariance is biased, the MCD covariance is closer to the true value
# The variances by both estimators are slightly affected. Let us see what happens with heavier contamination

# subquestion c
set.seed(123) 
results <- replicate(R, {
  # generate complete data
  x <- rmvnorm(N, mean = mu, sigma = sgm)
  eps <- rbinom(N, 1, 0.05)
  xc <- (1-eps)*x + eps*cbind(rep(5, N), rep(-5,N))
  sm <- colMeans(xc)
  sc <- cov(xc)
  mcd <- covMcd(xc)
  c(sc, mcd$cov, sm, mcd$center)
})
rowMeans(results)
# The MCD estimates are stable. The sample covariance matrix is completely broken down.
# Even the correlation has changed the sign.


