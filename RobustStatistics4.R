# Exercise 1
# a) 
# Hotelling's T^2
banknote = read.csv(file.choose(), header = T)
nominator = colMeans(banknote[101:200,1:6])-c(215,130,130,10,10,142)
100*nominator%*%solve(var(banknote[101:200,1:6]))%*%nominator
qf(0.95, 6, 94)*6*99/94 # critical value, where p=6, n=100
# reject H_0

# Using likelihood ratio test
n <- 100
u <- rep(1, n)
X <- banknote[101:200,1:6]
mu_null <- c(215,130,130,10,10,142)
cov.null <- as.matrix(t(X-u%*%t(mu_null)))%*%as.matrix((X-u%*%t(mu_null)))/n
cov.alt <- as.matrix(t(X-u%*%t(colMeans(X))))%*%as.matrix((X-u%*%t(colMeans(X))))/n
-n*log(det(cov.alt)/det(cov.null)) # LRT, i.e. the test statistic
qchisq(0.95, 6) # critical value
# also reject

# b)
g.bank <- banknote[1:100, c(6,1,2,3,4,5)] # cheap trick, change the order in the data matrix,
# to use the formula from the slide 23.
est.cor <- cor(g.bank)
R22 <- est.cor[2:6, 2:6]
R12 <- est.cor[1, 2:6]
-2*log((1 - R12%*%solve(R22)%*%(R12))^(100/2))
qchisq(0.95, 5) # reject H_0

# can we get the same result if we skip some algebra simplifications,
# which are made to get the result in slide 23?
# we write two likelihoods
n=100
X <- banknote[1:n, c(6,1,2,3,4,5)]

cov.null <- matrix(rep(NA, 6*6), 6, 6)
cov.null[1, 1] <- (n-1)/n*var(X[, 1])
cov.null[1,2:6] <- rep(0, 5)
cov.null[2:6,1] <- cov.null[1,2:6]
cov.null[2:6,2:6] <- (n-1)/n*var(X)[2:6,2:6]
cov.alt <- (n-1)/n*var(X)
tr <- function(matr) {sum(diag(matr))} 
nominator <- det(cov.null)^(-n/2)*exp(-1/2*n*tr(solve(cov.null)%*%cov.alt)) # This is L*_0, 
# I could have written the entire sum (x-bar(x))'\Sigma_0^{-1}(x-bar(x)), but using the trace is easier
denominator <- det(cov.alt)^(-n/2)*exp(-1/2*n*6) # This is L*_A
-2*log(nominator/denominator) # the same value of the test statistic

#------------------------------------------------------------------------
# Exercise 2
# In principle, we can implement the Holm procedure and Bonferroni correction ourselves,
# but I'll use a standard function p.adjust from the "stats" package (is loaded automatically).
# It adjusts the p-values such that they can be compared to the unadjusted level \alpha.

# how many false "discoveries" do we make?
R <- 1000
M <- 100 # number of tested hypotheses
N <- 100 # sample size
set.seed(123) # fix the random number generator for reproducibility

results <- replicate(R, {
  x <- matrix(rnorm(N*M, 0, 3), N, M) # generate data
  cm <- colMeans(x) # vector for \bar{X}'s
  csd <- apply(x, 2, sd) # vector for estimated standard deviations
  ct <- sqrt(N)*cm/csd # vector for test statistics
  pv <- 1 - pnorm(ct) # vector of unadjusted p-values

  r.holm <- sum(p.adjust(pv, "holm") < 0.025) # two-sided alternative
  r.bonf <- sum(p.adjust(pv, "bonferroni") < 0.025)
  r.no <- M - sum(abs(ct) < qnorm(0.975)) # no adjustment for multiple testing
  c(r.holm > 0, r.bonf > 0, r.no > 0) # report whether we made at least one false discovery
}
)
rowMeans(results) # report the average number of cases when we made at least one false discovery
# Less than 5% for Holm and Bonferroni,
# very close to 1 without adjustment
#------------------------------------------------------------------------
# what about power? 
R <- 1000
M <- 100 # number of tested hypotheses
N <- 100 # sample size
set.seed(123) # fix the random number generator for reproducibility

results <- replicate(R, {
  x <- matrix(rnorm(N*M, 1, 3), N, M) # generate data, notice that the mean is 1
  # you can check what happens with the power of the procedures for different values of the mean
  # and different sample sizes.
  cm = colMeans(x) # vector for \bar{X}'s
  csd = apply(x, 2, sd) # vector for estimated standard deviations
  ct = sqrt(N)*cm/csd # vector for test statistics, H_0 is that \mu=0
  pv = 1 - pnorm(ct) # vector of unadjusted p-values
  
  r.holm <- sum(p.adjust(pv, "holm") < 0.025) # two-sided alternative
  r.bonf <- sum(p.adjust(pv, "bonferroni") < 0.025)
  r.no <- M - sum(abs(ct)<1.96) # no adjustment
  c(r.holm/M, r.bonf/M, r.no/M) # percentage of correctly rejected false nulls
}
)
rowMeans(results) # Holm's procedure rejects false nulls more often than Bonferroni
