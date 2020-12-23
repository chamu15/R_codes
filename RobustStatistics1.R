## asymptotic variances at t-distribution

ASVmean <- function(df) df / (df - 2)

ASVmedian <- function(df) {
  M <- qt(0.5, df = df)
  1 / (4 * dt(M, df = df)^2)
}

ASVsd <- function(df) {
  mu4 <- 1 / (sqrt(pi) * gamma(df/2)) * (gamma(5/2) * gamma((df - 4)/2) * df^2)
  sigma2 <- df / (df - 2)
  1 / (4 * sigma2) * (mu4 - sigma2^2)
}

ASVmad <- function(df) {
  Q <- qt(0.75, df = df)
  1 / (16 * Q^2 * dt(Q, df = df)^2)
    
}

# 5 degrees of freedom
df <- 5
ASVmean(df) / ASVmedian(df)
ASVsd(df) / ASVmad(df)

# 10 degrees of freedom
df <- 10
ASVmean(df) / ASVmedian(df)
ASVsd(df) / ASVmad(df)

# 3 degrees of freedom
df <- 3
ASVmean(df) / ASVmedian(df)
# asymptotic variance of standard deviation cannot be computed as 4th central
# moment doesn't exist
ASVmad(df)


## correlation measures

library("mvtnorm")

# control parameters
# rho <- 0.5
rho <- 0.9
Sigma <- matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2)
n <- 100
R <- 500

# run simulation
set.seed(1)
results <- replicate(R, {
  xy <- rmvnorm(n, sigma = Sigma)
  x <- xy[, 1]
  y <- xy[, 2]
  c(Pearson = cor(x, y, method = "pearson"),
    Spearman = cor(x, y, method = "spearman"),
    Kendall = cor(x, y, method = "kendall"))
})

# compute average simulation results
avg <- rowMeans(results)
# compare to true correlation rho
abs(avg - rho)
# compare to true quantities
true <- c(Pearson = rho, Spearman = 6/pi*asin(rho/2), Kendall = 2/pi*asin(rho))
abs(avg - true)

# Fisher-consistency correction
corPearson <- function(...) cor(..., method = "pearson")
corSpearman <- function(...) 2 * sin(pi/6 * cor(..., method = "spearman"))
corKendall <- function(...) sin(pi/2 * cor(..., method = "kendall"))

# control parameters
rho <- 0.9
Sigma <- matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2)
n <- 100
R <- 500
epsilon <- seq(0, 0.1, by = 0.01)

# run simulation
resultsList <- lapply(epsilon, function(eps) {
  # reset random seed such that clean data is
  # the same for each contamination level
  set.seed(20170227)
  resultsEps <- replicate(R, {
    # generate clean data
    xy <- rmvnorm(n, sigma = Sigma)
    # replace the first observations
    nout <- n * eps
    x <- c(rep.int(-5, nout), xy[-(1:nout), 1])
    y <- c(rep.int(5, nout), xy[-(1:nout), 2])
    # compute correlation estimates
    c(Pearson = corPearson(x, y),
      Spearman = corSpearman(x, y),
      Kendall = corKendall(x, y))
  })
  data.frame(Epsilon = eps, t(resultsEps))
})
results <- do.call(rbind, resultsList)

# create plot of results
library("ggplot2")
library("tidyr")
library("plyr")
results <- gather(results, key = Method, value = Estimate, -Epsilon)
aggregated <- ddply(results, .(Epsilon, Method), summarize,
                    Average = mean(Estimate))
ggplot(aggregated, aes(x = Epsilon, y = Average, color = Method)) +
  geom_line() + geom_hline(yintercept = rho, color = "darkgrey")
