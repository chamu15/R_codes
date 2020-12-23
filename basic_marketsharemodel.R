# Advanced Marketing Models, assignment week 2
# Thanks to: Peter van den Brand

library('Matrix')

rm(list = ls())

ols <- function(y, X, n) {
  # Compute OLS betas
  beta <- solve(t(X) %*% X, t(X) %*% y)
  
  # Estimate covariance of errors
  e <- y - X %*% beta
  e <- matrix(e, ncol=n)
  S <- t(e) %*% e / nrow(e)
  
  return(list(beta=beta, S=S))
}

fgls <- function(y, X, n, eps=1e-6) {
  start <- ols(y, X, n)
  S <- start$S
  t <- nrow(X) / n
  prev <- NULL
  
  while (TRUE) {
    omega <- kronecker(solve(S), diag(t))  
    beta <- solve(t(X) %*% omega %*% X, t(X) %*% omega %*% y)
    
    e <- y - X %*% beta
    e <- matrix(e, ncol=n)
    S <- t(e) %*% e / nrow(e)
    
    sse <- sum(e ^ 2)
    if (!is.null(prev) && abs(sse - prev) < eps) {
      break
    }
    prev <- sse
  }
  
  return(list(beta=beta, S=S))
}

#------------------------------------------------------------------------------
# Read data and compute market shares

data <- read.csv('beerdata1.txt', sep='\t')
data <- as.matrix(subset(data, select=-c(OBS)))   #remove the OBS column

train <- data[1:201,]
test <- data[202:227,]

volumes <- train[,c('VOL1', 'VOL2', 'VOL3', 'VOL4')]
M <- volumes / apply(volumes, 1, sum)

volumes <- test[,c('VOL1', 'VOL2', 'VOL3', 'VOL4')]
M.test <- volumes / apply(volumes, 1, sum)

#------------------------------------------------------------------------------
# FE-MCI model

y <- c(log(M[,'VOL1']) - log(M[,'VOL4']), 
       log(M[,'VOL2']) - log(M[,'VOL4']), 
       log(M[,'VOL3']) - log(M[,'VOL4']))

X.fe <- cbind(
  1,  # intercept
  log(train[,'PRICE1']),
  log(train[,'PRICE2']),
  log(train[,'PRICE3']),
  log(train[,'PRICE4']),
  train[,c('PROMO11', 'PROMO12', 'PROMO13', 'PROMO14', 'PROMO21', 'PROMO22', 'PROMO23', 'PROMO24')])
X.fe <- as.matrix(bdiag(X.fe, X.fe, X.fe))  # create block diagonal matrix

# Use OLS to estimate the FE-MCI model
fe_mci <- ols(y, X.fe, 3)

#------------------------------------------------------------------------------
# RC model

X1 <- cbind(1, log(train[,'PRICE1']), train[,'PROMO11'], train[,'PROMO21'])
X2 <- cbind(1, log(train[,'PRICE2']), train[,'PROMO12'], train[,'PROMO22'])
X3 <- cbind(1, log(train[,'PRICE3']), train[,'PROMO13'], train[,'PROMO23'])

X4 <- cbind(log(train[,'PRICE4']), train[,'PROMO14'], train[,'PROMO24'])

X.rc <- as.matrix(bdiag(X1, X2, X3))
X.rc <- cbind(X.rc, rbind(X4, X4, X4))

# Estimate RC model with FGLS
rc_mci <- fgls(y, X.rc, 3)

# Perform LR test
lr <- nrow(train) * (log(det(rc_mci$S)) - log(det(fe_mci$S)))
p <- 1 - pchisq(lr, 3 * 4 * (4 - 2))
print(p)
# P-value of LR test (p) is 0, which means that we should use the unrestricted FE-MCI model.
# An RE model will also be rejected in favor of FE-MCI, because it's more restricted than RC.
# I haven't tested the RCM model, because it seems like it takes a lot of time to implement (and it only restricts 2 parameters).

#------------------------------------------------------------------------------
# Forecasting using the FE-MCI model

X.test <- cbind(
  1,  # intercept
  log(test[,'PRICE1']),
  log(test[,'PRICE2']),
  log(test[,'PRICE3']),
  log(test[,'PRICE4']),
  test[,c('PROMO11', 'PROMO12', 'PROMO13', 'PROMO14', 'PROMO21', 'PROMO22', 'PROMO23', 'PROMO24')])
X.test <- as.matrix(bdiag(X.test, X.test, X.test))  # create block diagonal matrix

yhat <- X.test %*% fe_mci$beta
yhat <- matrix(yhat, ncol=3)

# Naive forecast
half_sigma <- 0.5 * diag(fe_mci$S)
m_it <- exp(yhat + matrix(rep(half_sigma, nrow(yhat)), ncol=3, byrow=TRUE))
denom <- apply(m_it, 1, sum) + 1
M.naive <- m_it / denom
M.naive <- cbind(M.naive, 1 / denom)

RMSE.naive <- sqrt(colMeans((M.test - M.naive) ^ 2))
R2.naive <- sapply(1:4, function(i) { cor(M.test[,i], M.naive[,i])^2 })

# Simulation forecast
set.seed(12345)
P <- t(chol(fe_mci$S))
M.sim <- matrix(0, nrow=nrow(yhat), ncol=4)
n <- 10000
for (i in 1:n) {
  z <- t(P %*% matrix(rnorm(3 * nrow(yhat)), nrow=3))
  m_it <- exp(yhat + z)

  denom <- apply(m_it, 1, sum) + 1
  s <- m_it / denom
  s <- cbind(s, 1 / denom)
  M.sim <- M.sim + s
}
M.sim <- M.sim / n

RMSE.sim <- sqrt(colMeans((M.test - M.sim) ^ 2))
R2.sim <- sapply(1:4, function(i) { cor(M.test[,i], M.sim[,i])^2 })

print(RMSE.sim - RMSE.naive)

# Make a plot of all predictions and realizations
par(mfrow=c(2,2))
for (b in 1:4) {
  share = cbind(M.test[,b],M.naive[,b],M.sim[,b])
  matplot(share,type="b",pch=1)
  if(b==1){
    legend("topleft", legend = c("true","naive","sim"), col=1:3, pch=1)
  }
}
# Conclusion: simulation gives a slightly lower RMSE (more accurate forecast) for 3 out of 4 products in this category
