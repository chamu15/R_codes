library(MASS)
library(mvtnorm)
library(randomForest)
library(robustbase)
library(rrcov)

#------------------------------------------------------------------------------
# Exercise 1
bank <- read.csv(file.choose(), header = T)
head(bank)
attach(bank)
lda.bank <- lda(Y ~ Length+Left+Right+Bottom+Top+Diagonal, data = bank)
lda.bank
predict(lda.bank)$class
sum(predict(lda.bank)$class == bank$Y)
sum(predict(lda.bank)$class != bank$Y)
# one observation is miscalassified

# just out of curiosity, we can try quadratic discriminant analysis
qda.bank <- qda(Y ~ Length+Left+Right+Bottom+Top+Diagonal, data = bank)
qda.bank
predict(qda.bank)$class
sum(predict(qda.bank)$class == bank$Y)
sum(predict(qda.bank)$class != bank$Y) # same result, actually not surprising


#------------------------------------------------------------------------------
# Exercise 2
R <- 100 # number of replications
n0 <- 100
n1 <- 100
Sgm <- diag(1, 2)
mu1 <- rep(0, 2)
mu2 <- rep(2, 2)

# first, let's see what happens without outliers
set.seed(123)
results.lda <- replicate(R, {
  gr.t <- cbind(1,rmvnorm(n1, sigma = Sgm, mean = mu1)) # first mult. normal, label 1
  gr.n <- cbind(0,rmvnorm(n0, sigma = Sgm, mean = mu2)) # second mult. normal, label 0
  dtst <- data.frame(rbind(gr.t,gr.n)) # merge them
  sm.lda <- lda(X1 ~ ., data = dtst)
  lda00 <- sum(dtst[,1]==0 & predict(sm.lda)$class==0) # correct "no subscription"
  lda11 <- sum(dtst[,1]==1 & predict(sm.lda)$class==1) # correct "subscription"
  lda01 <- sum(dtst[,1]==0 & predict(sm.lda)$class==1) # error type I
  lda10 <- sum(dtst[,1]==1 & predict(sm.lda)$class==0) # error type II
  c(lda00, lda11, lda10, lda01)
})
rowMeans(results.lda)

set.seed(123) # to make sure that the same data is generated, to compare methods on the same samples
results.logit <- replicate(R, { # 
  gr.t <- cbind(1,rmvnorm(n1, sigma = Sgm, mean = mu1)) # first mult. normal, label 1
  gr.n <- cbind(0,rmvnorm(n0, sigma = Sgm, mean = mu2)) # second mult. normal, label 0
  dtst <- data.frame(rbind(gr.t,gr.n)) # merge them
  sm.glm <- glm(X1 ~ ., data = dtst, family=binomial(link = logit))
  logit00 <- sum(dtst[,1]==0 & fitted(sm.glm)<0.5)
  logit11 <- sum(dtst[,1]==1 & fitted(sm.glm)>0.5)
  logit01 <- sum(dtst[,1]==0 & fitted(sm.glm)>0.5)
  logit10 <- sum(dtst[,1]==1 & fitted(sm.glm)<0.5)
  c(logit00, logit11, logit10, logit01)
})
rowMeans(results.logit)

set.seed(123)
results.glmrob <- replicate(R, { # 
  gr.t <- cbind(1,rmvnorm(n1, sigma = Sgm, mean = mu1)) # first mult. normal, label 1
  gr.n <- cbind(0,rmvnorm(n0, sigma = Sgm, mean = mu2)) # second mult. normal, label 0
  dtst <- data.frame(rbind(gr.t,gr.n)) # merge them
  sm.glm <- glmrob(X1 ~ ., data = dtst, family=binomial(link = logit), weights.on.x = "hat")
  logit00 <- sum(dtst[,1]==0 & fitted(sm.glm)<0.5)
  logit11 <- sum(dtst[,1]==1 & fitted(sm.glm)>0.5)
  logit01 <- sum(dtst[,1]==0 & fitted(sm.glm)>0.5)
  logit10 <- sum(dtst[,1]==1 & fitted(sm.glm)<0.5)
  c(logit00, logit11, logit10, logit01)
})
rowMeans(results.glmrob)

set.seed(123)
results.ldarob <- replicate(R, { # 
  gr.t <- cbind(1,rmvnorm(n1, sigma = Sgm, mean = mu1)) # first mult. normal, label 1
  gr.n <- cbind(0,rmvnorm(n0, sigma = Sgm, mean = mu2)) # second mult. normal, label 0
  dtst <- data.frame(rbind(gr.t,gr.n)) # merge them
  sm.rlda <- Linda(X1 ~ ., data = dtst)
  rlda00 <- sum(dtst[,1]==0 & predict(sm.rlda)@classification==0)
  rlda11 <- sum(dtst[,1]==1 & predict(sm.rlda)@classification==1)
  rlda01 <- sum(dtst[,1]==0 & predict(sm.rlda)@classification==1)
  rlda10 <- sum(dtst[,1]==1 & predict(sm.rlda)@classification==0)
  c(rlda00, rlda11, rlda10, rlda01)
})
rowMeans(results.ldarob)

nms = c("lda", "ldarob", "glm", "glmrob")
par(mfrow = c(2,2))
boxplot(results.lda[1,], results.ldarob[1,], results.logit[1,], results.glmrob[1,], names = nms, main = "p00")
boxplot(results.lda[3,], results.ldarob[3,], results.logit[3,], results.glmrob[3,], names = nms, main = "p01")
boxplot(results.lda[4,], results.ldarob[4,], results.logit[4,], results.glmrob[4,], names = nms, main = "p10")
boxplot(results.lda[2,], results.ldarob[2,], results.logit[2,], results.glmrob[2,], names = nms, main = "p11")

#---------------------------------------------------------
# with outliers
set.seed(123)
results.lda <- replicate(R, {
  eps <- rbinom(n1, 1, 0.01)
  gr.t <- (1-eps)*cbind(1,rmvnorm(n1, sigma = Sgm, mean = mu1)) + eps*cbind(rep(1, n1), rep(3, n1), rep(3, n1))
  eps <- rbinom(n0, 1, 0.01)
  gr.n <- (1-eps)*cbind(0,rmvnorm(n0, sigma = Sgm, mean = mu2)) + eps*cbind(rep(0, n0), rep(-1, n0), rep(-1, n0))
  dtst <- data.frame(rbind(gr.t,gr.n)) # merge them
  sm.lda <- lda(X1 ~ ., data = dtst)
  lda00 <- sum(dtst[,1]==0 & predict(sm.lda)$class==0) # correct "no subscription"
  lda11 <- sum(dtst[,1]==1 & predict(sm.lda)$class==1) # correct "subscription"
  lda01 <- sum(dtst[,1]==0 & predict(sm.lda)$class==1) # error type I
  lda10 <- sum(dtst[,1]==1 & predict(sm.lda)$class==0) # error type II
  c(lda00, lda11, lda10, lda01)
})
rowMeans(results.lda)

set.seed(123) # 
results.logit <- replicate(R, { # 
  eps <- rbinom(n1, 1, 0.01)
  gr.t <- (1-eps)*cbind(1,rmvnorm(n1, sigma = Sgm, mean = mu1)) + eps*cbind(rep(1, n1), rep(3, n1), rep(3, n1))
  eps <- rbinom(n0, 1, 0.01)
  gr.n <- (1-eps)*cbind(0,rmvnorm(n0, sigma = Sgm, mean = mu2)) + eps*cbind(rep(0, n0), rep(-1, n0), rep(-1, n0))
  dtst <- data.frame(rbind(gr.t,gr.n)) # merge them
  sm.glm <- glm(X1 ~ ., data = dtst, family=binomial(link = logit))
  logit00 <- sum(dtst[,1]==0 & fitted(sm.glm)<0.5)
  logit11 <- sum(dtst[,1]==1 & fitted(sm.glm)>0.5)
  logit01 <- sum(dtst[,1]==0 & fitted(sm.glm)>0.5)
  logit10 <- sum(dtst[,1]==1 & fitted(sm.glm)<0.5)
  c(logit00, logit11, logit10, logit01)
})
rowMeans(results.logit)

set.seed(123)
results.glmrob <- replicate(R, { # 
  eps <- rbinom(n1, 1, 0.01)
  gr.t <- (1-eps)*cbind(1,rmvnorm(n1, sigma = Sgm, mean = mu1)) + eps*cbind(rep(1, n1), rep(3, n1), rep(3, n1))
  eps <- rbinom(n0, 1, 0.01)
  gr.n <- (1-eps)*cbind(0,rmvnorm(n0, sigma = Sgm, mean = mu2)) + eps*cbind(rep(0, n0), rep(-1, n0), rep(-1, n0))
  dtst <- data.frame(rbind(gr.t,gr.n)) # merge them
  sm.glm <- glmrob(X1 ~ ., data = dtst, family=binomial(link = logit), weights.on.x = "hat")
  logit00 <- sum(dtst[,1]==0 & fitted(sm.glm)<0.5)
  logit11 <- sum(dtst[,1]==1 & fitted(sm.glm)>0.5)
  logit01 <- sum(dtst[,1]==0 & fitted(sm.glm)>0.5)
  logit10 <- sum(dtst[,1]==1 & fitted(sm.glm)<0.5)
  c(logit00, logit11, logit10, logit01)
})
rowMeans(results.glmrob)

set.seed(123)
results.ldarob <- replicate(R, { # 
  eps <- rbinom(n1, 1, 0.01)
  gr.t <- (1-eps)*cbind(1,rmvnorm(n1, sigma = Sgm, mean = mu1)) + eps*cbind(rep(1, n1), rep(3, n1), rep(3, n1))
  eps <- rbinom(n0, 1, 0.01)
  gr.n <- (1-eps)*cbind(0,rmvnorm(n0, sigma = Sgm, mean = mu2)) + eps*cbind(rep(0, n0), rep(-1, n0), rep(-1, n0))
  dtst <- data.frame(rbind(gr.t,gr.n)) # merge them
  sm.rlda <- Linda(X1 ~ ., data = dtst)
  rlda00 <- sum(dtst[,1]==0 & predict(sm.rlda)@classification==0)
  rlda11 <- sum(dtst[,1]==1 & predict(sm.rlda)@classification==1)
  rlda01 <- sum(dtst[,1]==0 & predict(sm.rlda)@classification==1)
  rlda10 <- sum(dtst[,1]==1 & predict(sm.rlda)@classification==0)
  c(rlda00, rlda11, rlda10, rlda01)
})
rowMeans(results.ldarob)

nms = c("lda", "ldarob", "glm", "glmrob")
par(mfrow = c(2,2))
boxplot(results.lda[1,], results.ldarob[1,], results.logit[1,], results.glmrob[1,], names = nms, main = "p00")
boxplot(results.lda[3,], results.ldarob[3,], results.logit[3,], results.glmrob[3,], names = nms, main = "p01")
boxplot(results.lda[4,], results.ldarob[4,], results.logit[4,], results.glmrob[4,], names = nms, main = "p10")
boxplot(results.lda[2,], results.ldarob[2,], results.logit[2,], results.glmrob[2,], names = nms, main = "p11")

# Although, the performance of methods was reduced a bit, surprisingly, the difference is not dramatic. 
# But we will not make huge claims, that outliers do not matter.
# There are more things to try: sample sizes, contaminations by more heavy outliers, out-of-sample classification,
# the dimentionality, i.e. the number of predictors, etc..