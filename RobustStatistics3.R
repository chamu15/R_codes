# load packages
#install.packages("ssmrob")
library("ggplot2")
library("robustbase")

# Boston housing data
load(file.choose())

# formula for model in the literature
f <- logMedValue ~ CrimeRate + Lots + Industrial + CharlesRiver + NOxSq +
  RoomsSq + Pre1940 + logDistance + logHighways + Tax +
  PTRatio + Black + logStatus + Longitude + Latitude +
  LongitudeSq + LatitudeSq + LongitudeLatitude

# OLS
fitOLS <- lm(f, data = Boston)
summary(fitOLS)

# MM-estimator
fitMM <- lmrob(f, data = Boston)                   # gives warning
fitMM <- lmrob(f, data = Boston, cov = ".vcov.w")  # set argument accordingly
summary(fitMM)  # argument cov
# The parameter estimates by two methods are different by substantial margins. 
# As a rule of thumb consider a difference of one std.error.
# The significance of the parameters is different.

## forward or backward selection?

# The more variables, the bigger the probability of outliers.  Hence forward
# selection is preferable from a robustness point of view.


## robust forward stepwise variable selection via BIC

# lmrob() without overhead from formula interface
# also set some control arguments for the MM-algorithm to prevent warnings
fastLmrob <- function(x, y, max.it = 200, k.max = 1000, cov = ".vcov.w", ...) {
  control <- lmrob.control(max.it = max.it, k.max = k.max, cov = cov, ...)
  fit <- lmrob.fit(x, y, control = control)
  class(fit) <- "fastLmrob"
  fit
}

# BIC for results of fastLmrob()
BIC.fastLmrob <- function(object, ...) {
  n <- length(object$rweights)
  n * log(object$scale^2) + object$rank * log(n)
}

# set up response and candidate predictors
y <- Boston[, 1]
X <- as.matrix(Boston[, -1])
n <- length(y)
p <- ncol(X)

# start with constant term
x <- cbind("(Intercept)" = rep.int(1, n))
fit <- fastLmrob(x, y)
bic <- c(BIC(fit), rep.int(NA_real_, p))

# loop over candidate predictors
candidates <- colnames(X)
selected <- rep.int(NA_character_, p)
for (j in seq_len(p)) {
  # find not yet active candidate that yields smallest BIC
  values <- sapply(candidates, function(v) {
    x <- cbind(x, X[, v, drop = FALSE])  # add current variable
    fit <- fastLmrob(x, y, max.it = 200, k.max = 1000, cov = ".vcov.w")  # compute MM-estimator
    BIC(fit)                             # return BIC
  })
  best <- which.min(values)
  # add that candidate to explanatory variables
  selected[j] <- candidates[best]
  bic[j+1] <- values[best]
  x <- cbind(x, X[, selected[j], drop = FALSE])
  candidates <- candidates[-best]
}

# find optimal step
sOpt <- which.min(bic) - 1  # step count starts with 0

# plot BIC values
dfBIC <- data.frame(Step = 0:p, BIC = bic)
ggplot() +
  geom_vline(aes(xintercept = sOpt), color = "darkblue") +
  geom_line(aes(x = Step, y = BIC), dfBIC)

# construct formula for final model
lhs <- names(Boston)[1]
variables <- head(selected, sOpt)
rhs <- paste(variables, collapse = " + ")
f <- as.formula(paste(lhs, rhs, sep = "~"))

# estimate final forward stepwise model
fitFS <- lmrob(f, data = Boston, cov = ".vcov.w")
summary(fitFS)

# Note that standard errors are greatly overestimated as the explanatory
# variables are chosen that fit the given data best, hence significance is
# exagerrated.  Accurate standard errors for variable selection procedures
# with robust methods is still an unsolved problem in the literature.

# The only method at the moment that might be suitable is a procedure based on
# data splitting developed by Wasserman and Roeder (2009, Annals of Statistics).
# But this procedure could yield problems if outliers are oversampled in one of
# the subsets from data splitting.


#------------------------------------------------------------------------------
# Exercise 2
library("mvtnorm")
library("robustbase")
# control parameters
n <- 1000
R <- 500
bt <- c(1, 0.5, 0.25)
sgmX <- matrix(c(1,0.5, 0.25, 0.5, 1, 0.5, 0.25, 0.5, 1),3,3)
mu <- rep(0, 3)

# a 
results <- replicate(R, {
  # generate complete data
  x <- rmvnorm(n, mean = mu, sigma = sgmX)
  LP <- x%*%bt 
  y <- rbinom(n, 1, 1/(1+exp(-LP)))
  fit.ml = coef(glm(y~x, family = binomial(link = logit)))
  fit.rob = coef(glmrob(y~x, family = binomial(link = logit), weights.on.x = "covMcd"))
  c(fit.ml, fit.rob)
})

rowMeans(results)
par(mfrow = c(1, 4))
boxplot(results[1,], results[5,]); abline(h = 0)
boxplot(results[2,], results[6,]); abline(h = bt[1])
boxplot(results[3,], results[7,]); abline(h = bt[2])
boxplot(results[4,], results[8,]); abline(h = bt[3])

# b 
results <- replicate(R, {
  # generate complete data
  x <- rmvnorm(n, mean = mu, sigma = sgmX)
  LP <- x%*%bt 
  y <- rbinom(n, 1, 1/(1+exp(-LP)))
  eps <- rbinom(n, 1, 0.01) # contamination
  y <- (1-eps)*y + eps*1 # either original y or 1
  x <- (1-eps)*x + eps*matrix(rep(-2, n*3), n, 3)
  fit.ml = coef(glm(y~x, family = binomial(link = logit)))
  fit.rob = coef(glmrob(y~x, family = binomial(link = logit), weights.on.x = "covMcd"))
  c(fit.ml, fit.rob)
})

rowMeans(results)
par(mfrow = c(1, 4))
boxplot(results[1,], results[5,]); abline(h = 0)
boxplot(results[2,], results[6,]); abline(h = bt[1])
boxplot(results[3,], results[7,]); abline(h = bt[2])
boxplot(results[4,], results[8,]); abline(h = bt[3])
# MLE is clearly biased, Robust is close to the true value

#------------------------------------------------------------------------------
# Exercise 3
library("ssmrob")
data("MEPS2001")
modeleq <- dambexp ~ age + female + educ + blhisp + totchr + ins

fit.clas <- glm(modeleq, family = binomial(link = logit), data = MEPS2001)
summary(fit.clas)
fit.rob <- glmrob(modeleq, family = binomial(link = logit), weights.on.x = "covMcd", data = MEPS2001) #  Error: numerical issue
fit.rob <- glmrob(modeleq, family = binomial(link = logit), weights.on.x = "hat", data = MEPS2001) # now it is fine
summary(fit.rob)
# the fits are very close to each other. The robust estimator is slightly less efficient than the MLE,
# what corresponds to the theory.
# One significance level is different, but the difference is small and due to the higher variablity of 
# the robust estimator.
