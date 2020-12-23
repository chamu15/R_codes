# Estimates static and dynamically extended market share models for training set, and produces naive and simulated predictions for both models and plots those together and against each other.
# Thanks to Job Overbeek

rm(list=ls())
require(ggplot2)
require(magic)
require(fBasics) #install(fBasics)

repmat = function(X,m,n){
  ##R equivalent of repmat (matlab)
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)}

surReg <- function(X, y, k){
  #  k = number of brands estimated ( I )
  N <- NROW(X)
  time <- N / (k-1) # time = number of time periods ( T )
  
  kronSigmaInv <- kronecker(diag(k-1), diag(time))  # initial sigma estimate
  sseold <- -Inf
  ssenew <- 0
  iter <- 0
  while (abs(ssenew - sseold) > 0.000001){
    bSUR <- solve(t(X) %*% kronSigmaInv %*% X) %*% t(X) %*% kronSigmaInv %*% y
    
    e_i <- y - X %*% bSUR
    E_i <- matrix(e_i, ncol = (k-1))  #[ e_1; e_2; ... e_I-1] -> [ e_1 e_2 ... e_I-1]  -> e_t are now on the rows of the matrix
    E_t <- t(E_i)  # e_t on the columns
    sigmaHat <- E_t %*% t(E_t) / NCOL(E_t)
    
    # Update sigma
    sigma <- sigmaHat
    kronSigmaInv <- kronecker(solve(sigma), diag(time))
    
    #update sse
    sseold <- ssenew
    ssenew <- sum(e_i^2)
    iter <- iter + 1
  }
  cov_b <- solve(t(X) %*% kronSigmaInv %*% X)
  tval <- bSUR / sqrt(diag(cov_b))
  pval <- 2*pt(q = abs(tval), df = N, lower.tail = FALSE)
  return(list(beta=bSUR, detSigmaHat=sigmaHat, iterations=iter, tvalue=tval, pvalue=pval))
}

naiveMarketShareForecast <- function(Xpred, bHat, sigmaHat, k){
  #  k = # of brands
  yPred <- Xpred %*% bHat
  lnm_it <- matrix(yPred, ncol = (k-1))
  m_it <- lnm_it
  for(i in 1:(k-1)){
    m_it[, i] <- exp(lnm_it[, i] + 0.5*sigmaHat[i,i])
  }
  marketShare <- cbind(m_it /(1+rowSums(m_it)), 1/(1+rowSums(m_it)))
}

simulationMarketShareForecast <- function(Xpred, bHat, sigmaHat, k, nSim){
  #  k = # of brands
  #  nSim = # of simulation runs
  N <- NROW(Xpred)
  P <- t(chol(sigmaHat))
  
  simulatedMarketShares <- array(0, c(N/(k-1), k, nSim))
  for (i in 1:nSim){

    drawn_errors = vec(P %*% matrix(rnorm(N, 0, 1),nrow=(k-1)))
    yPred <- Xpred %*% bHat + drawn_errors
    
    lnm_it <- matrix(yPred, ncol = (k-1)) 
    m_it <- lnm_it
    for(j in 1:(k-1)){
      m_it[, j] <- exp(lnm_it[, j] + 0.5*sigmaHat[j, j])
    }
    simulatedMarketShares[,,i] <- cbind(m_it /(1+rowSums(m_it)), 1/(1+rowSums(m_it)))
  }
  averageSimulatedMarketShares <- apply(simulatedMarketShares, c(1,2), mean)  # average per time and brand over all simulations
}

#----------------------- Load data and select general model variables
filePath <- "beerdata1.txt"
rawData <- read.csv(filePath, sep='\t') 

# Determine M_i and m_i
S_t  <- rawData[, paste0("VOL", 1:4)] 
M_t <- S_t/rowSums(as.matrix(S_t))
logM_t <- log(M_t)
m_t <- as.matrix(logM_t[, 1:3] - logM_t[, 4])
m <- c(m_t)
logPrices <- log(rawData[, paste0("PRICE", 1:4)]) # log prices
indx.predict <- 206:227
M_t.predict <- M_t[indx.predict, ]


#---------------------------------------------------- FE static and dynamic (1 price lag, 2 market share lags) models
variableNamesPart <- c("constant", paste("logPrice", 1:4, sep=""), paste("logPriceLag1_", 1:4, sep=""))
variableNames <- c(variableNamesPart, "Market1Lag1", "Market1Lag2", variableNamesPart, "Market2Lag1", "Market2Lag2", variableNamesPart, "Market3Lag1", "Market3Lag2", "Market4Lag1", "Market4Lag2")

#  Training set matrices
wi.train <- as.matrix(cbind(rep(1, 203), logPrices[3:205,])) 
w1.train <- as.matrix(cbind(wi.train, logPrices[2:204,], logM_t[2:204, 1], logM_t[1:203, 1])) # dynamic
w2.train <- as.matrix(cbind(wi.train, logPrices[2:204,], logM_t[2:204, 2], logM_t[1:203, 2])) 
w3.train <- as.matrix(cbind(wi.train, logPrices[2:204,], logM_t[2:204, 3], logM_t[1:203, 3]))
W.train <- adiag(w1.train, w2.train, w3.train)  # block diagonal matrix
zi.train <- -cbind(logM_t[2:204, 4], logM_t[1:203, 4])  # common market share parameters

X.train.static <- kronecker(diag(3), wi.train)   # intercept and lnPrice
colnames(X.train.static) <- rep(c("constant", paste("logPrice", 1:4, sep="")),3)
X.train.dynamic <- cbind(W.train, repmat(zi.train, 3, 1))  # intercept, lnPrice, lnPrice_{-1} and M_i (lag 1 and lag 2)
colnames(X.train.dynamic) <- variableNames

y_train <- c(m_t[3:205, ])

#  Prediction set matrices
wi.pred <- as.matrix(cbind(rep(1, 22), logPrices[206:227,])) # non-dynamic
w1.pred <- as.matrix(cbind(wi.pred, logPrices[205:226,], logM_t[205:226, 1], logM_t[204:225, 1])) # dynamic
w2.pred <- as.matrix(cbind(wi.pred, logPrices[205:226,], logM_t[205:226, 2], logM_t[204:225, 2])) 
w3.pred <- as.matrix(cbind(wi.pred, logPrices[205:226,], logM_t[205:226, 3], logM_t[204:225, 3]))
W.pred <- adiag(w1.pred, w2.pred, w3.pred)  
zi.pred <- -cbind(logM_t[205:226, 4], logM_t[204:225, 4]) 

X.pred.static <- kronecker(diag(3), wi.pred)
colnames(X.pred.static) <- rep(c("constant", paste("logPrice", 1:4, sep="")),3)
X.pred.dynamic <- cbind(W.pred, repmat(zi.pred, 3, 1))
colnames(X.train.dynamic) <- variableNames

y_pred <- c(m_t[206:227, ])


#-------------------------------------------------  Estimate static and dymanic models
feEstimatesStatic <- surReg(X = X.train.static, y = y_train, k = 4)
feEstimatesDynamic <- surReg(X = X.train.dynamic, y = y_train, k = 4)


#-------------------------------- Compute naive forecasts
naiveMarketStatic <- naiveMarketShareForecast(Xpred = X.pred.static, bHat = feEstimatesStatic$beta, sigmaHat = feEstimatesStatic$detSigmaHat, k = 4)
naiveMarketDynamic <- naiveMarketShareForecast(Xpred = X.pred.dynamic, bHat = feEstimatesDynamic$beta, sigmaHat = feEstimatesDynamic$detSigmaHat, k = 4)


#-------------------------------- Compute simulated forecasts
simulationMarketStatic <- simulationMarketShareForecast(Xpred = X.pred.static, bHat = feEstimatesStatic$beta, sigmaHat = feEstimatesStatic$detSigmaHat, k = 4, nSim = 2000)
simulationMarketDynamic <- simulationMarketShareForecast(Xpred = X.pred.dynamic, bHat = feEstimatesDynamic$beta, sigmaHat = feEstimatesDynamic$detSigmaHat, k = 4, nSim = 2000)


#----------------------  plot realizations and between-naive predictions
ID <- c(rep("Brand 1", 22), rep("Brand 2", 22), rep("Brand 3", 22), rep("Brand 4", 22))
plotBetweenNaive <- data.frame(id = rep(1:22, 4*3), 
                               value = 100*c(c(as.matrix(naiveMarketStatic)), c(as.matrix(naiveMarketDynamic)), c(as.matrix(M_t.predict))), 
                               predictionType = c(rep("NaiveStatic",88), rep("NaiveDynamic", 88), rep("Realization",88)),
                               share = rep(ID, 3))

ggplot(plotBetweenNaive, aes(x = id, y = value)) + geom_point() + geom_line() + facet_grid( share ~ predictionType) + xlab("time") + ylab("market share in %")
ggplot(plotBetweenNaive, aes(x = id, y = value, colour=predictionType)) + geom_line() + geom_point() +facet_grid( share~. ) + xlab("time") + ylab("market share in %")


#----------------------  plot realizations and between-simulation predictions
ID <- c(rep("Brand 1", 22), rep("Brand 2", 22), rep("Brand 3", 22), rep("Brand 4", 22))
plotBetweenSimulation <- data.frame(id = rep(1:22, 4*3), 
                                    value = 100*c(c(as.matrix(simulationMarketStatic)), c(as.matrix(simulationMarketDynamic)), c(as.matrix(M_t.predict))), 
                                    predictionType = c(rep("StaticSimulation",88), rep("DynamicSimulation", 88), rep("Realization",88)),                         
                                    share = rep(ID, 3))

ggplot(plotBetweenSimulation, aes(x = id, y = value)) + geom_point() + geom_line() + facet_grid( share ~ predictionType) + xlab("time") + ylab("market share in %")
ggplot(plotBetweenSimulation, aes(x = id, y = value, colour=predictionType)) + geom_line() + geom_point() +facet_grid( share~. ) + xlab("time") + ylab("market share in %")

