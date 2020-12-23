# Install required packages
library(dplyr) # For data manipulation
library(numDeriv) # For computing Hessian

############################################################################
###### Question 2 ##########################################################
############################################################################


# Read in data
data <- read.csv('C:/Users/Sairam/Downloads/539385.csv')

# Create the X variable
data_week <- distinct(data, Week, .keep_all = TRUE)
X <- data_week$price
X <- cbind(1, X)

# Create the y variable
y <- matrix(data$y, nrow=250, ncol=25)

# Function to calculate the log likelihood
LogL <- function(theta, pi_k, y, X){
  N <- nrow(y) # Define nr of individuals
  K <- nrow(theta) # Define number of clusters
  loglik_tot <- 0 # Initialize loglikelihood value
  for (i in 1:N){
    lik_fi <- 0
    for (k in 1:K){
      # Define logit function
      logit <- function(theta, X){
        exp(theta%*%t(X))/(1+exp(theta%*%t(X)))
      }
      p <- logit(theta[k,],X)
      # Calculate likelihood contribution for individual i for segment k
      # Value of product of all densities f(y_{it}) (t=1,...,25)
      lik_fik <- pi_k[k] * prod(p^y[i,] * (1-p)^(1-y[i,]))
      
      # Calculate likelihood contribution for individual i
      lik_fi <- lik_fi + lik_fik
    }
    loglik_tot <- loglik_tot + log(lik_fi)
  }
  return(loglik_tot)
}

############################################################################
###### Question 3 ##########################################################
############################################################################

EStep <- function(theta, pi_k, y, X){
  N <- nrow(y) # Define nr of individuals
  K <- nrow(theta) # Define number of clusters
  W <- matrix(0, nrow = N, ncol = K) # Initialize the matrix with conditional cluster probabilities
  # We loop over each individual
  for (i in 1:N){
    # First we want to compute the f_yik value for each segment
    f_yik <- c() # Initialize empty array, which can be filled
    for (k in 1:K){
      # Define logit function
      logit <- function(theta, X){
        exp(theta%*%t(X))/(1+exp(theta%*%t(X)))
      }
      p <- logit(theta[k,],X)
      
      f_yik[k] <- prod(p^y[i,] * (1-p)^(1-y[i,]))
    }
    # Now that we have the values of f_yik for each segment we can compute the conditional cluster probabilities
    for (k in 1:K){
      W[i,k] <- pi_k[k]*f_yik[k]/(sum(pi_k*f_yik))
    }
  }
  return(W)
}

############################################################################
###### Question 5 ##########################################################
############################################################################

MStep <- function(W, y, X, theta_start){
  # Update the new estimate of pi
  N <- nrow(W)
  K <- ncol(W)
  pi_k <- c()
  for (k in 1:K){
    pi_k[k] <- 1/N * sum(W[,k])
  }
  # Update the new estimate of theta
  N <- nrow(W)
  K <- ncol(W)
  theta_new <- matrix(0, nrow = K, ncol = 2)
  for (k in 1:K){
    maxim <- function(theta){
      sum_tot <- 0
      for (i in 1:N){
        # Define logit function
        logit <- function(theta, X){
          exp(theta%*%t(X))/(1+exp(theta%*%t(X)))
        }
        p <- logit(theta,X)
        
        # Compute the value within the sum over i for each i
        sum_i <- W[i,k] * sum(y[i,]*log(p) + (1-y[i,])*log(1-p))

        sum_tot <- sum_tot + sum_i
      }
      return(-sum_tot)
    }
    theta_new[k,] <- optim(theta_start[k,], maxim)$par
  }
  params <- list("pi"=pi_k, "theta"=theta_new)
  return(params)
}

############################################################################
###### Question 6 ##########################################################
############################################################################

EM <- function(K, y, X){
  iter <- 1 # Initialize counter
  theta <- matrix(rnorm(2*K, 0, 1), nrow = K, ncol=2) # initialize theta with values random generated from a standard normal
  pi_k <- matrix(1/K, nrow = 1, ncol = K) # Initialize segment probabilities all equal
  # Initilize W matrices, such that the while loop is able to start
  W_0 <- matrix(0, nrow=nrow(y), ncol=K)
  W <- matrix(1, nrow = nrow(y), ncol=K)
  while(max(abs(W-W_0)) > 0.0001){
    # First do the E-step
    W_0 <- W
    W <- EStep(theta, pi_k, y, X)
    # Then the M-step
    params <- MStep(W,y,X,theta)
    theta <- params$theta
    pi_k <- params$pi
    print(iter)
    iter <- iter + 1
    #print(theta)
    #print(pi_k)
  }
  params <- list("pi"=pi_k, "theta"=theta)
  return(params)
}

############################################################################
###### Question 7 ##########################################################
############################################################################

Estimate <- function(K, y, X){
  loglik_best <- -10000 # Set a really low initial best loglik value
  for (i in 1:10){
    print(i)
    EM_algo <- EM(K, y, X)
    pi <- EM_algo$pi
    theta <- EM_algo$theta
    
    loglik_value <- LogL(theta, pi, y, X)
    if (loglik_value > loglik_best){
      pi_best <- pi
      theta_best <- theta
      loglik_best <- loglik_value
    }
    print(pi_best)
    print(theta_best)
    print(loglik_best)
  }
  best <-list("pi"=pi_best, "theta"=theta_best, "loglikelihood"=loglik_best)
  return(best)
}

############################################################################
###### Question 8 ##########################################################
############################################################################

Estimate <- function(K, y, X){
  loglik_best <- -10000 # Set a really low initial best loglik value
  for (i in 1:10){
    print(i)
    EM_algo <- EM(K, y, X)
    pi <- EM_algo$pi
    theta <- EM_algo$theta
    
    loglik_value <- LogL(theta, pi, y, X)
    if (loglik_value > loglik_best){
      pi_best <- pi
      theta_best <- theta
      loglik_best <- loglik_value
    }
    print(pi_best)
    print(theta_best)
    print(loglik_best)
  }
  
  ### Extend function calculate the variance matrix of estimated parameters
  
  # Compute the gamma vector
  gamma_best <- c()
  for (k in 1:K-1){
    gamma_best[k] <- log(pi_best[k]) - log(pi_best[K])
  }
  
  # Create new log likelihood function based on gamma input
  LogL2 <- function(theta){
    K <- length(gamma_best) + 1
    for (k in 1:K-1){
      pi[k] <- exp(gamma_best[k])/(1+sum(exp(gamma_best)))
    }
    pi[K] <- 1/(1+sum(exp(gamma_best)))
    LogL(theta, pi, y, X)
  }
  
  # Calculate the hessian
  hess <- hessian(LogL2, x=theta_best)
  
  # Estimate the covariance matrix
  cov_mat <- solve(-hess)
  
  # Obtain standard errors for parameters
  params_se <- sqrt(diag(cov_mat))
  
  best <-list("pi"=pi_best, "theta"=theta_best, "loglikelihood"=loglik_best, "standard erros"=params_se)
  return(best)
}

############################################################################
###### Question 9 ##########################################################
############################################################################

Estimate(2, y, X)
Estimate(3, y, X)
Estimate(4, y, X)
