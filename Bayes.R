install.packages("readxl")
install.packages("invgamma")
install.packages("mvtnorm")
library(mvtnorm)
install.packages("MASS")
library(readxl)
library(reshape2)
library(plyr)
library(dplyr)
library(invgamma)
library(xlsx)
library(MASS)


# Import  data base and select column
sales <- read_excel("C:/Users/Sairam/Downloads/data/sales.xls")
price <- read_excel("C:/Users/Sairam/Downloads/data/price.xls")
displ <- read_excel("C:/Users/Sairam/Downloads/data/displ.xls")
coupon <- read_excel("C:/Users/Sairam/Downloads/data/coupon.xls")

week = 1:nrow(sales)

sales_DF = melt(data = data.frame(week,sales), id.vars = "week", measure.vars = colnames(data.frame(week,sales))[-1],value.name = "sales")
price_DF = melt(data = data.frame(week,price), id.vars = "week", measure.vars = colnames(data.frame(week,price))[-1],value.name = "price")
displ_DF = melt(data = data.frame(week,displ), id.vars = "week", measure.vars = colnames(data.frame(week,displ))[-1],value.name = "displ")
coupon_DF = melt(data = data.frame(week,coupon), id.vars = "week", measure.vars = colnames(data.frame(week,coupon))[-1],value.name = "coupon")


brand_DF = join_all(list(sales_DF,price_DF,displ_DF,coupon_DF), by=c('week','variable'), type='left')
brand_DF = brand_DF[c('variable','week', 'sales',  'price', 'displ', 'coupon')]
brand_DF$sales = log(brand_DF$sales)
brand_DF$price = log(brand_DF$price)

#take only the Brand85
target_DF = brand_DF %>% filter(brand_DF$variable == 'brand85')

# add ones in X
brand_DF_ones = data.frame(replicate(nrow(brand_DF),1))

#Split X and y
X =target_DF[,4:6]
X =cbind(data.frame(replicate(nrow(X),1)), X)
X = data.matrix(X)
y = as.data.frame(target_DF[,3],col.names =c('sales'))
y= data.matrix(y)

draw_MVNormal <- function(y,X, sigmasq, B ,b) {
  
  
  sigma = sigmasq*(solve(t(X)%*%X + solve(B)))
  L <-chol(sigma)
  std_norm_draws <- rnorm(length(b),mean = 0 , sd =1)  
  draw = beta1 +t(L)%*%std_norm_draws
  
  return(draw )
}


#Draw values from inverse gamma
draw_InvGamma2 <- function(y,X,betas ,B, b) {
  dof= (dim(X)[1]  +  dim(X)[2] )
  res= t(y - X%*%betas)%*%(y - X%*%betas)  + t((b-betas))%*%solve(B)%*%(b-betas)
  draw_chi_sq = rchisq(1, df=dof) 
  
  return(as.numeric(res/draw_chi_sq))
}

set.seed(1402)

N = dim(X)[1]
k = dim(X)[2]

mu <- 10
delta <- 10

#Compare with OLS when B is large
lm <-lm(y ~ X)

#Prior
b = c(0,0,0,0)
B = diag(length(b))*1

#likelihood beta OLS
betas = solve(t(X)%*%X) %*% (t(X)%*%y)

#taking one initial value for sigma^2
sigmasq = rinvgamma(n=1, shape=10/2, rate=10/2)
draw_var = sigmasq


beta1 <- solve(t(X)%*%X + solve(B)) %*% (t(X)%*%y + solve(B)%*%b)

draw_beta0 <- numeric()
draw_beta1 <- numeric()
draw_beta2 <- numeric()
draw_beta3 <- numeric()
draw_var_Array <- numeric()

#Gibbs Sampling
#burn in
burnin = 10000
for (i in 1:burnin)
{
  draw_beta = draw_MVNormal(y,X,draw_var,B,b)
  draw_var =  draw_InvGamma2(y,X,draw_beta,B,b)
}

#Running

running = 100000 + burnin
for (i in (burnin+1):running)
{
  
  if( i %% 10000 == 0 ){print(i/100000)}
  draw_beta = draw_MVNormal(y,X,draw_var,B,b)
  draw_beta0<-c(draw_beta0 ,draw_beta[1])
  draw_beta1<-c(draw_beta1 ,draw_beta[2])
  draw_beta2<-c(draw_beta2 ,draw_beta[3])
  draw_beta3<-c(draw_beta3 ,draw_beta[4])
  
  
  draw_var = draw_InvGamma2(y,X,draw_beta,B,b)
  draw_var_Array <- c(draw_var_Array,draw_var)
}


#thin value = 5
range_index = seq.int(from=burnin+1, to=length(draw_beta0), by=5)
print(mean(draw_beta0[range_index]))
print(mean(draw_beta1[range_index]))
print(mean(draw_beta2[range_index]))
print(mean(draw_beta3[range_index]))
print(mean(draw_var_Array[range_index]))

print(sd(draw_beta0[range_index]))
print(sd(draw_beta1[range_index]))
print(sd(draw_beta2[range_index]))
print(sd(draw_beta3[range_index]))
print(sd(draw_var_Array[range_index]))

#posterior odds
sum(draw_beta1>0)/sum(draw_beta1<0)
