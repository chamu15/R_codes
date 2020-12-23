#------------------------------------------------------------------------------
# exercise 1

# I will make it for the first 100 observations. Making it together for 200 
# does not make sense, since there are two distinct groups
bank <- read.csv(file.choose(), header = T)

bank <- bank[1:100,1:6]
head(bank)
summary(bank)
cor(bank)
pairs(bank) 

factanal(bank, 2)
# something that we could expect, the 1st factor is about the size of the note
# the 2nd factor is about the size of the image on the note
#N.B.: You can do it with 3 factors as well, but this is an exact solution without dimension reduction

#------------------------------------------------------------------------------
# exercise 2

crimes <- read.table(file.choose(), header = F) # put header = F

head(crimes)
dim(crimes)
summary(crimes) 
options(digits = 4) # set number of digits = 4, to reduce the size of the matrix, so it becomes readable
cor(crimes)

pairs(crimes) # looks like there are some outliers..

par(mfrow = c(3,3)) # let's check normality, there might be some skewness in first two variables
for(i in 1:9) qqnorm(crimes[,i]) 
for(i in 1:9) hist(crimes[,i], breaks = 20) # variables 1,2,5 are skewed, let's make log transforms

t.crimes <- crimes
t.crimes[,1] <- log(t.crimes[,1])
t.crimes[,2] <- log(t.crimes[,2])
t.crimes[,5] <- log(t.crimes[,5])

par(mfrow=c(3,3)) # check again
for(i in 1:9) qqnorm(t.crimes[,i]) # not perfect but better
for(i in 1:9) hist(t.crimes[,i], breaks = 20) # the same

factanal(t.crimes[,-c(10,11)], 1) # we can clearly see that one factor is not enough (test p-value) 
# and V1 has very high \psi variance.

factanal(t.crimes[,-c(10,11)], 2) # difficult to interprete and two factors is not enough (test p-value)
# looks like V1 stands out, and can be possibly removed as irrelevant

factanal(t.crimes[,-c(1,10,11)], 4) #1st factor is related to violent crimes
# 2nd factor robbery-type non-violent crimes
# 3rd factor is dominated by autotheft and robbery appears everywhere, but a bit heavier in 3rd
# 4th is the population. 

# let's separate population and analyze only the crimes
factanal(t.crimes[,-c(1,2,10,11)], 3)
# same remarks as above, and the loadings are also close.

#------------------------------------------------------------------------------
# exercise 3
library(CCA)
cars <- read.csv(file.choose(), header = T)

head(cars)
summary(cars) # the first column is not a variable!
dim(cars)
n <- 23
  
row.names(cars)= cars[,1] # the first column become rownames attribute
cars = cars[,-1] # delete the first column

cov(cars)
pairs(cars) # value and price are strongly correlated
# other variables are also correlated with price and value
# let's see what CCA tells us

X <- cars[,3:4] # take price and value
Y <- cars[,-c(3:4)] # extract everything except columns 3 and 4, which are price and value

cor(X)
cor(Y)
cor(X,Y) # number of canonical correlations is 2
cca.cars <- cc(X,Y)
cca.cars$cor # very strong relation between groups 0.98

cca.cars$xcoef; cca.cars$ycoef
# Price is positively related to Economy and negatively to everything else
# Value is negatively related to economy and positively to everything else

# plot needs both commands
plot(cca.cars$scores$xscores[,1],cca.cars$scores$yscores[,1], type="n")
text(x = cca.cars$scores$xscores[,1], y = cca.cars$scores$yscores[,1],
     labels = row.names(cars), cex=.75) 
#we can see the group of cheap cars on one side and expensive cars on another
# \eta can be interpreted as value index and \phi can be interpreted as quality index
# note that eigenvectors can be multiplied by -1  for the graph to make sense (if necessary)

plot(cca.cars$scores$xscores[,2],cca.cars$scores$yscores[,2], type="n")
text(x = cca.cars$scores$xscores[,2], y = cca.cars$scores$yscores[,2],
     labels = row.names(cars), cex=.75) 
# \phi_2 consists of mainly Economy, Service and Easy (handling)

-(n-1/2*(2+6+3))*log((1-cca.cars$cor[1]^2)*(1-cca.cars$cor[2]^2)) # test statistic
qchisq(0.95, 12) # critical value
-(n-1/2*(2+6+3))*log((1-cca.cars$cor[2]^2))
qchisq(0.95, 5)
# reject twice

#------------------------------------------------------------------------------
# exercise 4
bank <- read.csv(file.choose(), header = T)

bank <- bank[101:200,1:6]
head(bank)
summary(bank)
cov(bank)
pairs(bank) 
n <- 100

X <- bank[,1:3]
Y <- bank[,4:6]
cor(X)
cor(Y)
cor(X,Y) 

cca.bank <- cc(X,Y)
cca.bank$cor # the relationship is quite weak
cca.bank$xcoef; cca.bank$ycoef

plot(cca.bank$scores$xscores[,1],cca.bank$scores$yscores[,1]) # actually looks almost like pure noise
plot(cca.bank$scores$xscores[,2],cca.bank$scores$yscores[,2]) # two noticeble groups, but hardly interpretable

-(n-1/2*(3+3+3))*log((1-cca.bank$cor[1]^2)*(1-cca.bank$cor[2]^2)*(1-cca.bank$cor[3]^2)) # test statistic
qchisq(0.95, 9) # significant
