# 3 psychological variables
# 4 academic performance variables
# gender variable: 1 stands for female
library(CCA)

psacad = read.csv(file.choose(), header = T) # open PsycAcad.csv
head(psacad)
summary(psacad)
dim(psacad)

X <- psacad[,c(1,2,3)]
Y <- psacad[,4:7]
cor(X)
cor(Y)
cor(X,Y)
pairs(X)
pairs(Y)
pairs(psacad)

cca.psacad = cc(X,Y)
cca.psacad$cor # canonical correlations

cca.psacad$xcoef; cca.psacad$ycoef # canonical coefficients

plot(cca.psacad$scores$xscores[,1],cca.psacad$scores$yscores[,1]) # noticeble correlation
plot(cca.psacad$scores$xscores[,2],cca.psacad$scores$yscores[,2]) # barely noticeble correlation
plot(cca.psacad$scores$xscores[,3],cca.psacad$scores$yscores[,3]) # no visible correlation

n <- length(X[,1]) 
-(n-1/2*(3+4+3))*log((1-cca.psacad$cor[1]^2)*(1-cca.psacad$cor[2]^2)*(1-cca.psacad$cor[3]^2)) # test statistic
qchisq(0.95, 12) # critical value
-(n-1/2*(3+4+3))*log((1-cca.psacad$cor[2]^2)*(1-cca.psacad$cor[3]^2))
qchisq(0.95, 6)
