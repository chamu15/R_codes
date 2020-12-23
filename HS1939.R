# data file is "HolzingerSwineford.csv"
library("psych")
library("robustbase")
HS = read.csv(file.choose(), header=T)
# Variables:
# Paragraph - Paragraph comprehension
# Sentence - Sentence completion
# Word - Word meaning
# Add - Speeded addition
# Dots - Speeded counting of dots
#
# the first 156 students are from Pasteur school
# the other 145 are from Grant-White school
dim(HS)
head(HS)
summary(HS)
var(HS)
cor(HS)

GW <- HS[157:301,]
fa.GW <- factanal(GW, 2, rotation="none") # this function is from basic stats package
fa.GW
load.GW <- fa.GW$loadings[,1:2]
plot(load.GW,type="n", xlim = c(-1,1), ylim = c(-1,1), main = "Factor plot without rotation"); abline(h=0); abline(v=0)
text(load.GW,labels=names(HS),cex=.7) 

fa.GWrd = factanal(GW, 2) # varimax is the default rotation
fa.GWrd
load.GWrd <- fa.GWrd$loadings[,1:2]
plot(load.GWrd,type="n", xlim = c(-1,1), ylim = c(-1,1), main = "Factor plot with rotation") ; abline(h=0); abline(v=0) 
text(load.GWrd,labels=names(HS),cex=.7) 

# there is another implementation in psych package with a lot of advanced options
# check yourself
cov.GW = cov(GW)
fa.GW = fa(cov.GW, nfactors = 2, rotate = "varimax", covar = T)
fa.GW

cor.GW = cor(GW)
fa.GW = fa(cor.GW, nfactors = 2, rotate = "varimax")
fa.GW

mcd.GW <- covMcd(GW)$cov
fa.GWmcd <- fa(mcd.GW, nfactors = 2, rotate = "varimax", covar = T)
fa.GWmcd

mcd.GW <- covMcd(GW, cor=T)$cor
fa.GWmcd <- fa(mcd.GW, nfactors = 2, rotate = "varimax")
fa.GWmcd

