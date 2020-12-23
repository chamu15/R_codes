# Course: Advanced Marketing Models 2014-2015
# Code for programming assignment 1

# -----------------------------------------------------------------------  
# -- 0 Clear workspace
# -----------------------------------------------------------------------

rm(list = ls())

# -----------------------------------------------------------------------  
# -- 1 Set user-defined functions
# -----------------------------------------------------------------------

# -- Function to obtain OLS estimates
estimateOLS <- function(y, x) {
  OLSestimates = solve(t(x)%*%x)%*%t(x)%*%y
  return(OLSestimates)
}

# Initialize key variables
numBrands = 4
depVar = "VOL"
explanatoryVarNames = c("PRICE","PROMO1","PROMO2")

# Optional path for the data files
path <- ""

categoryNames <- c("beer", "bjc", "cer", "che", "coo","cso","did","frd","frj","lnd", "ptw", "rfj","sdr", "sna","tpa")

# -----------------------------------------------------------------------  
# -- 2 Load data
# -----------------------------------------------------------------------
numCats = length(categoryNames)

# Note that we first converted the .xls files to .txt files so that the files can be easily (and quickly) loaded into R
# Put all dataframes into one list
dataSets <- list()
for (cat in 1:numCats) {
  tmp = read.delim(paste0(path,categoryNames[[cat]],"data1.txt"))
  # add the data set of a category with the category name attached to it
  dataSets = c(dataSets, setNames(list(tmp),categoryNames[cat]))
}


# -----------------------------------------------------------------------  
# -- 3 Estimate basic sales models
# -----------------------------------------------------------------------

# Initialize matrix in which the OLS estimates for each brand in each product category are stored together with explanatory variables at the brand level
OLSest_results = matrix(0, nrow = numCats*numBrands, ncol = 2+length(explanatoryVarNames))
colnames(OLSest_results) = c("PriceElas", "avg_share", paste0("avg_", explanatoryVarNames))

# Get the OLS estimates for each brand in each product category
index = 1;
for (cat in 1:numCats) {
  dataCat <- dataSets[[cat]] # get dataset of specific category
  for (brand in 1:numBrands) {
      
    # identify all "normal" observations
    usefulObs = as.logical((dataCat[paste0("VOL",brand)]>0)*(dataCat[paste0("PRICE",brand)]>0))
    usefulObs[!is.finite(usefulObs)] = FALSE
    
    # get the right matrices and perform OLS estimation
    y_brand = as.matrix(dataCat[paste0(depVar,brand)]);
    y_brand = log(y_brand[usefulObs])
    
    x_brand = as.matrix(dataCat[paste0(explanatoryVarNames,brand)])
    x_brand = x_brand[usefulObs,]
    x_brand[,1] = log(x_brand[,1])
    x_brand = cbind(1,x_brand)
  
    olsRes = estimateOLS(y=y_brand, x=x_brand) 
    OLSest_results[index][1] = olsRes[2]

    # add explanatory variables for this brand  
    share = as.matrix(dataCat[paste0("VOL",brand)]/dataCat["CSALES"])
    OLSest_results[index,2]   = mean(share[usefulObs])
    OLSest_results[index,3:5] = t(colMeans(dataCat[paste0(explanatoryVarNames,brand)][usefulObs,]))
    
    index = index + 1;
  }    
}

# -----------------------------------------------------------------------  
# -- 4 A model for the price elasticities of brands
# -----------------------------------------------------------------------
OLSest_priceElas = estimateOLS(y=OLSest_results[,1], x=cbind(1, OLSest_results[,2:5])) 
print(OLSest_priceElas)

