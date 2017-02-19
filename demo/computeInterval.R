##required packages
#to install, run: 
#install.packages(c("actuar", "LearnBayes"))
require(actuar)
require(LearnBayes)
# desired coverage for Prediction Interval
alpha = 0.95 
#observed upper.tail p-value
p.obt = 0.000108
# Sample Sizes
# assuming two-sample Z or t test
# for one-sample tests, set N <- N1 + N2
N1 <- 300
N2 <- 200
N <- 1 / (1/N1 + 1/N2)
sumN = N1 + N2

# prior parameters for the normal mean
m00 = 0 #mean
v00 <- 0.2 #variance

#computing interval through helper script
# setwd("~/Documents/summer2016/friendly script/")
source("helperScript.R")
#returns Baysian Discrete Prediction Interval according to OA Vsevolozhskaya, G Ruiz, DV Zaykin
#also returns Frequentist Prediction Interval
## statisticDistribution=c("normal","Student")
## p.obt = the observed upper-tail p-value
## N = sample size for test statistic
## meanPrior = the prior mean for the test statistic's mean or non-centrality parameter distribution
## varPrior = the prior variance for the test statistic's mean or non-centrality parameter distribution
summaryInfo<-computeInterval(statisticDistribution="Student", p.obt=p.obt,N=N, sumN=sumN,
    meanPrior=m00,varPrior=v00)#data frame object
#write summaryInfo to a file:
write.csv(summaryInfo,file="summaryInfo.csv",row.names=F)
