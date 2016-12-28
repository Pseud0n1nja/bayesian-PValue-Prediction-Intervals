##required packages
#to install, run: 
#install.packages(c("actuar", "LearnBayes"))
require(actuar)
require(LearnBayes)
# desired coverage for Prediction Interval
alpha = 0.95 
#observed upper.tail p-value
p.obt = 0.000216/2
#Sample Size for hypothesis test
N <- 1000

# prior parameters for Sqrt(N).Delta.over.Sigma ~ Normal(m0, s00)
m00 = 0 #mean
v00 <- 1 #variance
#computing interval through helper script
# setwd("~/Documents/summer2016/friendly script/")
source("helperScript.R")
#returns Baysian Discrete Prediction Interval according to OA Vsevolozhskaya, G Ruiz, DV Zaykin
#also returns Frequentist Confidence Interval according to Lazerroni
summaryInfo<-computeInterval(statisticDistribution="student's t",p.obt=p.obt,N=N,
    meanPrior=m00,varPrior=v00)#data frame object
#write summaryInfo to a txt file:
write.table(summaryInfo, sep="        \t",file="summaryInfo.txt",row.names=F)