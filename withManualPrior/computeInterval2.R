##importing weights and bin means
setwd("~/Documents/summer2016/friendly script/")
# install.packages("gdata")
require(gdata)
##example of discretized prior for z-,t-statistics. The Support is the real line.
# df = read.xls(xls = "manualPrior.xlsx", sheet = 1, header = TRUE)
##example of discretized prior for F-,Chisq-statistics. The Support is 0 to Inf.
df = read.xls(xls = "manualPrior2.xlsx", sheet = 1, header = TRUE)

# desired coverage for Prediction Interval
alpha = 0.95 
#observed upper.tail p-value
p.obt = 0.000216/2
#Sample Size for hypothesis test
N <- 1000

source("helperScript2.R")
#returns Baysian Discrete Prediction Interval according to OA Vsevolozhskaya, G Ruiz, DV Zaykin
#also returns Frequentist Confidence Interval according to Lazerroni as a data frame object
## statisticDistribution=c("Normal","Student's t", "F","Chi-Squared")
## p.obt = the observed upper-tail p-value
## priorBinMeans = the means or non-centrality values for the prior distribution
## priorBinWeights = the weights/probability that the mean/non-centrality parameter falls in the bin
## df1 = degrees of freedom associated with test if applicable; argument may be omitted if 
#### not applicable. If F-test, then df1=numerator degrees of freedom
## df2 = denominator degrees of freedom for F-test. If not applicable, argument may be omitted.
(summaryInfo<-computeInterval(statisticDistribution="F",p.obt=p.obt,
      priorBinMeans=df$bin.mean,priorBinWeights=df$bin.prob,df1=N-1,df2=1))
##
#write summaryInfo to a txt file:
write.table(summaryInfo, sep="        \t",file="summaryInfo.txt",row.names=F)