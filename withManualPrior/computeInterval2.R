##importing weights and bin means
#setwd("~/Documents/summer2016/friendly script/")
# install.packages("gdata")
require(gdata)
##example of discretized prior for z-,t-statistics. The Support is the real line.
# df = read.csv(file = "manualPrior.csv", header = TRUE)
##example of discretized prior for F-,Chisq-statistics. The Support is 0 to Inf.
# df = read.csv(file = "manualPrior2.csv", header = TRUE)

df = read.csv(file = "manualPrior2.csv", header = TRUE)
df$bin.mean = (df$lower.bin.bound + df$upper.bin.bound) / 2
df$bin.prob <- df$bin.prob/sum(df$bin.prob)

# desired coverage for Prediction Interval
alpha = 0.95 
# observed upper.tail p-value
p.obt = 0.0000108
# Sample sizes
# assuming two-sample Z or t test
# for one-sample tests, set N <- N1 + N2
N1 <- 300
N2 <- 200
N <- 1 / (1/N1 + 1/N2)
# Set degrees of freedom; geared toward F and t. May need to set manually
sumN = N1 + N2
if(N == sumN) {
    Df <- N-1
} else {
    Df <- sumN-2
}
source("helperScript2.R")
#returns Baysian Discrete Prediction Interval according to OA Vsevolozhskaya, G Ruiz, DV Zaykin
#also returns Frequentist Interval as a data frame object
## statisticDistribution=c("Normal","Student", "F","Chi-Squared")
## p.obt = the observed upper-tail p-value
## priorBinMeans = the means or non-centrality values for the prior distribution
## priorBinWeights = the weights/probability that the mean/non-centrality parameter falls in the bin
## df1 = degrees of freedom associated with test if applicable; argument may be omitted if 
#### not applicable. If F-test, then df1=numerator degrees of freedom
## df2 = denominator degrees of freedom for F-test. If not applicable, argument may be omitted.
(summaryInfo<-computeInterval(statisticDistribution="Normal",p.obt=p.obt, N=N,
      priorBinMeans=df$bin.mean,priorBinWeights=df$bin.prob,df1=1,df2=Df))
##
#write summaryInfo to a csv file:
write.csv(summaryInfo, file="summaryInfo.csv",row.names=F)

