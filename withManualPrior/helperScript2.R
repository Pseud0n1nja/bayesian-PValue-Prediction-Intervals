require(actuar)
require(LearnBayes)
####################
computeInterval<-function(statisticDistribution,p.obt,priorBinMeans,priorBinWeights,df1=NULL,df2=NULL){
  if(statisticDistribution=="Normal"){
    z.obt=qnorm(p=p.obt,lower.tail=F)
    return(computeInterval_normDistr(z.obt=z.obt,fud=priorBinMeans,px=priorBinWeights))
  }else if(statisticDistribution=="Student's t"){
    #observed t-statistic; note: must be upper tail t-statistic
    t.obt=qt(p=p.obt,df=df1,lower.tail=F)
    return(computeInterval_tDistr(t.obt=t.obt,fud=priorBinMeans,px=priorBinWeights,df=df1))
  }else if(statisticDistribution=="F"){
    if(min(priorBinMeans)<0){
      print("ERROR: bin means must be >= 0. Support of any F-distribution is 0 to Infty.")
      return(NULL)
    }
    #observed F-statistic; note: must be upper tail F-statistic
    F.obt=qf(p=p.obt,df1=df1,df2=df2,lower.tail=F)
    return(computeInterval_FDistr(F.obt=F.obt,fud=priorBinMeans,px=priorBinWeights,df1=df1,df2=df2))
  }else if(statisticDistribution=="Chi-Squared"){
    if(min(priorBinMeans)<0){
      print("ERROR: bin means must be >= 0. Support of any Chisq-distribution is 0 to Infty.")
      return(NULL)
    }
    #observed Chisq-statistic; note: must be upper tail Chisq-statistic
    Chisq.obt=qchisq(p=p.obt,df=df1,lower.tail=F)
    return(computeInterval_ChisqDistr(Chisq.obt=Chisq.obt,fud=priorBinMeans,px=priorBinWeights,df=df1))
  }else{
    print("Error: Unkown test statstic distribution speciefied.")
    return(NULL)
  }
}

##################################
###Stuff Related to t-Statistic###
##################################
##
F_tDistr = function(x,w,df,ncp){ 
  prob<- Inf;
  if(x>=sum(w*ncp)){
    prob<-sum( w*(1-pt(q=x,df=df,ncp=ncp,lower.tail=F)) ) 
  }else{#split into cases to avoid precision warnings
    prob<-sum( w*pt(q=x,df=df,ncp) ) 
  }
  return(prob)
}
# provide an initial bracket for the quantile. default is c(-10000,10000). 
F_inv_tDistr = function(p,w,df,ncp,br=c(-100000,100000))
{
  G = function(x) {F_tDistr(x,w,df,ncp) - p}
  return( uniroot(G,br)$root ) 
}
##
#returns Baysian Discrete Prediction Interval according to OA Vsevolozhskaya, G Ruiz, DV Zaykin
#also returns Frequentest Confidence Interval according to Lazerroni
computeInterval_tDistr<-function(t.obt,fud,px,df){
  ##posterior weights and mean estimate given observed p-value/t-statistic
  wts_post<-px*dt(x=t.obt,ncp=fud,df=df)/sum(px*dt(x=t.obt,ncp=fud,df=df))
  #Note: posterior mean is a weighted average using updated weights
  th.discr<-sum(fud*wts_post)
  ##
  #For Discrete Bayesian PI##
  Lower.Discr<-F_inv_tDistr(p=(1-alpha)/2,w=wts_post,df=df,ncp=fud)
  Upper.Discr<-F_inv_tDistr(p=1-(1-alpha)/2,w=wts_post,df=df,ncp=fud)
  ##################
  #For Lazerroni PI#
  Lower.L<-NA#qt((1-alpha)/2, ncp=t.obt, df=df)
  Upper.L<-NA#qt(1-(1-alpha)/2, ncp=t.obt, df=df)
  ##################
  
  ###################
  ###Summary Table###
  ###################
  
  
  summary_table<-data.frame(
    pInterval.type=c("Discrete Bayesian","Lazerroni"),
    lower.tStat.bound=c(Lower.Discr,Lower.L),
    upper.tStat.bound=c(Upper.Discr,Upper.L),
    lower.pInterval.bound=
      pt(c(Upper.Discr,Upper.L),df=df,lower.tail=F),
    upper.pInterval.bound=
      pt( c(Lower.Discr,Lower.L),df=df,lower.tail=F)
  )
  return(summary_table)
}

###################################
###Stuff Related to Z-Statistics###
###################################
##
# taken from: http://stats.stackexchange.com/questions/14481/quantiles-from-the-combination-of-normal-distributions
# evaluate the function at the point x, where the components 
# of the mixture have weights w, means stored in u, and std deviations
# stored in s - all must have the same length.
F_normDistr = function(x,w,u,s){ 
  prob<- Inf;
  if(x>=sum(w*u)){
    prob<-sum( w*(1-pnorm(q=x,mean=u,sd=s,lower.tail=F)) ) 
  }else{#split into cases to avoid precision warnings
    prob<-sum( w*pnorm(q=x,mean=u,sd=s) ) 
  }
  return(prob)
}
# provide an initial bracket for the quantile. default is c(-10000,10000). 
F_inv_normDistr = function(p,w,u,s,br=c(-100000,100000))
{
  G = function(x) {F_normDistr(x,w,u,s) - p}
  return( uniroot(G,br)$root ) 
}
##

#returns Baysian Discrete Prediction Interval for Z-statistics according to OA Vsevolozhskaya, G Ruiz, DV Zaykin
#also returns Frequentest Confidence Interval according to Lazerroni
computeInterval_normDistr<-function(z.obt,fud,px){
  ##posterior weights and mean estimate given observed p-value/t-statistic
  wts_post<-px*dnorm(x=z.obt,mean=fud,sd=1)/sum(px*dnorm(x=z.obt,mean=fud,sd=1))
  #Note: posterior mean is a weighted average using updated weights
  th.discr<-sum(fud*wts_post)
  ##
  #For Discrete Bayesian PI##
  Lower.Discr<-F_inv_normDistr(p=(1-alpha)/2,w=wts_post,u=fud,s=1)
  Upper.Discr<-F_inv_normDistr(p=1-(1-alpha)/2,w=wts_post,u=fud,s=1)
  ##################
  #For Lazerroni PI#
  Lower.L<-qnorm((1-alpha)/2, mean=z.obt, sd=sqrt(2))
  Upper.L<-qnorm(1-(1-alpha)/2, mean=z.obt, sd=sqrt(2))
  ##################
  
  ###################
  ###Summary Table###
  ###################
  
  
  summary_table<-data.frame(
    pInterval.type=c("Discrete Bayesian","Lazerroni"),
    lower.zStat.bound=c(Lower.Discr,Lower.L),
    upper.zStat.bound=c(Upper.Discr,Upper.L),
    lower.pInterval.bound=
      pnorm(c(Upper.Discr,Upper.L),lower.tail=F),
    upper.pInterval.bound=
      pnorm( c(Lower.Discr,Lower.L),lower.tail=F)
    
  )
  return(summary_table)
}

#####################################
###Funcions related to F-statistic###
#####################################
F_FDistr = function(x,w,df1,df2,ncp){ 
  prob<- Inf;
  if(x>=sum(w*ncp)){
    prob<-sum( w*(1-pf(q=x,df1=df1,df2=df2,ncp=ncp,lower.tail=F)) ) 
  }else{#split into cases to avoid precision warnings
    prob<-sum( w*pf(q=x,df1=df1,df2=df2,ncp=ncp) )
  }
  return(prob)
}
# provide an initial bracket for the quantile. default is c(-10000,10000). 
F_inv_FDistr = function(p,w,df1,df2,ncp,br=c(-100000,100000))
{
  G = function(x){F_FDistr(x,w,df1,df2,ncp) - p}
  return( uniroot(G,br)$root ) 
}
##
#returns Baysian Discrete Prediction Interval according to OA Vsevolozhskaya, G Ruiz, DV Zaykin
#also returns Frequentest Confidence Interval according to Lazerroni
computeInterval_FDistr<-function(F.obt,fud,px,df1,df2){
  ##posterior weights and mean estimate given observed p-value/t-statistic
  wts_post<-px*df(x=F.obt,ncp=fud,df1=df1,df2=df2)/sum(px*df(x=F.obt,ncp=fud,df1=df1,df2=df2))
  #Note: posterior mean is a weighted average using updated weights
  th.discr<-sum(fud*wts_post)
  ##
  #For Discrete Bayesian PI##
  Lower.Discr<-F_inv_FDistr(p=(1-alpha)/2,w=wts_post,df1=df1,df2=df2,ncp=fud)
  Upper.Discr<-F_inv_FDistr(p=1-(1-alpha)/2,w=wts_post,df1=df1,df2=df2,ncp=fud)
  ##################
  #For Lazerroni PI#
  Lower.L<-NA#qt((1-alpha)/2, ncp=t.obt, df=df)
  Upper.L<-NA#qt(1-(1-alpha)/2, ncp=t.obt, df=df)
  ##################
  
  ###################
  ###Summary Table###
  ###################
  summary_table<-data.frame(
    pInterval.type=c("Discrete Bayesian","Lazerroni"),
    lower.FStat.bound=c(Lower.Discr,Lower.L),
    upper.FStat.bound=c(Upper.Discr,Upper.L),
    lower.pInterval.bound=
      pf(c(Upper.Discr,Upper.L),df1=df1,df2=df2,lower.tail=F),
    upper.pInterval.bound=
      pf( c(Lower.Discr,Lower.L),df1=df1,df2=df2,lower.tail=F)
  )
  return(summary_table)
}
##

#####################################
###Funcions related to F-statistic###
#####################################
F_ChisqDistr = function(x,w,df,ncp){ 
  prob<- Inf;
  if(x>=sum(w*ncp)){
    prob<-sum( w*(1-pchisq(q=x,df=df,ncp=ncp,lower.tail=F)) ) 
  }else{#split into cases to avoid precision warnings
    prob<-sum( w*pchisq(q=x,df=df,ncp=ncp) )
  }
  return(prob)
}
# provide an initial bracket for the quantile. default is c(-10000,10000). 
F_inv_ChisqDistr=function(p,w,df,ncp,br=c(-100000,100000))
{
  G = function(x){F_ChisqDistr(x,w,df,ncp) - p}
  return( uniroot(G,br)$root ) 
}
##
#returns Baysian Discrete Prediction Interval according to OA Vsevolozhskaya, G Ruiz, DV Zaykin
#also returns Frequentest Confidence Interval according to Lazerroni
computeInterval_ChisqDistr<-function(Chisq.obt,fud,px,df){
  ##posterior weights and mean estimate given observed p-value/t-statistic
  wts_post<-px*dchisq(x=Chisq.obt,ncp=fud,df=df)/sum(px*dchisq(x=Chisq.obt,ncp=fud,df=df))
  #Note: posterior mean is a weighted average using updated weights
  th.discr<-sum(fud*wts_post)
  ##
  #For Discrete Bayesian PI##
  Lower.Discr<-F_inv_ChisqDistr(p=(1-alpha)/2,w=wts_post,df=df,ncp=fud)
  Upper.Discr<-F_inv_ChisqDistr(p=1-(1-alpha)/2,w=wts_post,df=df,ncp=fud)
  ##################
  #For Lazerroni PI#
  Lower.L<-NA#qt((1-alpha)/2, ncp=t.obt, df=df)
  Upper.L<-NA#qt(1-(1-alpha)/2, ncp=t.obt, df=df)
  ##################
  
  ###################
  ###Summary Table###
  ###################
  summary_table<-data.frame(
    pInterval.type=c("Discrete Bayesian","Lazerroni"),
    lower.ChisqStat.bound=c(Lower.Discr,Lower.L),
    upper.FStat.bound=c(Upper.Discr,Upper.L),
    lower.pInterval.bound=
      pchisq(c(Upper.Discr,Upper.L),df=df,lower.tail=F),
    upper.pInterval.bound=
      pchisq( c(Lower.Discr,Lower.L),df=df,lower.tail=F)
  )
  return(summary_table)
}
####################