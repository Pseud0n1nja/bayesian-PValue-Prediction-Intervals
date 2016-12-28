##########
require(actuar)
require(LearnBayes)
#

####################
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
####################
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
####################
computeInterval<-function(statisticDistribution="Normal",p.obt,N,meanPrior,varPrior){
  if(statisticDistribution=="normal"){
    z.obt=qnorm(p=p.obt,lower.tail=F)
    return(computeInterval_normDistr(z.obt=z.obt,N=N,meanPrior=meanPrior,varPrior=varPrior))
  }else if(statisticDistribution=="student's t"){
    #observed t-statistics; note: must be upper tail t-statistic
    t.obt=qt(p=p.obt,df=N-1,lower.tail=F)
    return(computeInterval_tDistr(t.obt=t.obt,N=N,meanPrior=meanPrior,varPrior=varPrior))
  }else{
    print("Error: Unkown test statstic distribution speciefied.")
  }
}


#returns Baysian Discrete Prediction Interval according to OA Vsevolozhskaya, G Ruiz, DV Zaykin
#also returns Frequentest Confidence Interval according to Lazerroni
computeInterval_tDistr<-function(t.obt,N,meanPrior,varPrior){
  #How wide every bin should be, except leftmost and rightmost bins
  #the narrower, the more precision; 
  #However, negligible precision increase for smaller bin widths than below
  stp = sqrt(varPrior)/1000 
  #because we cannot get quantiles (the bin means) for the entire 
  #domain of the standard normal (the real line), we settle with bin means 
  #from the (1-trunc)th quantile to the (trunc)th quantile
  #the sum of our bin weights "px" will equal approximately "trunc"
  trunc <- 0.99999999999999
  #the probability associated with every bin
  px <- discretize(pnorm(x, mean=meanPrior, sd=sqrt(varPrior)),
                   method = "rounding",
                   from = qnorm(1-trunc, mean=meanPrior, sd=sqrt(varPrior)),
                   to = qnorm(trunc, mean=meanPrior, sd=sqrt(varPrior)), step = stp)
  #the mean associated with every bin
  fud <- seq(from=qnorm(1-trunc,mean=meanPrior, sd=sqrt(varPrior)), 
             to=qnorm(trunc, mean=meanPrior, sd=sqrt(varPrior))-stp, by=stp)
  ###
  ###########
  ###
  ##posterior weights and mean estimate given observed p-value/t-statistic
  wts_post<-px*dt(x=t.obt,ncp=fud,df=N-1)/sum(px*dt(x=t.obt,ncp=fud,df=N-1))
  #Note: posterior mean is a weighted average using updated weights
  th.discr<-sum(fud*wts_post)
  ##
  #For Discrete Bayesian PI##
  Lower.Discr<-F_inv_tDistr(p=(1-alpha)/2,w=wts_post,df=N-1,ncp=fud)
  Upper.Discr<-F_inv_tDistr(p=1-(1-alpha)/2,w=wts_post,df=N-1,ncp=fud)
  ##################
  #For Lazerroni PI#
  Lower.L<-NA#qt((1-alpha)/2, ncp=t.obt, df=N-1)
  Upper.L<-NA#qt(1-(1-alpha)/2, ncp=t.obt, df=N-1)
  ##################
  
  ###################
  ###Summary Table###
  ###################
  
  
  summary_table<-data.frame(
    pInterval.type=c("Discrete Bayesian","Lazerroni"),
    lower.tStat.bound=c(Lower.Discr,Lower.L),
    upper.tStat.bound=c(Upper.Discr,Upper.L),
    lower.pInterval.bound=
      pt(c(Upper.Discr,Upper.L),df=N-1,lower.tail=F),
    upper.pInterval.bound=
      pt( c(Lower.Discr,Lower.L),df=N-1,lower.tail=F)
  )
  return(summary_table)
}

#returns Baysian Discrete Prediction Interval for Z-statistics according to OA Vsevolozhskaya, G Ruiz, DV Zaykin
#also returns Frequentest Confidence Interval according to Lazerroni
computeInterval_normDistr<-function(z.obt,N,meanPrior,varPrior){
  #How wide every bin should be, except leftmost and rightmost bins
  #the narrower, the more precision; 
  #However, negligible precision increase for smaller bin widths than below
  stp = sqrt(varPrior)/1000 
  #because we cannot get quantiles (the bin means) for the entire 
  #domain of the standard normal (the real line), we settle with bin means 
  #from the (1-trunc)th quantile to the (trunc)th quantile
  #the sum of our bin weights "px" will equal approximately "trunc"
  trunc <- 0.99999999999999
  #the probability associated with every bin
  px <- discretize(pnorm(x, mean=meanPrior, sd=sqrt(varPrior)),
                   method = "rounding",
                   from = qnorm(1-trunc, mean=meanPrior, sd=sqrt(varPrior)),
                   to = qnorm(trunc, mean=meanPrior, sd=sqrt(varPrior)), step=stp)
  #the mean associated with every bin
  fud <- seq(from=qnorm(1-trunc,mean=meanPrior, sd=sqrt(varPrior)), 
             to=qnorm(trunc, mean=meanPrior, sd=sqrt(varPrior))-stp, by=stp)
  ###
  ###########
  ###
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