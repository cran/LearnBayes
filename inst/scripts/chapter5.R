# chapter5.R

library(LearnBayes)

#####################################################
# Section 5.5 Approximations Based on Posterior Modes
#####################################################

 data(cancermortality)
 mycontour(betabinexch0,c(.0001,.003,1,20000),cancermortality)

######################################################
# Section 5.6  The Example
######################################################

 fit=laplace(betabinexch,c(-7,6),cancermortality)
 fit

 npar=list(m=fit$mode,v=fit$var)
 mycontour(lbinorm,c(-8,-4.5,3,16.5),npar)
 title(xlab="logit p", ylab="log K")

 se=sqrt(diag(fit$var))
 fit$mode-1.645*se
 fit$mode+1.645*se


#########################################################
# Section 5.7 Monte Carlo Method for Computing Integrals
#########################################################

 p=rbeta(1000, 14.4, 23.4)
 est=mean(p^2)
 se=sd(p^2)/sqrt(1000)
 c(est,se)

#########################################################
# Section 5.8 Rejection Sampling
#########################################################

betabinT=function(theta,datapar)
{
data=datapar$data
tpar=datapar$par
d=betabinexch(theta,data)-dmt(theta,mean=c(tpar$m),
  S=tpar$var,df=tpar$df,log=TRUE)
return(d)
}

tpar=list(m=fit$mode,var=2*fit$var,df=4)
datapar=list(data=cancermortality,par=tpar)

 start=c(-6.9,12.4)
 fit1=laplace(betabinT,start,datapar)
 fit1$mode

 betabinT(fit1$mode,datapar)

 theta=rejectsampling(betabinexch,tpar,-569.2813,10000,cancermortality)
 dim(theta)

 mycontour(betabinexch,c(-8,-4.5,3,16.5),cancermortality)
 points(theta[,1],theta[,2])

#############################################
# Section 5.9 Importance Sampling
#############################################

 tpar=list(m=fit$mode,var=2*fit$var,df=4)
 myfunc=function(theta)
   return(theta[2])
 s=impsampling(betabinexch,tpar,myfunc,10000,cancermortality)
 cbind(s$est,s$se)

##############################################
# Section 5.10 Sampling Importance Resampling
##############################################

tpar=list(m=fit$mode,var=2*fit$var,df=4)

theta.s=sir(betabinexch,tpar,10000,cancermortality)

 S=bayes.influence(theta.s,cancermortality)

 plot(c(0,0,0),S$summary,type="b",lwd=3,xlim=c(-1,21),
  ylim=c(5,11), xlab="Observation removed",ylab="log K")
 for (i in 1:20)
  lines(c(i,i,i),S$summary.obs[i,],type="b")

