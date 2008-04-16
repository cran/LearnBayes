# chapter9.R

library(LearnBayes)

################################
# Section 9.2.6 An Example
################################

 data(birdextinct)
 attach(birdextinct)
 logtime=log(time)
 plot(nesting,logtime)
 identify(nesting,logtime,label=species,n=5)
 plot(jitter(size),logtime,xaxp=c(0,1,1))
 plot(jitter(status),logtime,xaxp=c(0,1,1))

##### Least-squares fit

 fit=lm(logtime~nesting+size+status,data=birdextinct,x=TRUE,y=TRUE)
 summary(fit)

##### Sampling from posterior

 theta.sample=blinreg(fit$y,fit$x,5000)

 par(mfrow=c(2,2))
 hist(theta.sample$beta[,2],main="NESTING",
  xlab=expression(beta[1]))
 hist(theta.sample$beta[,3],main="SIZE",
  xlab=expression(beta[2]))
 hist(theta.sample$beta[,4],main="STATUS",
  xlab=expression(beta[3]))
 hist(theta.sample$sigma,main="ERROR SD",
  xlab=expression(sigma))

 apply(theta.sample$beta,2,quantile,c(.05,.5,.95))

 quantile(theta.sample$sigma,c(.05,.5,.95))

###### Estimating mean extinction times

 cov1=c(1,4,0,0)
 cov2=c(1,4,1,0)
 cov3=c(1,4,0,1)
 cov4=c(1,4,1,1)
 X1=rbind(cov1,cov2,cov3,cov4)
 mean.draws=blinregexpected(X1,theta.sample)
 par(mfrow=c(2,2))
 hist(mean.draws[,1],main="Covariate set A",xlab="log TIME")
 hist(mean.draws[,2],main="Covariate set B",xlab="log TIME")
 hist(mean.draws[,3],main="Covariate set C",xlab="log TIME")
 hist(mean.draws[,4],main="Covariate set D",xlab="log TIME")

######## Predicting extinction times

 cov1=c(1,4,0,0)
 cov2=c(1,4,1,0)
 cov3=c(1,4,0,1)
 cov4=c(1,4,1,1)
 X1=rbind(cov1,cov2,cov3,cov4)
 pred.draws=blinregpred(X1,theta.sample)
 par(mfrow=c(2,2))
 hist(pred.draws[,1],main="Covariate set A",xlab="log TIME")
 hist(pred.draws[,2],main="Covariate set B",xlab="log TIME")
 hist(pred.draws[,3],main="Covariate set C",xlab="log TIME")
 hist(pred.draws[,4],main="Covariate set D",xlab="log TIME")

######### Model checking via posterior predictive distribution

 pred.draws=blinregpred(fit$x,theta.sample)
 pred.sum=apply(pred.draws,2,quantile,c(.05,.95))
 par(mfrow=c(1,1))
 ind=1:length(logtime)
 matplot(rbind(ind,ind),pred.sum,type="l",lty=1,col=1,
  xlab="INDEX",ylab="log TIME")
 points(ind,logtime,pch=19)

######### Model checking via bayes residuals

 prob.out=bayesresiduals(fit,theta.sample,2)
 par(mfrow=c(1,1))
 plot(nesting,prob.out)
 identify(nesting,prob.out,label=species,n=4)

##############################################
# Section 9.3 Survival Modeling
##############################################

 data(chemotherapy)
 attach(chemotherapy)
 library(survival)
 survreg(Surv(time,status)~factor(treat)+age,dist="weibull")

 start=c(-.5,9,.5,-.05)
 d=cbind(time,status,treat-1,age)
 fit=laplace(weibullregpost,start,d)
 fit

 proposal=list(var=fit$var,scale=1.5)
 bayesfit=rwmetrop(weibullregpost,proposal,fit$mode,10000,d)
 bayesfit$accept

 par(mfrow=c(2,2))
 sigma=exp(bayesfit$par[,1])
 mu=bayesfit$par[,2]
 beta1=bayesfit$par[,3]
 beta2=bayesfit$par[,4]
 hist(beta1,xlab="treatment",main="")
 hist(beta2,xlab="age",main="")
 hist(sigma,xlab="sigma",main="")

