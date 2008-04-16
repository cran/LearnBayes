#chapter8.R

library(LearnBayes)

###############################################
# Section 8.3 A One-Sided Test of a Normal Mean
###############################################

 pmean=170; pvar=25
 probH=pnorm(175,pmean,sqrt(pvar))
 probA=1-probH
 prior.odds=probH/probA
 prior.odds

 weights=c(182, 172, 173, 176, 176, 180, 173, 174, 179, 175)
 xbar=mean(weights)
 sigma2=3^2/length(weights)

 post.precision=1/sigma2+1/pvar
 post.var=1/post.precision

 post.mean=(xbar/sigma2+pmean/pvar)/post.precision
 c(post.mean,sqrt(post.var))

 post.odds=pnorm(175,post.mean,sqrt(post.var))/
  (1-pnorm(175,post.mean,sqrt(post.var)))
 post.odds

 BF = post.odds/prior.odds
 BF

 postH=probH*BF/(probH*BF+probA)
 postH

 z=sqrt(length(weights))*(mean(weights)-175)/3
 1-pnorm(z)

 weights=c(182, 172, 173, 176, 176, 180, 173, 174, 179, 175)
 data=c(mean(weights),length(weights),3)
 prior.par=c(170,1000)
 mnormt.onesided(175,prior.par,data)

#################################################
# Section 8.5  A Two-Sided Test of a Normal Mean
#################################################

 weights=c(182, 172, 173, 176, 176, 180, 173, 174, 179, 175)
 data=c(mean(weights),length(weights),3)
 t=c(.5,1,2,4,8)
 mnormt.twosided(170,.5,t,data)

#################################################
# Section 8.6 Models for Soccer Goals
#################################################

 data(soccergoals)
 attach(soccergoals)
 datapar=list(data=goals,par=c(4.57,1.43))
 fit1=laplace(logpoissgamma,.5,datapar)
 datapar=list(data=goals,par=c(1,.5))
 fit2=laplace(logpoissnormal,.5,datapar)
 datapar=list(data=goals,par=c(2,.5))
 fit3=laplace(logpoissnormal,.5,datapar)
 datapar=list(data=goals,par=c(1,2))
 fit4=laplace(logpoissnormal,.5,datapar)

 postmode=c(fit1$mode,fit2$mode,fit3$mode,fit4$mode)
 postsd=sqrt(c(fit1$var,fit2$var,fit3$var,fit4$var))
 logmarg=c(fit1$int,fit2$int,fit3$int,fit4$int)
 cbind(postmode,postsd,logmarg)

###################################################
# Section 8.7  Is a Baseball Hitter Really Streaky?
###################################################

 data(jeter2004)
 attach(jeter2004)
 data=cbind(H,AB)
 data1=regroup(data,5)

 logK=seq(2,6)
 logBF=0*logK
 for (j in 1:length(logK))
 {
 s=laplace(bfexch,0,list(data=data1,K=exp(logK[j])))
 logBF[j]=s$int
 }
 cbind(logK,exp(logK),logBF,exp(logBF))

###################################################################
# Section 8.8 A Test of Independence in a Two-Way Contingency Table
###################################################################

 data=matrix(c(11,9,68,23,3,5),c(2,3))
 data

 chisq.test(data)

 a=matrix(rep(1,6),c(2,3))
 a

 ctable(data,a)

 logK=seq(2,7,by=.2)
 logBF=0*logK
 for (j in 1:length(logK))
 {x=bfindep(data,exp(logK[j]),100000); logBF[j]=log(x$bf)}
 cbind(logK,logBF,exp(logBF))[seq(1,26,by=5),]

