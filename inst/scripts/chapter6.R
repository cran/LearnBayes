# chapter6.R

library(LearnBayes)

####################################################
# Section 6.2 Introduction to Discrete Markov Chains
####################################################

 T=matrix(c(.5,.5,0,0,0,0,.25,.5,.25,0,0,0,0,.25,.5,.25,0,0,
           0,0,.25,.5,.25,0,0,0,0,.25,.5,.25,0,0,0,0,.5,.5),
           nrow=6,ncol=6,byrow=TRUE)
 T

 s=array(0,c(50000,1))

 s[1]=3
 for (j in 2:50000)
   s[j]=sample(1:6,size=1,prob=T[s[j-1],])

 m=c(500,2000,8000,50000)
 for (i in 1:4)
   print(table(s[1:m[i]])/m[i])

 w=matrix(c(.1,.2,.2,.2,.2,.1),nrow=1,ncol=6)
 w%*%T

##################################################################
# Section 6.7 Learning about a Normal Population from Grouped Data
##################################################################

LO=c(-Inf,66,68,70,72,74)
HI=c(66,68,70,72,74,Inf)
F=c(14,30,49,70,33,15)
d=list(int.lo=LO,int.hi=HI,f=F)

y=c(rep(65,14),rep(67,30),rep(69,49),rep(71,70),rep(73,33),
  rep(75,15))
 mean(y)

 log(sd(y))

 start=c(70,1)
 fit=laplace(groupeddatapost,start,d)
 fit

 modal.sds=sqrt(diag(fit$var))

 proposal=list(var=fit$var,scale=2)
 fit2=rwmetrop(groupeddatapost,proposal,start,10000,d)

 fit2$accept

 post.means=apply(fit2$par,2,mean)
 post.sds=apply(fit2$par,2,sd)

 cbind(c(fit$mode),modal.sds)

 cbind(post.means,post.sds)

 mycontour(groupeddatapost,c(69,71,.6,1.3),d)
 points(fit2$par[5001:10000,1],fit2$par[5001:10000,2])
 title(xlab="mu",ylab="log sigma")

##################################################
# Section 6.8 Example of Output Analysis
##################################################

 start=array(c(65,1),c(1,2))
 proposal=list(var=fit$var,scale=0.2)

 bayesfit=rwmetrop(groupeddatapost,proposal,start,10000,d)


###################################################
# Section 6.9 Modeling Data with Cauchy Errors
###################################################

 data(darwin)
 attach(darwin)
 mean(difference)

 log(sd(difference))

 laplace(cauchyerrorpost,c(21.6,3.6),difference)

 c(24.7-4*sqrt(34.96),24.7+4*sqrt(34.96))
 c(2.77-4*sqrt(.138),2.77+4*sqrt(.138))

 mycontour(cauchyerrorpost,c(-10,60,1,4.5),difference)
 title(xlab="mu",ylab="log sigma")

 fitlaplace=laplace(cauchyerrorpost,c(21.6,3.6),
   difference)
 mycontour(lbinorm,c(-10,60,1,4.5),list(m=fitlaplace$mode,
   v=fitlaplace$var))
 title(xlab="mu",ylab="log sigma")

 proposal=list(var=fitlaplace$var,scale=2.5)
 start=array(c(20,3),c(1,2))
 m=1000
 s=rwmetrop(cauchyerrorpost,proposal,start,m,difference)
 mycontour(cauchyerrorpost,c(-10,60,1,4.5),difference)
 title(xlab="mu",ylab="log sigma")
 points(s$par[,1],s$par[,2])

 fitgrid=simcontour(cauchyerrorpost,c(-10,60,1,4.5),difference,
  50000)
 proposal=list(var=fitlaplace$var,scale=2.5)
 start=array(c(20,3),c(1,2))
 fitrw=rwmetrop(cauchyerrorpost,proposal,start,50000,
  difference)
 proposal2=list(var=fitlaplace$var,mu=t(fitlaplace$mode))
 fitindep=indepmetrop(cauchyerrorpost,proposal2,start,50000,
  difference)
 fitgibbs=gibbs(cauchyerrorpost,start,50000,c(12,.75),
  difference)

 apply(fitrw$par,2,mean)

 apply(fitrw$par,2,sd)

#############################################################
# Section 6.10 Analysis of the Stanford Heart Transplant Data
#############################################################

data(stanfordheart)

 start=c(0,3,-1)
 laplacefit=laplace(transplantpost,start,stanfordheart)
 laplacefit

 proposal=list(var=laplacefit$var,scale=2)
 s=rwmetrop(transplantpost,proposal,start,10000,stanfordheart)
 s$accept

 tau=exp(s$par[,1])
 plot(density(tau),main="TAU")

 apply(exp(s$par),2,quantile,c(.05,.5,.95))

 p=exp(s$par[,1])
 lambda=exp(s$par[,2])
 t=seq(1,240)
 p5=0*t; p50=0*t; p95=0*t
 for (j in 1:240)
 { S=(lambda/(lambda+t[j]))^p
   q=quantile(S,c(.05,.5,.95))
   p5[j]=q[1]; p50[j]=q[2]; p95[j]=q[3]}
 plot(t,p50,type="l",ylim=c(0,1),ylab="Prob(Survival)",
   xlab="time")
 lines(t,p5,lty=2)
 lines(t,p95,lty=2)

