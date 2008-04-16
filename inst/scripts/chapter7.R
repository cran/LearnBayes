# chapter7.tex

library(LearnBayes)


##############################################
# Section 7.3 Individual or Combined Estimates
##############################################

 data(hearttransplants)
 attach(hearttransplants)

 plot(log(e), y/e, pch = as.character(y), xlim=c(6,9.7))

##############################################
# Section 7.4 Equal Mortality Rates?
##############################################

sum(y)
sum(e)

 lambda=rgamma(1000,shape=277,rate=294681)
 ys94=rpois(1000,e[94]*lambda)

 hist(ys94,breaks=seq(1.5,26.5,by=1))
 lines(c(y[94],y[94]),c(0,120),lwd=3)

 pout=0*y
 lambda=rgamma(1000,shape=277,rate=294681)
 for (i in 1:94){
   ysi=rpois(1000,e[i]*lambda)
   pleft=sum(ysi<=y[i])/1000
   pright=sum(ysi>=y[i])/1000
   pout[i]=min(pleft,pright)
 }

 plot(log(e),pout,ylab="Prob(extreme)")

########################################################
# Section 7.5 Modeling a Prior Belief of Exchangeability
########################################################

 par(mfrow = c(2, 2))
 m = 500
 alphas = c(5, 20, 80, 400)
 for (j in 1:4) {
     mu = rgamma(m, shape = 10, rate = 10)
     lambda1 = rgamma(m, shape=alphas[j], rate=alphas[j]/mu)
     lambda2 = rgamma(m, shape=alphas[j], rate=alphas[j]/mu)
     plot(lambda1, lambda2)
     title(main=paste("alpha=",as.character(alphas[j])))
 }

#########################################################
# Section 7.7 Simulating from the Posterior
#########################################################

 datapar = list(data = hearttransplants, z0 = 0.53)
 start=c(2, -7)
 fit = laplace(poissgamexch, start, datapar)
 fit

 par(mfrow = c(1, 1))
 mycontour(poissgamexch, c(0, 8, -7.3, -6.6), datapar)
 title(xlab="log alpha",ylab="log mu")

 start = c(4, -7)
 fitgibbs = gibbs(poissgamexch, start, 1000, c(1,.15), datapar)
 fitgibbs$accept

 mycontour(poissgamexch, c(0, 8, -7.3, -6.6), datapar)
 points(fitgibbs$par[, 1], fitgibbs$par[, 2])

 plot(density(fitgibbs$par[, 1], bw = 0.2))

 alpha = exp(fitgibbs$par[, 1])
 mu = exp(fitgibbs$par[, 2])
 lam1 = rgamma(1000, y[1] + alpha, e[1] + alpha/mu)

 alpha = exp(fitgibbs$par[, 1])
 mu = exp(fitgibbs$par[, 2])
 plot(log(e), y/e, pch = as.character(y))
 for (i in 1:94) {
     lami = rgamma(1000, y[i] + alpha, e[i] + alpha/mu)
     probint = quantile(lami, c(0.05, 0.95))
     lines(log(e[i]) * c(1, 1), probint)
 }

 shrinkage = 0 * e
 for (i in 1:94) shrinkage[i] = mean(alpha/(alpha + e[i] * mu))

 plot(log(e), shrinkage)

 hospital=1:94
 meanrate=array(0,c(94,1))
 for (i in 1:94)
 meanrate[i]=mean(rgamma(1000, y[i] + alpha, e[i] + alpha/mu))
 hospital[meanrate==min(meanrate)]

 better=array(0,c(94,94))
 for (i in 1:94){
   for (j in (i+1):94){
   if (j <=94) {
   lami=rgamma(1000,y[i]+alpha,e[i]+alpha/mu)
   lamj=rgamma(1000,y[j]+alpha,e[j]+alpha/mu)
   better[i,j]=mean(lami<lamj)
   better[j,i]=1-better[i,j]
 }}}

 better[1:24,85]

#################################################
# Section 7.9 Posterior Predictive Model Checking
#################################################

 lam94=rgamma(1000,y[94]+alpha,e[94]+alpha/mu)

 ys24=rpois(1000,e[94]*lam94)

 hist(ys94,breaks=seq(1.5,39.5,by=1))
 lines(y[94]*c(1,1),c(0,100),lwd=3)

 pout.exchange=0*y
 for (i in 1:94){
   lami=rgamma(1000,y[i]+alpha,e[i]+alpha/mu)
   ysi=rpois(1000,e[i]*lami)
   pleft=sum(ysi<=y[i])/1000
   pright=sum(ysi>=y[i])/1000
   pout.exchange[i]=min(pleft,pright)
 }

 plot(pout,pout.exchange,xlab="P(extreme), equal means",
 ylab="P(extreme), exchangeable")
 abline(0,1)

