# chapter2.R

 library(LearnBayes)

####################################
# Section 2.3 Using a Discrete Prior
####################################

 p = seq(0.05, 0.95, by = 0.1)
 prior = c(2, 4, 8, 8, 4, 2, 1, 1, 1, 1)
 prior = prior/sum(prior)
 windows()
 plot(p, prior, type = "h", ylab="Prior Probability")

 data = c(11, 16)
 post = pdisc(p, prior, data)
 cbind(p, prior, post)
 
 windows()
 plot(p, post, type = "h", ylab="Posterior Probability")

################################
# Section 2.4 Using a Beta Prior
################################

 p = seq(0, 1, length = 500)
 a = 3.4
 b = 7.4
 s = 11
 f = 16
 prior=dbeta(p,a,b)
 like=dbeta(p,s+1,f+1)
 post=dbeta(p,a+s,b+f)
 windows()
 plot(p,post,type="l",ylab="Density",lty=2,lwd=3)
 lines(p,like,lty=1,lwd=3)
 lines(p,prior,lty=3,lwd=3)
 legend(.7,4,c("Prior","Likelihood","Posterior"),
     lty=c(3,1,2),lwd=c(3,3,3))

 1 - pbeta(0.5, a + s, b + f)

 qbeta(c(0.05, 0.95), a + s, b + f)

 ps = rbeta(1000, a + s, b + f)
 windows()
 hist(ps,xlab="p")

 sum(ps >= 0.5)/1000

 quantile(ps, c(0.05, 0.95))

#####################################
# Section 2.5 Using a Histogram Prior
#####################################

 midpt = seq(0.05, 0.95, by = 0.1)
 prior = c(2, 4, 8, 8, 4, 2, 1, 1, 1, 1)
 prior = prior/sum(prior)
 p = seq(0, 1, length = 500)

 plot(p,histprior(p,midpt,prior),type="l",
  ylab="Prior density",ylim=c(0,.25))

 like = dbeta(p, s + 1, f + 1)
 post = like * histprior(p, midpt, prior)

 plot(p, post, type = "l",ylab="Posterior density")

 post = post/sum(post)

 ps = sample(p, replace = TRUE, prob = post)

 hist(ps, xlab="p")

########################
# Section 2.6 Prediction
########################

 p=seq(0.05, 0.95, by=.1)
 prior=c(2, 4, 8, 8, 4, 2, 1, 1, 1, 1)
 prior=prior/sum(prior)
 m=20; ys=0:20
 pred=pdiscp(p, prior, m, ys)
 cbind(0:20,pred)

 ab=c(3.4, 7.4)
 m=20; ys=0:20
 pred=pbetap(ab, m, ys)

 p=rbeta(1000,3.4, 7.4)

 y = rbinom(1000, 20, p)

 table(y)

 freq=table(y)
 ys=ys=as.integer(names(freq))
 predprob=freq/sum(freq)
 plot(ys,predprob,type="h",xlab="y",
   ylab="Predictive Probability")

 dist=cbind(ys,predprob)

 covprob=.9
 discint(dist,covprob)

