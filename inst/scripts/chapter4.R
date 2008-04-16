# chapter4.R

library(LearnBayes)

######################################################
# Section 4.2 Normal Data with Both Parameters Unknown
######################################################

 data(marathontimes)
 attach(marathontimes)
 d = mycontour(normchi2post, c(220, 330, 500, 9000), time)
 title(xlab="mean",ylab="variance")

 S = sum((time - mean(time))^2) 
 n = length(time)
 sigma2 = S/rchisq(1000, n - 1)
 mu = rnorm(1000, mean = mean(time), sd = sqrt(sigma2)/sqrt(n))

 points(mu, sigma2)

 quantile(mu, c(0.025, 0.975))

 quantile(sqrt(sigma2), c(0.025, 0.975))

###################################################
# Section 4.3 A Multinomial Model
###################################################

 alpha = c(728, 584, 138)
 theta = rdirichlet(1000, alpha)

 hist(theta[, 1] - theta[, 2], main="")

###################################################
# Section 4.4 A Bioassay Experiment
###################################################

 x = c(-0.86, -0.3, -0.05, 0.73)
 n = c(5, 5, 5, 5)
 y = c(0, 1, 3, 5)
 data = cbind(x, n, y)

 glmdata = cbind(y, n - y)
 results = glm(glmdata ~ x, family = binomial)
 summary(results)

 mycontour(logisticpost, c(-4, 8, -5, 39), data)
 title(xlab="beta0",ylab="beta1")

 s = simcontour(logisticpost, c(-4, 8, -5, 39), data, 1000)
 points(s$x, s$y)

 plot(density(s$y),xlab="beta1")

 theta=-s$x/s$y
 hist(theta,xlab="LD-50")

 quantile(theta,c(.025,.975))

###########################################
# Section 4.5 Comparing Two Proportions
###########################################

 sigma=c(2,1,.5,.25)
 plo=.0001;phi=.9999
 par(mfrow=c(2,2))
 for (i in 1:4)
 {
 mycontour(howardprior,c(plo,phi,plo,phi),c(1,1,1,1,sigma[i]))
 title(main=paste("sigma=",as.character(sigma[i])),
   xlab="p1",ylab="p2")
 }

 sigma=c(2,1,.5,.25)
 par(mfrow=c(2,2))
 for (i in 1:4)
 {
 mycontour(howardprior,c(plo,phi,plo,phi),
   c(1+3,1+15,1+7,1+5,sigma[i]))
 lines(c(0,1),c(0,1))
 title(main=paste("sigma=",as.character(sigma[i])),
   xlab="p1",ylab="p2")
 }

 s=simcontour(howardprior,c(plo,phi,plo,phi),
   c(1+3,1+15,1+7,1+5,2),1000)
 sum(s$x>s$y)/1000

