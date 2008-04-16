# chapter3.R

library(LearnBayes)

######################################################################
# Section 3.2 Normal Distribution with Known Mean but Unknown Variance
######################################################################
 
 data(footballscores)
 attach(footballscores)
 d = favorite - underdog - spread
 n = length(d)
 v = sum(d^2)

 P = rchisq(1000, n)/v
 s = sqrt(1/P)
 hist(s)

 quantile(s, probs = c(0.025, 0.5, 0.975))

##########################################################
# Section 3.3 Estimating a Heart Transplant Mortality Rate
##########################################################

 alpha=16;beta=15174
 yobs=1; ex=66
 y=0:10
 lam=alpha/beta
 py=dpois(y, lam*ex)*dgamma(lam, shape = alpha,
   rate = beta)/dgamma(lam, shape= alpha + y,
   rate = beta + ex)
 cbind(y, round(py, 3))

 lambdaA = rgamma(1000, shape = alpha + yobs, rate = beta + ex)

 ex = 1767; yobs=4
 y = 0:10
 py = dpois(y, lam * ex) * dgamma(lam, shape = alpha, 
     rate = beta)/dgamma(lam, shape = alpha + y,
     rate = beta + ex)
 cbind(y, round(py, 3))

 lambdaB = rgamma(1000, shape = alpha + yobs, rate = beta + ex)

 lambda = seq(0, max(c(lambdaA,lambdaB)), length = 500)
 par(mfrow = c(2, 1))
 hist(lambdaA, freq = FALSE, main="", ylim=c(0,1500))
 lines(lambda, dgamma(lambda, shape = alpha, rate = beta))
 hist(lambdaB, freq = FALSE, main="", ylim=c(0,1500))
 lines(lambda, dgamma(lambda, shape = alpha, rate = beta))

####################################################
# Section 3.4 An Illustration of Bayesian Robustness
####################################################

 mu = 100
 tau = 12.16
 sigma = 15
 n = 4
 se = sigma/sqrt(4)
 ybar = c(110, 125, 140)
 tau1 = 1/sqrt(1/se^2 + 1/tau^2)
 mu1 = (ybar/se^2 + mu/tau^2) * tau1^2
 summ1=cbind(ybar, mu1, tau1)
 summ1

 tscale = 20/qt(0.95, 2)
 tscale

 theta = seq(60, 140, length = 200)
 plot(theta,1/tscale*dt((theta-mu)/tscale,2),
   type="l",ylab="Prior Density")
 lines(theta,1/10*dnorm((theta-mu)/tau),lwd=3)
 legend("topright",legend=c("t density","normal density"),
   lwd=c(1,3))

 summ2 = c()
 for (i in 1:3) {
     theta = seq(60, 180, length = 500)
     like = dnorm((theta - ybar[i])/7.5)
     prior = dt((theta - mu)/tscale, 2)
     post = prior * like
     post = post/sum(post)
     m = sum(theta * post)
     s = sqrt(sum(theta^2 * post) - m^2)
     summ2 = rbind(summ2, c(ybar[i], m, s))
 }
 summ2

 cbind(summ1,summ2)

 normpost = dnorm(theta, mu1[3], tau1)
 normpost = normpost/sum(normpost)
 plot(theta,normpost,type="l",lwd=3,ylab="Posterior Density")
 lines(theta,post)
 legend("topright",legend=c("t prior","normal prior"),lwd=c(1,3))

#######################################################
# Section 3.5 A Bayesian Test of the Fairness of a Coin
#######################################################

 pbinom(5, 20, 0.5)

 n = 20
 y = 5
 a = 10
 p = 0.5
 m1 = dbinom(y, n, p) * dbeta(p, a, a)/dbeta(p, a + y, a + n - 
     y)
 lambda = dbinom(y, n, p)/(dbinom(y, n, p) + m1)
 lambda

 pbetat(p,.5,c(a,a),c(y,n-y))

 loga = seq(-4, 5, length = 100)
 a = exp(loga)
 m2 = dbinom(y, n, p) * dbeta(p, a, a)/dbeta(p, a + y, a + n - 
     y)
 lambda = dbinom(y, n, p)/(dbinom(y, n, p) + m2)

par(mfrow=c(1,1))
plot(loga,lambda,type="l",xlab="log(a)",ylab="Prob(coin is fair)")

 n=20
 y=5
 a=10
 p=.5
 m2=0
 for (k in 0:y)
 {m2=m2+dbinom(k,n,p)*dbeta(p,a,a)/dbeta(p,a+k,a+n-k)}
 lambda=pbinom(y,n,p)/(pbinom(y,n,p)+m2)
 lambda

