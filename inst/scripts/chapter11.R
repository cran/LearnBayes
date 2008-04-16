# chapter11.R

#################################################
# R script for Section 11.4 A Change-Point Model
#################################################

require(arm)

N=112
D=c(4,5,4,1,0,4,3,4,0,6,
3,3,4,0,2,6,3,3,5,4,5,3,1,4,4,1,5,5,3,4,2,5,2,2,3,4,2,1,3,2,
1,1,1,1,1,3,0,0,1,0,1,1,0,0,3,1,0,3,2,2,
0,1,1,1,0,1,0,1,0,0,0,2,1,0,0,0,1,1,0,2,
2,3,1,1,2,1,1,1,1,2,4,2,0,0,0,1,4,0,0,0,
1,0,0,0,0,0,1,0,0,1,0,0)
data=list("N","D")
inits = function() {list(b=c(0,0),changeyear=50)}
parameters <- c("changeyear","b")
coalmining.sim <- bugs (data, inits, parameters, "coalmining.bug", n.chains=3, n.iter=1000)
print(coalmining.sim)
attach.bugs(coalmining.sim)
plot(density(changeyear))

#####################################################

#######################################################
# R script for 11.5 A Robust Regression Model
#######################################################

require(arm)
require(LearnBayes)

data(election)
attach(election)
y=sqrt(buchanan)
x=sqrt(perot)
N=length(y)

data=list("N","y","x")
inits = function() {list(b=c(0,0),tau=1)}
parameters <- c("tau","lam","b")
robust.sim <- bugs (data, inits, parameters, "robust.bug", n.chains=3, n.iter=1000)
print(robust.sim)
attach.bugs(robust.sim)
mean(tau)
plot(tau)

#######################################################

#######################################################################
#  R script for Section 11.6  Estimating Career Trajectories
#######################################################################

require(arm)
require(LearnBayes)

data(sluggerdata)
s=careertraj.setup(sluggerdata)
N=s$N; T=s$T; y=s$y; n=s$n; x=s$x

mean = c(0, 0, 0)
Omega=diag(c(.1,.1,.1))
prec=diag(c(1.0E-6,1.0E-6,1.0E-6))

beta0=matrix(c(-7.69,.350,-.0058),nrow=10,ncol=3,byrow=TRUE)
mu.beta0=c(-7.69,.350,-.0058)
R0=diag(c(.1,.1,.1))

data=list("N","T","y","n","x","mean","Omega","prec")
inits = function() {list(beta=beta0,mu.beta=mu.beta0,R=R0)}
parameters <- c("beta")
career.sim <- bugs (data, inits, parameters, "career.bug", n.chains=1, n.iter=10000, n.thin=1)

###########################################################################
