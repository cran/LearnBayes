\name{betabinexch}
\alias{betabinexch}
\title{Log posterior of logit mean and log precision for Binomial/beta exchangeable model}
\description{
 Computes the log posterior density of logit mean and log precision for a Binomial/beta exchangeable model
}
\usage{
betabinexch(theta,data)
}
\arguments{
  \item{theta}{matrix of parameter values where each row represents a value of (logit eta, log K)}
  \item{data}{a matrix with columns y (counts) and n (sample sizes)}
}
}
\value{
vector of values of the log posterior where each value corresponds to each row of the parameters in theta
}

\author{Jim Albert}

\examples{
n=c(20,20,20,20,20)
y=c(1,4,3,6,10)
data=cbind(y,n)
par1=c(-1,0)
par2=c(-1,1)
theta=rbind(par1,par2)
betabinexch(theta,data)
}

\keyword{models}
