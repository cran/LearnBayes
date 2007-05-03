\name{betabinexch0}
\alias{betabinexch0}
\title{Log posterior of mean and precision for Binomial/beta exchangeable model}
\description{
 Computes the log posterior density of mean and precision for a Binomial/beta exchangeable model
}
\usage{
betabinexch0(theta,data)
}
\arguments{
  \item{theta}{matrix of parameter values where each row represents a value of (eta, K)}
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
par1=c(.1, 10)
par2=c(.2, 20)
theta=rbind(par1,par2)
betabinexch0(theta,data)
}

\keyword{models}
