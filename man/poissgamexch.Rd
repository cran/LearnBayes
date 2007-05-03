\name{poissgamexch}
\alias{poissgamexch}
\title{Log posterior of Poisson/gamma exchangeable model}
\description{
 Computes the log posterior density of log alpha and log mu for a Poisson/gamma exchangeable model
}
\usage{
poissgamexch(theta,datapar)
}
\arguments{
  \item{theta}{matrix of parameter values where each row represents a value of (log alpha, log mu)}
  \item{datapar}{list with components data, a matrix with columns e and y, and z0, prior hyperparameter}
}
}
\value{
vector of values of the log posterior where each value corresponds to each row of the parameters in theta
}

\author{Jim Albert}

\examples{
e=c(532,584,672,722,904)
y=c(0,0,2,1,1)
data=cbind(e,y)
par1=c(-4,0)
par2=c(-4,1)
z0=.5
datapar=list(data=data,z0=z0)
theta=rbind(par1,par2)
poissgamexch(theta,datapar)
}

\keyword{models}
