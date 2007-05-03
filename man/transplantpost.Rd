\name{transplantpost}
\alias{transplantpost}
\title{Log posterior of a Pareto model for survival data}
\description{
 Computes the log posterior density of (log tau, log lambda, log p)  for a Pareto
model for survival data
}
\usage{
transplantpost(theta,data)
}
\arguments{
  \item{theta}{matrix of parameter values where each row represents a value of (log tau, log lambda, log p)}
  \item{data}{data matrix with columns survival time, transplant indicator, time to transplant, and
censoring indicator}
}
}
\value{
vector of values of the log posterior where each value corresponds to each row of the parameters in theta
}

\author{Jim Albert}

\examples{
data(stanfordheart)
par1=c(0,3,-1)
par2=c(0,2,0)
theta=rbind(par1,par2)
transplantpost(theta,stanfordheart)
}

\keyword{models}
