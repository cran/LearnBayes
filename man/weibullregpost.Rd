\name{weibullregpost}
\alias{weibullregpost}
\title{Log posterior of a Weibull proportional odds model for survival data}
\description{
 Computes the log posterior density of (log sigma, mu, beta) for a Weibull
proportional odds regression model
}
\usage{
weibullregpost(theta,data)
}
\arguments{
  \item{theta}{matrix of parameter values where each row represents a value of (log sigma, mu, beta)}
  \item{data}{data matrix with columns survival time, censoring variable, and covariate matrix}
}
}
\value{
vector of values of the log posterior where each value corresponds to each row of the parameters in theta
}

\author{Jim Albert}

\examples{
data(chemotherapy)
attach(chemotherapy)
d=cbind(time,status,treat-1,age)
par1=c(-.6,11,.6,0)
par2=c(-.6,11,.6,.1)
theta=rbind(par1,par2)
weibullregpost(theta,d)
}

\keyword{models}
