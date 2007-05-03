\name{groupeddatapost}
\alias{groupeddatapost}
\title{Log posterior of normal parameters when data is in grouped form}
\description{
 Computes the log posterior density of (M,log S) for normal sampling where the data is
observed in grouped form
}
\usage{
groupeddatapost(theta,data)
}
\arguments{
  \item{theta}{matrix of parameter values where each row represents a value of (M, log S)}
  \item{data}{list with components b, a vector of midpoints and f, the corresponding bin frequencies}
}
}
\value{
vector of values of the log posterior where each value corresponds to each row of the parameters in theta
}

\author{Jim Albert}

\examples{
b=c(10,15,20,25,30)
f=c(2,5,8,4,2)
data=list(b=b,f=f)
par1=c(20,1)
par2=c(22,0)
theta=rbind(par1,par2)
groupeddatapost(theta,data)
}

\keyword{models}
