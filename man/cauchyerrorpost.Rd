\name{cauchyerrorpost}
\alias{cauchyerrorpost}
\title{Log posterior of median and log scale parameters for Cauchy sampling}
\description{
 Computes the log posterior density of (M,log S)  when
a sample is taken from a Cauchy density with location M and scale S and a uniform
prior distribution is taken on (M, log S)
}
\usage{
cauchyerrorpost(theta,data)
}
\arguments{
  \item{theta}{matrix of parameter values where each row represents a value of (M, log S)}
  \item{data}{vector containing sample of observations}
}
}
\value{
vector of values of the log posterior where each value corresponds to each row of the parameters in theta
}

\author{Jim Albert}

\examples{
data=c(108, 51, 7, 43, 52, 54, 53, 49, 21, 48)
theta1=array(c(40,0),c(1,2))
theta2=array(c(50,1),c(1,2))
theta=rbind(theta1,theta2)
cauchyerrorpost(theta,data)
}

\keyword{models}
