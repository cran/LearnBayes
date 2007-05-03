\name{logctablepost}
\alias{logctablepost}
\title{Log posterior of difference and sum of logits in a 2x2 table}
\description{
 Computes the log posterior density for the difference and sum of logits
in a 2x2 contingency table for independent binomial samples and uniform
prior placed on the logits
}
\usage{
logctablepost(theta,data)
}
\arguments{
  \item{theta}{matrix of parameter values where each row represents (difference of logits, sum of logits)}
  \item{data}{vector containing number of successes and failures for first sample, and then second sample}
}
}
\value{
vector of values of the log posterior where each value corresponds to each row of the parameters in theta
}

\author{Jim Albert}

\examples{
s1=6; f1=2; s2=3; f2=10
data=c(s1,f1,s2,f2)
theta1=c(2,4); theta2=c(1,1)
theta=rbind(theta1,theta2)
logctablepost(theta,data)
}

\keyword{models}
