\name{howardprior}
\alias{howardprior}
\title{Logarithm of Howard's dependent prior for two proportions}
\description{
Computes the logarithm of a dependent prior on two proportions
proposed by Howard in a Statistical Science paper in 1998.
}
\usage{
howardprior(xy,par)
}
\arguments{
  \item{xy}{matrix of parameter values where each row represents a value of the proportions (p1, p2)}
  \item{par}{vector containing parameter values alpha, beta, gamma, delta, sigma}
}
}
\value{
vector of values of the log posterior where each value corresponds to each row of the parameters in xy
}

\author{Jim Albert}

\examples{
par=c(1,1,1,1,2)
p1=array(c(.1, .5),c(1,2))
p2=array(c(.3, .8),c(1,2))
xy=cbind(p1,p2)
howardprior(xy,par)
}

\keyword{models}
