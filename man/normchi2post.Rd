\name{normchi2post}
\alias{normchi2post}
\title{Log posterior density for mean and variance for normal sampling}
\description{
 Computes the log of the posterior density of a mean M and a variance S2 when a sample is taken from a normal density and a standard noninformative prior is used.
}
\usage{
normchi2post(theta,data)
}
\arguments{
  \item{theta}{matrix of parameter values where each row is a value of (M, S2)}
  \item{data}{vector containing the sample observations}
}
}
\value{
a vector of values of the log posterior where the values corresond to the rows in theta
}
\author{Jim Albert}

\examples{
parameter1=c(20,5)
parameter2=c(25,5)
parameter=rbind(parameter1,parameter2)
data=c(20, 32, 21, 43, 33, 21, 32)
normchi2post(parameter,data)
}

\keyword{models}
