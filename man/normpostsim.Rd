\name{normpostsim}
\alias{normpostsim}
\title{Simulation from Bayesian normal sampling model}
\description{
 Gives a simulated sample from the joint posterior distribution of the mean and
variance for a normal sampling prior with a noninformative prior.
}
\usage{
normpostsim(data,m=1000)
}
\arguments{
  \item{data}{vector of observations}
  \item{m}{number of simulations desired}
}

\value{
\item{mu}{vector of simulated draws of normal mean}
\item{sigma2}{vector of simulated draws of normal variance}
}
\author{Jim Albert}

\examples{
data(darwin)
s=normpostsim(darwin$difference)
}

\keyword{models}
