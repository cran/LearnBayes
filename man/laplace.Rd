\name{laplace}
\alias{laplace}
\title{Summarization of a posterior density by the Laplace method}
\description{
 For a general posterior density, computes the posterior mode, 
the associated variance-covariance matrix, and an estimate at the
logarithm at the normalizing constant.
}
\usage{
laplace(logpost,mode,iter,par)
}
\arguments{
  \item{logpost}{function that defines the logarithm of the posterior density}
  \item{mode}{vector that is a guess at the posterior mode}
  \item{iter}{number of iterations of Newton-Raphson algorithm}
  \item{par}{vector or list of parameters associated with the function logpost}
}
}
\value{
\item{mode}{current estimate at the posterior mode}
\item{var}{current estimate at the associated variance-covariance matrix}
\item{int}{estimate at the logarithm of the normalizing constant}
}

\author{Jim Albert}

\examples{
logpost=function(theta,data)
{
s=5
val=0*theta
for (i in 1:length(data))
{
val=val-log(1+(data[i]-theta)^2/s^2)
}
return(val)
}
data=c(10,12,14,13,12,15)
laplace(logpost,10,5,data)
}

\keyword{models}
