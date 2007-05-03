\name{lbinorm}
\alias{lbinorm}
\title{Logarithm of bivariate normal density}
\description{
Computes the logarithm of a bivariate normal density
}
\usage{
lbinorm(xy,par)
}
\arguments{
  \item{xy}{matrix of values where each row corresponds to a value of (x, y)}
  \item{par}{list with components m, a vector of means, and v, a variance-covariance matrix}
}
}
\value{
vector of values of the kernel of the log density
}

\author{Jim Albert}

\examples{
mean=c(0,0)
varcov=diag(c(1,1))
values=rbind(c(0,0),c(1,1))
par=list(m=mean,v=varcov)
lbinorm(values,par)
}

\keyword{models}
