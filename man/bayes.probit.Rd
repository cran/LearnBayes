\name{bayes.probit}
\alias{bayes.probit}
\title{Simulates from a probit binary response regression model using data augmentation and Gibbs sampling}
\description{
 Gives a simulated sample from the joint posterior distribution of the regression
vector for a binary response regression model with a probit link and a noninformative prior.
}
\usage{
bayes.probit(y,X,m)
}
\arguments{
  \item{y}{vector of binary responses}
  \item{X}{covariate matrix}
  \item{m}{number of simulations desired}
}
}
\value{
matrix of simulated draws of regression vector beta where each row corresponds to one draw
}
\author{Jim Albert}

\examples{
response=c(0,1,0,0,0,1,1,1,1,1)
covariate=c(1,2,3,4,5,6,7,8,9,10)
X=cbind(1,covariate)
m=1000
s=bayes.probit(response,X,m)
}

\keyword{models}
