\name{logisticpost}
\alias{logisticpost}
\title{Log posterior for a binary response model with a logistic link and a uniform prior}
\description{
 Computes the log posterior density of (beta0, beta1)  when
yi are independent binomial(ni, pi) and logit(pi)=beta0+beta1*xi and
a uniform prior is placed on (beta0, beta1)
}
\usage{
logisticpost(beta,data)
}
\arguments{
  \item{beta}{matrix of parameter values where each row represents a value of (beta0, beta1)}
  \item{data}{matrix of columns of covariate values x, sample sizes n, and number of successes y}
}
}
\value{
vector of values of the log posterior where each value corresponds to each row of the parameters in beta
}

\author{Jim Albert}

\examples{
x = c(-0.86,-0.3,-0.05,0.73)
n = c(5,5,5,5)
y = c(0,1,3,5)
data = cbind(x, n, y)
beta1=c(0,0)
beta2=c(2,10)
beta=rbind(beta1,beta2)
logisticpost(beta,data)
}

\keyword{models}
