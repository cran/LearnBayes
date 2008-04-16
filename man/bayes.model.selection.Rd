\name{bayes.model.selection}
\alias{bayes.model.selection}
\title{Bayesian regression model selection using G priors}
\description{
Using Zellner's G priors, computes the log marginal density for all possible regression models
}
\usage{
bayes.model.selection(y, X, c, constant=TRUE)
}
\arguments{
  \item{y}{vector of response values}
  \item{X}{matrix of covariates}
  \item{c}{parameter of the G prior}
  \item{constant}{logical variable indicating if a constant term is in the matrix X}
}

\value{
\item{mod.prob}{matrix specifying the model and the associated value of the log marginal density}
\item{converge}{logical vector indicating if the laplace algorithm converged for each model}
}

\author{Jim Albert}

\examples{
data(birdextinct)
logtime=log(birdextinct$time)
X=cbind(1,birdextinct$nesting,birdextinct$size,birdextinct$status)
bayes.model.selection(logtime,X,100)
}

\keyword{models}
