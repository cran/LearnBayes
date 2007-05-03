\name{mycontour}
\alias{mycontour}
\title{Contour plot of a bivariate density function}
\description{
 For a general two parameter density, draws a contour graph where the
contour lines are drawn at 10 percent, 1 percent, and .1 percent of 
the height at the mode.
}
\usage{
mycontour(logf,limits,data)
}
\arguments{
  \item{logf}{function that defines the logarithm of the  density}
  \item{limits}{limits (xlo, xhi, ylo, yhi) where the graph is to be drawn}
  \item{data}{vector or list of parameters associated with the function logpost}
}
}
\value{
A contour graph of the density is drawn
}

\author{Jim Albert}

\examples{
lbinorm=function(xy,par)
{
	m=par$m
	v=par$v
	x=xy[,1]
	y=xy[,2]
	zx=(x-m[1])/sqrt(v[1,1])
	zy=(y-m[2])/sqrt(v[2,2])
	r=v[1,2]/sqrt(v[1,1]*v[2,2])
	return(-.5/(1-r^2)*(zx^2-2*r*zx*zy+zy^2))
}
m=array(c(0,0),c(2,1))
v=array(c(1,.6,.6,1),c(2,2))
normpar=list(m=m,v=v)
mycontour(lbinorm,c(-4,4,-4,4),normpar)
}

\keyword{models}
