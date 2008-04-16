normpostsim=function(data,m=1000)
{
S=sum((data-mean(data))^2)
xbar=mean(data)
n=length(data)
sigma2=S/rchisq(m,n-1)
mu=rnorm(m,mean=xbar,sd=sqrt(sigma2)/sqrt(n))
return(list(mu=mu, sigma2=sigma2))
}