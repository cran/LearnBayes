weibullregpost=function(theta,data)
{
s=dim(data); k=s[2]; p=k-2
sp=dim(theta); N=sp[1]
t=data[,1]; c=data[,2]; X=data[,3:k]
sigma=exp(theta[,1])
mu=theta[,2]
beta=array(theta[,3:k],c(N,p))
n=length(t)
o=0*mu
for (i in 1:n)
{
  lp=0
  for (j in 1:p) lp=lp+beta[,j]*X[i,j]
  zi=(log(t[i])-mu-lp)/sigma
  fi=1/sigma*exp(zi-exp(zi))
  Si=exp(-exp(zi))
  o=o+c[i]*log(fi)+(1-c[i])*log(Si)
}
return(o)
}