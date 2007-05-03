poissgamexch=function(theta,datapar)
{
y=datapar$data[,2]; e=datapar$data[,1]
z0=datapar$z0
alpha=exp(theta[,1]); mu=exp(theta[,2])
beta=alpha/mu
N=length(y)
val=0*alpha;
for (i in 1:N)
{
val=val+lgamma(alpha+y[i])-(y[i]+alpha)*log(e[i]+beta)+alpha*log(beta)
}
val=val-N*lgamma(alpha)+log(alpha)-2*log(alpha+z0)
return(val)
}
