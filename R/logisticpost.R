logisticpost=function(beta,data)
{
x=data[,1]; n=data[,2]; y=data[,3]
beta0=beta[,1]; beta1=beta[,2]

N=length(x)
z=0*beta0
for (i in 1:N)
{
lp=beta0+beta1*x[i]
pi=exp(lp)/(1+exp(lp))
z=z+y[i]*log(pi)+(n[i]-y[i])*log(1-pi)
}
return(z)
}
