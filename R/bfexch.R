bfexch=function(theta,datapar)
{
y=datapar$data[,1]; n=datapar$data[,2]; K=datapar$K
eta=exp(theta)/(1+exp(theta))
N=length(y)
z=0*theta;
for (i in 1:N)
  z=z+lbeta(K*eta+y[i],K*(1-eta)+n[i]-y[i])
z=z-N*lbeta(K*eta,K*(1-eta))+log(eta*(1-eta))
z=z-lbeta(sum(y)+1,sum(n-y)+1)
return(z)
}