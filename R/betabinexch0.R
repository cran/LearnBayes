betabinexch0=function(theta,data)
{
eta=theta[,1]                    
K=theta[,2]                     
y=data[,1]; n=data[,2]
N=length(y)
val=0*K;
for (i in 1:N)
   val=val+lbeta(K*eta+y[i],K*(1-eta)+n[i]-y[i])
val=val-N*lbeta(K*eta,K*(1-eta))
val=val-2*log(1+K)-log(eta)-log(1-eta)
return(val)
}

