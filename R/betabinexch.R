betabinexch=function(theta,data)
{
theta1=theta[,1]                     
theta2=theta[,2]                     
eta=exp(theta1)/(1+exp(theta1))
K=exp(theta2)
y=data[,1]; n=data[,2]      
N=length(y);                    
val=0*K;
for (i in 1:N)
   val=val+lbeta(K*eta+y[i],K*(1-eta)+n[i]-y[i])
val=val-N*lbeta(K*eta,K*(1-eta))
val=val+theta2-2*log(1+exp(theta2))
return(val)
}

