transplantpost=function(theta,data)
{
x=data[,1]   #  survival time
y=data[,3]   #  time to transplant
t=data[,2]   #  transplant indicator
d=data[,4]   #  censoring indicator (d = 0 if died)

tau=exp(theta[,1])
lambda=exp(theta[,2])
p=exp(theta[,3])
val=0*tau

xnt=x[t==0]; dnt=d[t==0]
z=x[t==1]; y=y[t==1]; dt=d[t==1]
N=length(xnt)
M=length(z)

for (i in 1:N)
  val=val+(dnt[i]==0)*(p*log(lambda)+log(p)-(p+1)*log(lambda+xnt[i]))+
          (dnt[i]==1)*p*log(lambda/(lambda+xnt[i]))

for (i in 1:M)
  val=val+(dt[i]==0)*(p*log(lambda)+log(p*tau)-(p+1)*log(lambda+y[i]+tau*z[i]))+
          (dt[i]==1)*p*log(lambda/(lambda+y[i]+tau*z[i]))

val=val+theta[,1]+theta[,2]+theta[,3]
return(val)
}