bayes.probit=function(y,X,m)
{
N=length(y)
fit=glm(y~X-1,family=binomial(link=probit))  		   
beta=fit$coef
p=length(beta)
beta=array(beta,c(p,1))

Mb=array(0,dim=c(m,p))
aa=chol(solve(t(X)%*%X))

for (i in 1:m)
{
	lp=X%*%beta					           
	bb=pnorm(-lp)				            
	tt=(bb*(1-y)+(1-bb)*y)*runif(N)+bb*y	
	z=qnorm(tt)+lp			                  
						
	mn=solve(t(X)%*%X)%*%(t(X)%*%z)			                   
	beta=t(aa)%*%array(rnorm(p),c(p,1))+mn			            		
	Mb[i,]=t(beta) 
	
}
return(Mb)
}
