impsampling=function(logf,tpar,h,n,data)
{
theta=rmt(n,mean=c(tpar$m),S=tpar$var,df=tpar$df)

lf=logf(theta,data)       
lp=dmt(theta,mean=c(tpar$m),S=tpar$var,df=tpar$df,log=TRUE)
md=max(lf-lp)
wt=exp(lf-lp-md)
 
est=sum(wt*h(theta))/sum(wt)         
SEest=sqrt(sum((h(theta)-est)^2*wt^2))/sum(wt)  

return(list(est=est,se=SEest,theta=theta,wt=wt))
}
