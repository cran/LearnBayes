cauchyerrorpost=function(theta,data)
{
mu=theta[,1]; lsigma=theta[,2]
sigma=exp(lsigma)
val=0*mu
for (i in 1:length(data))
	{val=val+log(dt((data[i]-mu)/sigma,df=1)/sigma)}
return(val)
}

