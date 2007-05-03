groupeddatapost=function(theta,data)
{
cpoints=data$b
nbins=length(cpoints)
freq=data$f
m=theta[,1]; logsigma=theta[,2]

z=0*m; s=exp(logsigma)
z=freq[1]*log(pnorm(cpoints[1],m,s))
for (j in 1:(nbins-1))
  z=z+freq[j+1]*log(pnorm(cpoints[j+1],m,s)-pnorm(cpoints[j],m,s))
z=z+freq[nbins]*log(1-pnorm(cpoints[nbins],m,s))

return(z)
}