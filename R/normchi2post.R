normchi2post=function(theta,data)
{
mu=theta[,1]; sig2=theta[,2]
n=length(data)
z=0*mu
for (i in 1:n)
{
z=z-(data[i]-mu)^2/2/sig2
}
z=z-(n/2+1)*log(sig2)

return(z)
}
