laplace=function(logpost,mode,iter,par)
#------------------------
# Written by Jim Albert
# albert@bgnet.bgsu.edu
# November 2004
#------------------------
{
{
nr=function(f,old,par)
{
f1=function(f,x,par)
{
h=.0001
val=array(0,dim=c(1,length(x)))
s=array(c(-h/2,h/2),dim=c(2,1))
x2=array(1,dim=c(2,1))%*%x
for (i in 1:length(x))
{
   y=x2; y[,i]=y[,i]+s
   v=diff(f(y,par))/h
   val[i]=v
}
return(val)
}
f2=function(f,x,par)
{
h=.0001; n=length(x)
val=array(0,dim=c(n,n))
s=array(c(-h,0,h),dim=c(3,1))
x2=array(1,dim=c(3,1))%*%x
for (i in 1:n)
{
   y=x2; y[,i]=y[,i]+s
   t=f(y,par)
   val[i,i]=(t[1]-2*t[2]+t[3])/h^2
}
s=array(c(h/2,-h/2,h/2,-h/2,h/2,h/2,-h/2,-h/2),dim=c(4,2))
x2=array(1,dim=c(4,1))%*%x
if (n>1)
{
for (i in 1:(n-1))
{
   for (j in (i+1):n)
   {  
      if (j<=n)
      {
	y=x2;y[,c(i,j)]=y[,c(i,j)]+s
      t=f(y,par)
	v=(t[1]-t[2]-t[3]+t[4])/h^2
      val[i,j]=v; val[j,i]=v
      }
   }
}
}
return(val)
}

h=-solve(f2(f,old,par))
new=old+f1(f,old,par)%*%h
p=length(old)
int=p/2*log(2*pi)+.5*log(det(h))+f(new,par)
stuff=list(mode=new,var=h,int=int)
return(stuff)
}


}
stuff=nr(logpost,mode,par)
for (i in 1:iter)
 {
  stuff=nr(logpost,stuff$mode,par)
  }
return(stuff)
}
