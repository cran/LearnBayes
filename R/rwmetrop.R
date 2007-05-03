rwmetrop=function(logpost,proposal,start,m,par)
{

pb=length(start)
Mpar=array(0,c(m,pb))                   # Sets up storage
b=t(start)
lb=logpost(start,par)      # Evaluates log likelihood at current value
a=chol(proposal$var)
scale=proposal$scale

accept=0;

for (i in 1:m)
{   
  bc=b+scale*t(a)%*%array(rnorm(pb),c(pb,1))    # Candidate value selected
  lbc=logpost(t(bc),par)  # Evaluates log likelihood at candidate value
  prob=exp(lbc-lb)           # Computes acceptance probability
  if (is.na(prob)==FALSE)
  {
  if (runif(1)<prob)
  {                
    lb=lbc; b=bc; accept=accept+1
  }
  }

  Mpar[i,]=b                # Store simulated value
}

accept=accept/m;               # Computes acceptance proportion

stuff=list(par=Mpar,accept=accept)
return(stuff)
}






