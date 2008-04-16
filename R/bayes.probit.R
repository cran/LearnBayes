bayes.probit=function (y, X, m, prior=list(beta=0,P=0)) 
{
   # probit regression with a N(beta0, P) prior on beta
   # P is the precision matrix

   rtruncated=function(n,lo,hi,pf,qf,...)
     qf(pf(lo,...)+runif(n)*(pf(hi,...)-pf(lo,...)),...)
   
    beta0=prior$beta; BI=prior$P
    N = length(y)
    beta = glm(y ~ X - 1, family = binomial(link = probit))$coef
    p = length(beta)
    beta = array(beta, c(p, 1))
    beta0=array(beta0, c(p, 1)); BI=array(BI, c(p, p))
    Mb = array(0, dim = c(m, p))
    lo=c(-Inf,0); hi=c(0,Inf); LO=lo[y+1]; HI=hi[y+1]; 
    aa = chol(solve(BI+t(X) %*% X)); BIbeta0=BI%*%beta0

    for (i in 1:m) {
	  z=rtruncated(N,LO,HI,pnorm,qnorm,X%*%beta,1)
        mn=solve(BI+t(X)%*%X,BIbeta0+t(X)%*%z)
        beta = t(aa) %*% array(rnorm(p), c(p, 1)) + mn
        Mb[i, ] = beta
    }
    return(Mb)
}
