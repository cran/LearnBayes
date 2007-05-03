rejectsampling=function (logf, tpar, dmax, n, data) 
{
    d=length(tpar$m)
    theta = rmt(n, mean = c(tpar$m), S = tpar$var, df = tpar$df)
    lf = logf(theta, data)
    lg = dmt(theta, mean = c(tpar$m), S = tpar$var, df = tpar$df, 
        log = TRUE)
    if (d==1)
    {
     prob = exp(c(lf) - lg - dmax)
     return(theta[runif(n)<prob])
    }
    else
    {
     prob = exp(lf - lg - dmax)
     return(theta[runif(n) < prob, ])
    }
}
