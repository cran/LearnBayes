rwmetrop=function (logpost, proposal, start, m, par) 
{
    pb = length(start)
    Mpar = array(0, c(m, pb))
    b = matrix(t(start))
    lb = logpost(start, par)
    a = chol(proposal$var)
    scale = proposal$scale
    accept = 0
    for (i in 1:m) {
        bc = b + scale * t(a) %*% array(rnorm(pb), c(pb, 1))
        lbc = logpost(t(bc), par)
        prob = exp(lbc - lb)
        if (is.na(prob) == FALSE) {
            if (runif(1) < prob) {
                lb = lbc
                b = bc
                accept = accept + 1
            }
        }
        Mpar[i, ] = b
    }
    accept = accept/m
    stuff = list(par = Mpar, accept = accept)
    return(stuff)
}
