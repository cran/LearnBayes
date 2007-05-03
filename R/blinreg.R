blinreg=function(y,X,m)
{
    fit = lm(y ~ 0 + X)
    bhat = matrix(fit$coef, c(1, fit$rank))
    s2 = sum(fit$residuals^2)/fit$df.residual
    shape = fit$df.residual/2
    rate = fit$df.residual/2 * s2
    sigma = sqrt(1/rgamma(m, shape = shape, rate = rate))
    vbeta = vcov(fit)/s2
    a = chol(vbeta)
    beta = array(1, c(m, 1)) %*% bhat + array(sigma, c(m, 1)) %*% 
        array(1, c(1, fit$rank)) * t(t(a) %*% matrix(rnorm(m * 
        fit$rank), c(fit$rank, m)))
    return(list(beta = beta, sigma = sigma))
}






