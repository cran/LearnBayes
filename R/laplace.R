laplace=function (logpost, mode, par) 
{
    {
        nr = function(f, old, par) {
            f1 = function(f, x, par) {
                h = 1e-04
                val = array(0, dim = c(1, length(x)))
                s = array(c(-h/2, h/2), dim = c(2, 1))
                x2 = array(1, dim = c(2, 1)) %*% x
                for (i in 1:length(x)) {
                  y = x2
                  y[, i] = y[, i] + s
                  v=(f(y[2,],par)-f(y[1,],par))/h
                  val[i] = v
                }
                return(val)
            }
            f2 = function(f, x, par) {
                h = 1e-04
                n = length(x)
                val = array(0, dim = c(n, n))
                s = array(c(-h, 0, h), dim = c(3, 1))
                x2 = array(1, dim = c(3, 1)) %*% x
                for (i in 1:n) {
                  y = x2
                  y[, i] = y[, i] + s
                  val[i, i]=(f(y[1,],par) - 2*f(y[2,],par) +
                           f(y[3,],par))/h^2
                }
                s = array(c(h/2, -h/2, h/2, -h/2, h/2, h/2, -h/2, 
                  -h/2), dim = c(4, 2))
                x2 = array(1, dim = c(4, 1)) %*% x
                if (n > 1) {
                  for (i in 1:(n - 1)) {
                    for (j in (i + 1):n) {
                      if (j <= n) {
                        y = x2
                        y[, c(i, j)] = y[, c(i, j)] + s
                        v=(f(y[1,],par)-f(y[2,],par)-f(y[3,],par)+
                          f(y[4,],par))/h^2
                        val[i, j] = v
                        val[j, i] = v
                      }
                    }
                  }
                }
                return(val)
            }
            h = -solve(f2(f, old, par))
            new = old + f1(f, old, par) %*% h
            p = length(old)
            int = p/2 * log(2 * pi) + 0.5 * log(det(h)) + f(new, 
                par)
            stuff = list(mode = new, var = h, int = int)
            return(stuff)
        }
    }
    stuff = nr(logpost, mode, par)
    converge = 0
    niter = 1
    mode.old = stuff$mode
    while (converge == 0 & niter < 10) {
        stuff = nr(logpost, stuff$mode, par)
        relative.error = sum(abs((stuff$mode - mode.old)/mode.old))
        if (relative.error < 1e-04) 
            converge = 1
        niter = niter + 1
        mode.old = stuff$mode
    }
    stuff$converge = as.logical(converge)
    stuff$iter = niter
    return(stuff)

}
