lplot <- function(...) plot(..., type='l')
## stochastic Lorenz 96 system
## dx_i = ((x_{i+1} - x_{i-2}) * x_{i-1} - x_i + F)dt + sig_p * dB_t
## i=1, ..., N in circular connection (x_0=x_N)
StoLorenz96step <- function(x, dt, inc.t, sig_p, F) {
    ## simulate stochastic Lorenz 96 for dt
    ## assumes that dt is an integer multiple of inc.t
    ## x: current state; dt: simulation length; inc.t: discretization timestep
    ## sig_p: process noise; F: external forcing
    dim <- length(x)
    newx <- numeric(dim)
    for (n in 1:ceiling(dt/inc.t)) {
        newx[1] <- x[1] + ((x[2]-x[dim-1])*x[dim]-x[1]+F) * inc.t + sqrt(inc.t)*sig_p*rnorm(1)
        newx[2] <- x[2] + ((x[3]-x[dim])*x[1]-x[2]+F) * inc.t + sqrt(inc.t)*sig_p*rnorm(1)
        for (i in 3:(dim-1))
            newx[i] <- x[i] + ((x[i+1]-x[i-2])*x[i-1]-x[i]+F) * inc.t + sqrt(inc.t)*sig_p*rnorm(1)
        newx[dim] <- x[dim] + ((x[1]-x[dim-2])*x[dim-1]-x[dim]+F) * inc.t + sqrt(inc.t)*sig_p*rnorm(1)
        x <- newx
    }
    return(x)
}

StoLorenz96 <- function(x0, dt, inc.t, ndata, sig_p, F) {
    ## simulate stochastic Lorenz 96 for ndata points
    ## x0: initial point; other parameters are the same as defined in StoLorenz96step
    xmat <- x0
    x <- x0
    for (m in 1:(ndata-1)) {
        x <- StoLorenz96step(x, dt, inc.t, sig_p, F)
        xmat <- rbind(xmat, x)
    }
    return(xmat)
}

## simulate a sample path
##dim <- 4
##seed <- 827358
##set.seed(seed)
##ndata <- 25
##sig_p <- 1
##F <- 8
##x0 <- numeric(dim); x0[dim] <- .01
##x0 <- rnorm(dim,0,5)
##dt <- .1; inc.t <- .01
##sL <- replicate(10, StoLorenz96(x0=x0, dt=dt, inc.t=inc.t, ndata=ndata, sig_p=sig_p, F=F), simplify="array")
##nsL <- StoLorenz96(x0=x0, dt=dt, inc.t=inc.t, ndata=ndata, sig_p=0, F=F)
##par(mfrow=c(2,2), mar=c(2,2,1,1))
##for (co in 1:dim) {
##    plot(nsL[,co], ylim=c(-10,10), lty=1, col='red')
##    for (m in 1:10) { lines(sL[,co,m], lty=2) }
##}
