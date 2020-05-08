ntimes <- 200
ncity <- 50

args = commandArgs(trailingOnly=TRUE)
##args = c(2000, 0) # J, expno

F <- 8; sig_p <- 1; sig_m <- 1
J <- as.integer(args[1]); expno <- as.integer(args[2])
obs_int <- 0.5; obs_int_string <- "0.5"

state <- as.matrix(read.table(paste('../OBSINT_',obs_int_string,'/state_dim',ncity,'F',format(F,nsmall=6),'sig_p',format(sig_p,nsmall=6),'sig_m',format(sig_m,nsmall=6),'.txt', sep='')))
obs <- as.matrix(read.table(paste('../OBSINT_',obs_int_string,'/obs_dim',ncity,'F',format(F,nsmall=6),'sig_p',format(sig_p,nsmall=6),'sig_m',format(sig_m,nsmall=6),'.txt', sep='')))

state0 <- rep(0, ncity); state0[ncity] <- 0.01 # initial state



## Run EnKF using pomp with Csnippet
library(pomp)
identity <- function(x) { x }

obs.data.frame <- data.frame(time=obs_int*1:ntimes, Y=obs[1:ntimes,])
colnames(obs.data.frame) <- c("time", sprintf("Y%d", 1:ncity))

StoLorenz96step_pomp_C <- Csnippet("
    double *x = &X1;
    double newx[U];
    newx[0] = x[0] + ((x[1]-x[U-2])*x[U-1]-x[0]+F)*dt + sqrt(dt)*sig_p*rnorm(0,1);
    newx[1] = x[1] + ((x[2]-x[U-1])*x[0]-x[1]+F)*dt + sqrt(dt)*sig_p*rnorm(0,1);
    for (int i=3; i<U; i++) {
        newx[i-1] = x[i-1] + ((x[i]-x[i-3])*x[i-2]-x[i-1]+F)*dt + sqrt(dt)*sig_p*rnorm(0,1);
    }
    newx[U-1] = x[U-1] + ((x[0]-x[U-3])*x[U-2]-x[U-1]+F)*dt + sqrt(dt)*sig_p*rnorm(0,1);
    for (int i=0; i<U; i++) { x[i] = newx[i]; }
")

rinitC <- Csnippet("
    double *x = &X1;
    for (int i=0; i<U-1; i++) { x[i] = 0.0; }
    x[U-1] = 0.01;
")    

lorenz_globals <- Csnippet(paste0("#define U ", ncity, "\n"))

pomp.enkf.obj_C <- pomp(data=obs.data.frame,
                        t0=0,
                        times="time",
                        globals=lorenz_globals,
                        rinit=rinitC,
                        rprocess=euler(StoLorenz96step_pomp_C, delta.t=0.01),
                        obsnames=sprintf("Y%d",1:ncity),
                        statenames=sprintf("X%d",1:ncity),
                        paramnames=c("F", "sig_p")
                        )

set.seed(9879345+ 834581*expno)

time.elapsed <- system.time(EnKF.result_C <- enkf(data=pomp.enkf.obj_C, Np=J, h=identity, R=diag(rep(sig_m^2,ncity)), params=c(F=F, sig_p=sig_p)))

(ll <- logLik(EnKF.result_C))

filename <- paste0('EnKF_K', ncity, '_F', format(round(F,6),nsmall=6), '_sig_p', format(round(sig_p,6),nsmall=6), 'J', J, 'T', ntimes, '_', expno, '_dtObs', obs_int, '.RData') 

save(time.elapsed, EnKF.result_C, ll, file=filename)

##filter.mean(EnKF.result_C)








## run a separately defined EnKF
source('../Lorenz.R')

EnKF <- function(J) {
    ##J: number of particles
    require('mvtnorm')
    x <- outer(state0, rep(1,J)) # initial particles
    R <- diag(rep(sig_m^2, ncity)) # observation error covariance
    ## likelihood of data is estimated in two different ways
    ## a: hat{l}(y_n|y_{1:n-1}) = (1/J) sum_{j=1}^J phi(y_n | X_{t_n}^j, sigma_m^2 I) --> the average of the measurement density for each particle
    ## b: hat{l}(y_n|y_{1:n-1}) = phi(y_n | mean(X_{t_n}^j), var(X_{t_n}^j) + sigma_m^2 I) --> the mean and the variance of the particles are used.
    lla <- 0 #log likelihood estimate using method a.
    llb <- 0 #log likelihood estimate using method b.
    xmean <- c()
    for (n in 1:ntimes) {
        x <- apply(x, 2, function(xx) StoLorenz96step(xx, obs_int, 0.01, sig_p, F))
        lla <- lla + log(mean(apply(x, 2, function(xx) dmvnorm(obs[n,], mean=xx, sigma=R))))
        C <- var(t(x)) # sample variance of x
        mu <- apply(x, 1, mean) # sample mean of x
        llb <- llb + dmvnorm(obs[n,], mean=mu, sigma=C+R, log=TRUE)
        D <- outer(obs[n,], rep(1,J)) + t(rmvnorm(n=J, sigma=R))
        x <- x + C %*% solve(C+R, D-x)
        xmean <- rbind(xmean, apply(x, 1, mean)) 
    }
    return(list(lla=lla, llb=llb, xmean=xmean))
}

##time.elapsed <- system.time({ result <- EnKF(J) })

