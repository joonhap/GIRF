ntimes <- 200
ncity <- 50

args = commandArgs(trailingOnly=TRUE)
#args = c(2000, 0, 3, 1.1) # J, expno, B, var_infl_factor

F <- 8; sig_p <- 1; sig_m <- 1
J <- as.integer(args[1]); expno <- as.integer(args[2])
obs_int <- 0.5; obs_int_string <- "0.5"

state <- as.matrix(read.table(paste('../OBSINT_',obs_int_string,'/state_dim',ncity,'F',format(F,nsmall=6),'sig_p',format(sig_p,nsmall=6),'sig_m',format(sig_m,nsmall=6),'.txt', sep='')))
obs <- as.matrix(read.table(paste('../OBSINT_',obs_int_string,'/obs_dim',ncity,'F',format(F,nsmall=6),'sig_p',format(sig_p,nsmall=6),'sig_m',format(sig_m,nsmall=6),'.txt', sep='')))

state0 <- rep(0, ncity); state0[ncity] <- 0.01 # initial state




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
    filmean <- matrix(NA, ntimes, ncity)
    for (n in 1:ntimes) {
        x <- apply(x, 2, function(xx) StoLorenz96step(xx, obs_int, 0.01, sig_p, F))
        lla <- lla + log(mean(apply(x, 2, function(xx) dmvnorm(obs[n,], mean=xx, sigma=R))))
        C <- var(t(x)) # sample variance of x
        mu <- apply(x, 1, mean) # sample mean of x
        llb <- llb + dmvnorm(obs[n,], mean=mu, sigma=C+R, log=TRUE)
        D <- outer(obs[n,], rep(1,J)) + t(rmvnorm(n=J, sigma=R))
        x <- x + C %*% solve(C+R, D-x)
        filmean[n,] <- apply(x, 1, mean)
    }
    return(list(lla=lla, llb=llb, filmean=filmean))
}

##time.elapsed <- system.time({ result <- EnKF(J=J) })


EnKF_local_varinfl <- function(J, B, rho=1) {
    ##J: number of particles
    ##B: number of neighborhood components on one side (k-B:k+B is the neighborhood, where a ring-like structure is assumed)
    ##rho (>1): variance inflation parameter for the proposed particle distribution
    require('mvtnorm')
    x <- outer(state0, rep(1,J)) # initial particles
    R <- diag(rep(sig_m^2, ncity)) # observation error covariance
    ll <- 0 #log likelihood estimate
    filmean <- matrix(NA, ntimes, ncity)
    for (n in 1:ntimes) {
        x <- apply(x, 2, function(xx) StoLorenz96step(xx, obs_int, 0.01, sig_p, F))
        C <- var(t(x)) # sample variance of x
        mu <- apply(x, 1, mean) # sample mean of x
        x_infl <- sqrt(rho)*x - (sqrt(rho)-1)*outer(mu, rep(1,J))
        D <- outer(obs[n,], rep(1,J)) + t(rmvnorm(n=J, sigma=R))
        ll <- ll + dmvnorm(obs[n,], mean=mu, sigma=rho*C+R, log=TRUE)
        for (city in 1:ncity) {
            nbhd <- (city-B):(city+B)
            nbhd <- nbhd + (nbhd<=0) * ncity - (nbhd>ncity) * ncity # impose a circular structure
            ## ll <- ll + dmvnorm(obs[n,], mean=mu, sigma=C+R, log=TRUE)
            x[city,] <- x_infl[city,] + rho*C[city,nbhd] %*% solve(rho*C[nbhd,nbhd]+R[nbhd,nbhd], D[nbhd,]-x_infl[nbhd,])
        }
        filmean[n,] <- apply(x, 1, mean)
    }
    return(list(ll=ll, filmean=filmean))
}

B = as.integer(args[3]); rho = as.numeric(args[4])
filename <- paste0('EnKF_local_varinfl_K', ncity, '_F', format(round(F,6),nsmall=6), '_sig_p', format(round(sig_p,6),nsmall=6), 'J', J, 'T', ntimes, '_', expno, '_dtObs', obs_int, '_B', B, '_varinfl', format(round(rho,2),nsmall=2), '.RData')

lEnKF_results <- EnKF_local_varinfl(J=J, B=B, rho=rho)
ll <- lEnKF_results$ll

save(lEnKF_results, ll, file=filename)

##time.elapsed <- system.time({ result_local <- EnKF_local(J=J,B=3) })

#bb1 = EnKF_local(J=3000,B=3,rho=1) #-859.38
#bb2 = EnKF_local(J=3000,B=3,rho=1) #-857.71
#bb3 = EnKF_local(J=3000,B=3,rho=1) #-858.39
#bb4 = EnKF_local(J=3000,B=3,rho=1) #-858.96
#bbi1 = EnKF_local(J=3000,B=3,rho=1.1) #-860.87
#bbi2 = EnKF_local(J=3000,B=3,rho=1.1) #-859.86
#bbi3 = EnKF_local(J=3000,B=3,rho=1.1) #-862.14
#bbi4 = EnKF_local(J=3000,B=3,rho=1.1) #-860.38

#aa1 = EnKF(J=3000) #-855.1
#aa2 = EnKF(J=3000) #-852.54
#pomp EnKF results (J=3000): -852.37 -854.10 -853.19 -854.18
