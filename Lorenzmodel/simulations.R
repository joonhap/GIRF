ncity <- 50

F <- 8; sig_p <- 1; sig_m <- 1
J <- 2000; 
obs_int <- 0.5; obs_int_string <- "0.5"

source('../Lorenz.R')

state0 <- rep(0, ncity); state0[ncity] <- 0.01 # initial state

### in order to assess the degree to which the one-step forecast distribution (the distribution of X_{t_{n+1}} given fixed X_{t_n}=x) is non-Gaussian, simulate a number of stochastic Lorenz paths
## simulate a bunch of times
set.seed(123120)
x0 <- StoLorenz96(x0=state0, dt=100, inc.t=0.01, ndata=2, sig_p=sig_p, F=F)[2,]
x0.1 <- replicate(n=2000, StoLorenz96(x0=x0, dt=0.1, inc.t=0.01, ndata=2, sig_p=sig_p, F=F)[2,])
x0.5 <- replicate(n=2000, StoLorenz96(x0=x0, dt=0.5, inc.t=0.01, ndata=2, sig_p=sig_p, F=F)[2,])
x1.0 <- replicate(n=2000, StoLorenz96(x0=x0, dt=1.0, inc.t=0.01, ndata=2, sig_p=sig_p, F=F)[2,])
x2.0 <- replicate(n=2000, StoLorenz96(x0=x0, dt=2.0, inc.t=0.01, ndata=2, sig_p=sig_p, F=F)[2,])

## deterministic simulations
x0.1det <- StoLorenz96(x0, dt=0.1, inc.t=0.01, ndata=2, sig_p=0, F=F)[2,]
x0.5det <- StoLorenz96(x0, dt=0.5, inc.t=0.01, ndata=2, sig_p=0, F=F)[2,]
x1.0det <- StoLorenz96(x0, dt=1.0, inc.t=0.01, ndata=2, sig_p=0, F=F)[2,]
x2.0det <- StoLorenz96(x0, dt=2.0, inc.t=0.01, ndata=2, sig_p=0, F=F)[2,]

## dt=0.5
par(mfrow=c(3,3),mar=c(0,0,0,0))
for (component in c(1:8,50)) { 
    hist(x0.5[component,]); abline(v=x0.5det[component], col='red')
}

par(mfrow=c(1,1),mar=c(2,2,0,0))
comps <- c(26,50)
plot(x0.5[comps[1],],x0.5[comps[2],]); points(x0.5det[comps[1]],x0.5det[comps[2]],pch='x',cex=2,col='red')

## dt=1.0
par(mfrow=c(3,3),mar=c(0,0,0,0))
for (component in c(1:4,46:50)) {
    hist(x1.0[component,]); abline(v=x1.0det[component], col='red')
}

par(mfrow=c(1,1),mar=c(2,2,0,0))
comps <- c(1,2)
plot(x1.0[comps[1],],x1.0[comps[2],]); points(x1.0det[comps[1]],x1.0det[comps[2]],pch='x',cex=2,col='red')

## dt=0.1
par(mfrow=c(3,3),mar=c(0,0,0,0))
for (component in c(1:8,50)) {
    hist(x0.1[component,]); abline(v=x0.1det[component], col='red')
}

par(mfrow=c(1,1),mar=c(2,2,0,0))
comps <- c(1,2)
plot(x0.1[comps[1],],x0.1[comps[2],]); points(x0.1det[comps[1]],x0.1det[comps[2]],pch='x',cex=2,col='red')

## dt=2.0
par(mfrow=c(3,3),mar=c(0,0,0,0))
for (component in c(1:4,46:50)) {
    hist(x2.0[component,]); abline(v=x2.0det[component], col='red')
}

par(mfrow=c(1,1),mar=c(2,2,0,0))
comps <- c(23,25)
plot(x2.0[comps[1],],x2.0[comps[2],]); points(x2.0det[comps[1]],x2.0det[comps[2]],pch='x',cex=2,col='red')

