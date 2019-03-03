ntimes <- 415

ncity <- 40
measles <- read.csv('../../UKmeasles.csv')
citynames <- colnames(measles)[-1]

case <- TRUE
R <- 5; J <- 600; S <- ncity; runif <- TRUE; weighted <- TRUE; max_la <- 3; n <- 0; M <- ifelse(runif, 10, 1)

lth <- function(G, n, bo=0) { ## bo: burn-out: as the initial condition is not known, the likelihood estimate of the first data points will likely have higher variability. In order to obtain more reliable comparison between different G values, we sum the likelihood estimates only for the time steps after certain point.
    filenamecore <- paste('_K', ncity, ifelse(case, '_case_', ''), 'G', format(round(G,6),nsmall=6), ifelse(case, '', paste('rep', format(round(0.5,6),nsmall=6), sep='')), 'R', R, 'J', J, 'S', S, 'T', ntimes, 'if_', ifelse(runif,'T','F'), '_wi_', ifelse(weighted,'T','F'), '_la', max_la, '_', n, sep='') 
    ## compute likelihood
    lfilename <- paste('im_est_likelihood', filenamecore, '.txt', sep='')
    if(file.exists(lfilename)) { lest <- as.matrix(read.table(lfilename))
    } else { print(paste("Warning: no lest file exits.", "G =",G, 'n =', n)); lest <- NA }
    ll <- apply(cbind(matrix(c(lest), nrow=ntimes)[(bo+1):ntimes,]), 2, mean, na.rm=FALSE) * (ntimes-bo) # ll for each m
    nNaN <- apply(matrix(c(lest), nrow=ntimes), 2, function(x) sum(is.na(x)))
    if ((lenll <- length(ll)) < M) { ll[(lenll+1):M] <- NA }
    if ((lennn <- length(nNaN)) < M) { nNaN[(lennn+1):M] <- NA }
    if (runif) {
        thfilename <- paste('im_theta', filenamecore, '.txt', sep='')
        if(file.exists(thfilename)) { th <- read.table(thfilename)
        } else { print(paste('Warning: no theta file exists.', "G =", G, 'n =', n)); th <- matrix(NA, ntimes*M, thsize) }
        if ((lenth <- dim(th)[1]) < M*ntimes) { th[(lenth+1):(M*ntimes),] <- NA }
        return(list(ll=ll, nNaN=nNaN, th=th))
    }
    return(list(ll=ll, nNaN=nNaN))
}


G <- 100; n <- 0; R <- 5; J <- 600; S <- ncity; runif <- TRUE; weighted <- TRUE; max_la <- 3
filenamecore <- paste('_K', ncity, '_case_', 'G', format(round(G,6),nsmall=6), 'R', R, 'J', J, 'S', S, 'T', ntimes, 'if_', ifelse(runif,'T','F'), '_wi_', ifelse(weighted,'T','F'), '_la', max_la, '_', n, sep='') 


## ivp estimation
ivplfilename <- paste('im_ivp_est_likelihood', filenamecore, '.txt', sep='')
ivplraw <- unlist(read.table(ivplfilename, sep=' '))
ivpl <- apply(matrix(ivplraw, 3, 80), 2, sum)
plot(ivpl)

ivpfilename <- paste('im_ivp', filenamecore, '.txt', sep='')
ivpraw <- read.table(ivpfilename, sep=' ')

ivpthetaname <- paste('ivp_thetaswarm_evallik', filenamecore, '.txt', sep='')
logistic <- function(x) { 1/(1+exp(-x)) }
ivpthetaraw <- logistic(read.table(ivpthetaname, sep=' '))[,18:137]
compsd <- t(matrix(apply(ivpthetaraw, 2, sd),3,ncity))
plot(compsd[,1])

city <- 13; SEI <- 1; comp <- (city-1)*3+SEI
plot(ivpthetaraw[, comp])
abline(h=ivpraw[240, comp])
abline(h=ivpraw[240,comp]-compsd[city, SEI], lty=2, col='red')
abline(h=ivpraw[240,comp]+compsd[city, SEI], lty=2, col='red')

## read estimated state means
sfilename <- paste('im_est_state_mean', filenamecore, '.txt', sep='')
ss <- read.table(sfilename, sep=' ')
city <- 22
startyear <- 1949
startrow <- which.max(measles[,1]>startyear)
plot(ss[,5*(city-1)+4], type='l', ylab='estimated biweekly counts')
lines(measles[-(1:(startrow-1)),city+1], lty=2, col='red')




## plot theta ##
thfilename <- paste('im_theta', filenamecore, '.txt', sep='')
theta <- read.table(thfilename, sep=' ')
logitamp <- function(x) { .68/(exp(-x)+1) }
iden <- function(x) { x }
logit <- function(x) { 1/(exp(-x)+1) }
trfun <- c(exp, logitamp, iden, iden, exp, logit, exp, logit, exp, exp, iden, iden, iden, logit)
parname <- c("R0", "amp", "alpha", "mu", "gen", "infPerProp", "sigma2", "rep", "repOD", "G", "sourcePow", "destPow", "distPow", "cohortEntry")
param <- sapply(1:length(trfun), function(j) { sapply(theta[,j], trfun[[j]]) } )
par(mfcol=c(7,2), mar=c(0,3,0,1), oma=c(1,1,1,0))
for (j in 1:length(trfun)) {
    plot(param[,j], type='l', xlab='', ylab='')
    mtext(side=2, line=2, text=parname[j])
}

## plot Monte Carlo MLE
Gvalues <- c(20,50,100,200,300,400)
Gn <- cbind(rep(Gvalues, each=6), rep(0:5,length(Gvalues)))
runif <- FALSE; R <- 10; J <- 1000; M <- 1; ll <- matrix(sapply(1:dim(Gn)[1], function(x) { do.call(lth, as.list(Gn[x,]))$ll }), nrow=6) # likelihood estimates
top3 <- apply(ll, 2, function(x) { which(rank(x) > 3) }) # pick highest 3 likelihood estimates
top3 <- c(top3) + 6*rep(0:(length(Gvalues)-1), each=3)
lltop3 <- c(ll)[top3]
runif <- TRUE; R <- 5; J <- 600; M <- 10; MCMLE <- sapply(1:dim(Gn)[1], function(x) { th=do.call(lth, as.list(Gn[x,]))$th; th[dim(th)[1],] }) # likelihood estimates
MCMLE <- t(MCMLE[1:length(trfun),]) # take only regular (non-IVP) parameters
MCMLEtop3 <- MCMLE[top3,] # pick points with top 3 likelihood estimates
MCMLEtop3 <- sapply(1:length(trfun), function(j) { sapply(MCMLEtop3[,j], trfun[[j]]) } ) # transform parameters into natural scale
library('latex2exp')

pdf('plots/pairs_MCMLE_ll_top3.pdf')
pairs(cbind(lltop3, MCMLEtop3[,-c(4,8,11,12,13)]), labels=c(TeX('$\\hat{l}$'), TeX('$R_0$'), 'a', TeX('$\\alpha$'), TeX('$\\nu_{EI}^{-1} + \\nu_{IR}^{-1}$'), TeX('$\\frac{\\nu_{IR}^{-1}}{\\nu_{EI}^{-1} + \\nu_{IR}^{-1}}$'), TeX('$\\sigma^2$'), TeX('$\\psi$'), 'G', 'c'), cex=.6, cex.labels=1 ) # pairs plot
dev.off()

## plot likelihood estimates ##
lth100_0 <- lth(100,0)$ll

## plot likelihood estimates for ivp estimation ##
ivplfilename <- paste('im_ivp_est_likelihood', filenamecore, '.txt', sep='')
ivpl <- apply(matrix(c(as.matrix(read.table(ivplfilename, sep=' '))), 3, 300), 2, sum)
plot(ivpl, type='l')
abline(v=50*1:6, lty=2)

## IF2 plot (beta, G, dist_pow)
##pdf(paste('plots/', date, '/IF2', filenamecore, '.pdf', sep=''))

expand.grid2 <- function(...) { al <- list(...); y <- expand.grid(rev(al)); retval <- y[,dim(y)[2]:1]; colnames(retval) <- NULL; return(retval) }

nvec <- c(0:5); gvec <- c(20, 50,100*c(1,2,3,4)); evalid <- c(0)

param <- expand.grid2(gvec,nvec)
thList <- list(); ll <- matrix(NA,nrow=M,ncol=dim(param)[1]); nNaN <- matrix(NA,nrow=M,ncol=dim(param)[1])
for (i in 1:dim(param)[1]) {
    lthread <- do.call("lth", as.list(param[i,]))
    thList[[i]] <- lthread$th
    ll[,i] <- lthread$ll
    nNaN[,i] <- lthread$nNaN
}
color <- as.numeric(factor(param[,1]))

estcomp <- c(1, 2, 3, 5, 6, 7, 8, 9, 10, 14) # components estimated via IF2
compnames <- c('R0', 'amp', 'alpha', 'mu', 'gen', 'infec period\nprop', 'sigma2', 'rep', 'repOD', 'G', 'source_pow', 'dest_pow', 'dist_pow', 'cohort entry\nfrac', 'school start day', 'entryage', 'start year') # component names 
logistic <- function(x) { 1/(1+exp(-x)) }; identity <- function(x) {x}; logisticAmp <- function(x) { .68*1/(1+exp(-x)) }
ftrans <- list(exp, logisticAmp, identity, identity, exp, logistic, exp, logistic, exp, exp, identity, identity, identity, logistic, identity, identity, identity, logistic)


mcap <- function(lp,parameter,confidence=0.95,lambda=0.75,Ngrid=1000){
    smooth_fit <- loess(lp ~ parameter,span=lambda)
    parameter_grid <- seq(min(parameter), max(parameter), length.out = Ngrid)
    smoothed_loglik <- predict(smooth_fit,newdata=parameter_grid)
    smooth_arg_max <- parameter_grid[which.max(smoothed_loglik)]
    dist <- abs(parameter-smooth_arg_max)
    included <- dist <= sort(dist)[trunc(lambda*length(dist))]
    maxdist <- max(dist[included])
    weight <- rep(0,length(parameter))
    weight[included] <- (1-(dist[included]/maxdist)^3)^3
    quadratic_fit <- lm(lp ~ a + b, weight=weight,data = data.frame(lp=lp,b=parameter,a=-parameter^2))
    b <- unname(coef(quadratic_fit)["b"] )
    a <- unname(coef(quadratic_fit)["a"] )
    m <- vcov(quadratic_fit)
    var_b <- m["b","b"]
    var_a <- m["a","a"]
    cov_ab <- m["a","b"]
    se_mc_squared <- (1 / (4 * a^2)) * (var_b - (2 * b/a) * cov_ab + (b^2 / a^2) * var_a)
    se_stat_squared <- 1/(2*a)
    se_total_squared <- se_mc_squared + se_stat_squared
    delta <- qchisq(confidence,df=1) * ( a * se_mc_squared + 0.5)
    loglik_diff <- max(smoothed_loglik) - smoothed_loglik
    ci <- range(parameter_grid[loglik_diff < delta])
    list(lp=lp,parameter=parameter,confidence=confidence,
         quadratic_fit=quadratic_fit, quadratic_max=b/(2*a),
         smooth_fit=smooth_fit,
         fit=data.frame(
             parameter=parameter_grid,
             smoothed=smoothed_loglik,
             quadratic=predict(quadratic_fit, list(b = parameter_grid, a = -parameter_grid^2))
         ),
         mle=smooth_arg_max, ci=ci, delta=delta,
         se_stat=sqrt(se_stat_squared), se_mc=sqrt(se_mc_squared), se=sqrt(se_total_squared)
    )
}


plot.profile <- function(lp,parameter,confidence=0.95,np=500,ylab,lambda, ...){
  mcap1 <- mcap(lp,parameter,confidence,lambda)
  fit <- mcap1$fit

  ##if(missing(ylab))  ylab <- expression("l"^"P")

  plot(lp ~ parameter, ylab='', xlab='', xaxt='n', yaxt='n')
  ##axis(side=1, at=sqrt(gvec), labels=gvec)
  ##mtext(side=3, text=substitute(paste(lambda,'=',lmb,',  ', 'square-root scale'), list(lmb=lambda)))

  ##lines(fit$parameter, fit$quadratic, col = "blue", lwd = 1.5)
  lines(fit$parameter, fit$smoothed, col = "red", lwd = 1.5)

  diff_max <- max(fit$smoothed,na.rm=T) - fit$smoothed
  lower <- min(fit$parameter[diff_max < mcap1$delta],na.rm=T)
  upper <- max(fit$parameter[diff_max < mcap1$delta],na.rm=T)
  abline(v=c(lower,upper),col="red")
  abline(h=max(fit$smoothed,na.rm=T)-mcap1$delta,col="red")
}


## real data MCAP for G 
llp <- apply(ll, 2, function(lvec) lvec[max(which(!is.na(lvec)))]) # log likelihood estimate from the last completed iteration
llmat <- matrix(llp, 6, 6)
lltop3 <- c(apply(llmat, 2, function(lvec) sort(lvec, decreasing=TRUE)[1:3]))
gvec <- c(20, 50,100*c(1,2,3,4)); reps <- 3
gs <- rep(gvec, each=reps)
mcapresult <- mcap(lltop3, sqrt(gs), confidence=0.95, lambda=1.0)
    
#pdf('plots/mcap_G.pdf', width=8, height=5)
cex=1.35
par(mar=c(3+(cex-1)*2,3+(cex-1)*2+.1,0.5+(cex-1),0.5))
plot.profile(lltop3,sqrt(gs),lambda=1.0)
axis(side=1, at=sqrt(gvec), labels=NA)
axis(side=1, at=sqrt(gvec), labels=gvec, tick=FALSE, line=-0.25+(cex-1), cex.axis=cex)
axis(side=2, labels=NA)
axis(side=2, tick=FALSE, line=-0.25+(cex-1), cex.axis=cex)
mtext(side=1, line=2+(cex-1)*2, text='G', cex=cex)
mtext(side=2, line=2+(cex-1)*2, text="log likelihood estimate", cex=cex)
#dev.off()


## plot log likelihood estimates over IF2 iterations
rll = range(ll[is.finite(ll)])
plot(1:10, rep(NA,10), ylim=c(-90000,rll[2]))
for(i in 1:24) {
    lines(ll[,i], col=color[i])
}
legend('bottomright', legend=c(50,100,300,500), lty=1, col=1:4)
