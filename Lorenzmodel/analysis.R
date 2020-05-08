ntimes <- 200
ncity <- 50

lth <- function(F, n) { 
    filenamecore <- paste0('_K', ncity, '_F', format(round(F,6),nsmall=6), ifelse(profile,'_profile_','_slice_'), 'R', R, 'J', J, 'S', S, 'T', ntimes, 'if_', ifelse(runif,'T','F'), '_wi_', ifelse(weighted,'T','F'), '_la', max_la, '_', n, '_dtObs', format(dtObs,nsmall=1)) 
    ## compute likelihood
    lfilename <- paste('im_est_likelihood', filenamecore, '.txt', sep='')
    if(file.exists(lfilename)) { lest <- as.matrix(read.table(lfilename))
    } else { print(paste("Warning: no lest file exits.", "F =",F, 'n =', n)); lest <- NA }
    ##ll <- apply(cbind(matrix(c(lest), nrow=ntimes)[1:ntimes,]), 2, mean, na.rm=FALSE) * ntimes # ll for each m
    ll <- apply(cbind(matrix(c(lest), nrow=ntimes)[1:ntimes,]), 2, sum) # ll for each m
    nNaN <- apply(matrix(c(lest), nrow=ntimes), 2, function(x) sum(is.na(x)))
    if ((lenll <- length(ll)) < M) { ll[(lenll+1):M] <- NA }
    if ((lennn <- length(nNaN)) < M) { nNaN[(lennn+1):M] <- NA }
    if (runif) {
        thfilename <- paste('im_theta', filenamecore, '.txt', sep='')
        if(file.exists(thfilename)) { th <- read.table(thfilename)
        } else { print(paste('Warning: no theta file exists.', "F =", F, 'n =', n)); th <- matrix(NA, ntimes*M, thsize) }
        if ((lenth <- dim(th)[1]) < M*ntimes) { th[(lenth+1):(M*ntimes),] <- NA }
        return(list(ll=ll, nNaN=nNaN, th=th))
    }
    return(list(ll=ll, nNaN=nNaN))
}
lplot <- function(...) plot(..., type='l')

F <- 8; sig_p <- 1; sig_m <- 1
profile <- FALSE; runif <- FALSE; weighted <- TRUE; M <- ifelse(runif, 20, 1)
dtObs <- 0.5 # observation interval

## read GIRF results
options(stringsAsFactors=FALSE)
S <- ncity; max_la <- 2
startrow=0
R <- 5; J <- 2000
results <- data.frame(R=rep(5,5),J=rep(2000,5),n=0:4,ll=NA,method='GIRF')
results$ll <- sapply(0:4, function(n) lth(8,n)$ll/(ntimes*ncity))
R <- 1; J <- 2000
results[startrow+6:10,] <- data.frame(R=rep(1,5),J=rep(2000,5),n=0:4,ll=NA,method='GIRF')
results$ll[startrow+6:10] <- sapply(0:4, function(n) lth(8,n)$ll/(ntimes*ncity))
R <- 1; J <- 400
results[startrow+11:15,] <- data.frame(R=rep(1,5),J=rep(400,5),n=0:4,ll=NA,method='GIRF')
results$ll[startrow+11:15] <- sapply(0:4, function(n) lth(8,n)$ll/(ntimes*ncity))
results$RJid[startrow+1:15] <- rep(3:1,each=5)

## read EnKF results
startrow=15
results[startrow+1:15,] <- data.frame(R=1,J=rep(c(400,2000,10000),each=5),n=0:4,ll=NA, method='EnKF',RJid=rep(1:3,each=5))
rowno=startrow
for (J_enkf in c(400,2000,10000)) {
    for (expno_enkf in 0:4) {
        rowno <- rowno+1
        load(paste0('EnKF_K', ncity, '_F', format(round(F,6),nsmall=6), '_sig_p', format(round(sig_p,6),nsmall=6), 'J', J_enkf, 'T', ntimes, '_', expno_enkf, '_dtObs', dtObs, '.RData'))
        results$ll[rowno] <- ll/(ntimes*ncity)
    }
}

## read local EnKF results (only for d=50)
##if (ncity==50) {
##startrow=30
##results[startrow+1:15,] <- data.frame(R=1,J=rep(c(400,2000,10000),each=5),n=0:4,ll=NA, method='LEnKF',RJid=rep(1:3,each=5))
##rowno=startrow
##for (J_enkf in c(400,2000,10000)) {
##    for (expno_enkf in 0:4) {
##        rowno <- rowno+1
##        B = 3; rho = 1.1
##        load(paste0('EnKF_local_varinfl_K', ncity, '_F', format(round(F,6),nsmall=6), '_sig_p', format(round(sig_p,6),nsmall=6), 'J', J_enkf, 'T', ntimes, '_', expno_enkf, '_dtObs', dtObs, '_B', B, '_varinfl', format(round(rho,2),nsmall=2), '.RData'))
##        results$ll[rowno] <- ll/(ntimes*ncity)
##    }
##}
##}

## read boostrap PF results (only used for d=4)
if (ncity==4) {
S <- 1; max_la <- 1
startrow=30
R <- 5; J <- 2000
results[startrow+1:5,] <- data.frame(R=rep(5,5),J=rep(2000,5),n=0:4,ll=NA,method='BPF')
results$ll[startrow+1:5] <- sapply(0:4, function(n) lth(8,n)$ll/(ntimes*ncity))
R <- 1; J <- 2000
results[startrow+6:10,] <- data.frame(R=rep(1,5),J=rep(2000,5),n=0:4,ll=NA,method='BPF')
results$ll[startrow+6:10] <- sapply(0:4, function(n) lth(8,n)$ll/(ntimes*ncity))
R <- 1; J <- 400
results[startrow+11:15,] <- data.frame(R=rep(1,5),J=rep(400,5),n=0:4,ll=NA,method='BPF')
results$ll[startrow+11:15] <- sapply(0:4, function(n) lth(8,n)$ll/(ntimes*ncity))
results$RJid[startrow+1:15] <- rep(3:1,each=5)
}

## plot likelihood estimates
if (ncity==4) { results$method <- factor(results$method, levels=c('GIRF','EnKF','BPF')) }
if (ncity==50) { results$method <- factor(results$method, levels=c('GIRF','EnKF'))} #,'LEnKF')) }
library(ggplot2)
ggplot(results, aes(RJid, ll)) +
    geom_point(aes(color=method, shape=method)) +
    scale_y_continuous(name='log likelihood estimate /Nd') + 
    scale_x_continuous(name='total number of particles',breaks=1:3,labels=c(400,2000,10000)) +
    ggtitle(bquote("d="*.(ncity)*","*~~~~Delta[obs]*"="*.(dtObs))) + theme(plot.title = element_text(hjust = 0.5))

give_pdf=TRUE
if(give_pdf) { ggsave(paste0('figures/GIRF_EnKF_K',ncity,'_F',format(round(F,1),nsmall=1),'_sig_p',format(round(sig_p,2),nsmall=2),'T',ntimes,'_dtObs',dtObs,'.pdf'), width=4, height=3) }



## read estimated state means
state <- read.table(paste('../OBSINT_0.5/state_dim',ncity,'F',format(F,nsmall=6),'sig_p',format(sig_p,nsmall=6),'sig_m',format(sig_m,nsmall=6),'.txt', sep=''))
obs <- read.table(paste('../OBSINT_0.5/obs_dim',ncity,'F',format(F,nsmall=6),'sig_p',format(sig_p,nsmall=6),'sig_m',format(sig_m,nsmall=6),'.txt', sep=''))
filenamecore <- paste0('_K', ncity, '_F', format(round(F,6),nsmall=6), ifelse(profile,'_profile_','_slice_'), 'R', R, 'J', J, 'S', S, 'T', ntimes, 'if_', ifelse(runif,'T','F'), '_wi_', ifelse(weighted,'T','F'), '_la', max_la, '_', n)#, '_dtObs', format(dtObs,nsmall=1)) 
ss <- read.table(paste('im_est_state_mean', filenamecore, '.txt', sep=''))
load(paste0('EnKF_K', ncity, '_F', format(round(F,6),nsmall=6), '_sig_p', format(round(sig_p,6),nsmall=6), 'J', J_enkf, 'T', ntimes, '_', expno_enkf, '_dtObs', format(dtObs_enkf,nsmall=1), '.RData'))
library(pomp)
fmean_enkf <- t(filter.mean(EnKF.result_C))
##pdf(paste('plots/GIRF_EnKF_dim',ncity,'_F8_OBSINT0.5_la',max_la,'GIRF_J_2000_EnKF_J_16000_0.pdf',sep=''), width=12, height=5)
par(mar=c(4,4,.5,.5))
city <- 1
rn <- 1:100 # plotting range
cs <- 1.5
lplot(0.5*rn, state[1+rn, city], col='red', lty=1, xaxt='n', yaxt='n', xlab='',ylab='', ylim=c(-10,13))
points(0.5*rn, obs[rn, city], col='blue')
lines(0.5*rn, ss[rn,city], col='purple', lty=1)
lines(0.5*rn, fmean_enkf[rn,city], lty=2, col='black')
axis(side=1, cex.axis=cs)
axis(side=2, cex.axis=cs)
mtext(text='State estimates', side=2, line=2.5, cex=cs)
mtext(text='time', side=1, line=2.5, cex=cs)
##mtext(text=paste(ncity,'dimensional stochastic Lorenz model'), side=3)
legend('bottomright', legend=c('observation','true state', 'GIRF estimation', 'EnKF estimation'), pch=c(1,NA,NA,NA), lty=c(NA,1,1,2), col=c('blue','red', 'purple','black'), bg='white', cex=cs, horiz=TRUE, x.intersp=1, box.lty=0)
##legend('bottomright', title='Log likelihood estimate', legend=c(paste('GIRF: ', round(lth(F,n)$ll)), paste('EnKF:', round(result$ll))), lty=1, col=c('white','white'), box.lty=0)
##dev.off()

ESS <- unlist(read.table(paste('im_ESS', filenamecore, '.txt', sep='')))




## plot theta ##
F <- 8; sig_p <- 1; sig_m <- 1
R <- 5; J <- 2000; S <- ncity; runif <- TRUE; weighted <- TRUE; max_la <- 2; n <- 1; M <- ifelse(runif, 20, 1)
filenamecore <- paste('_K', ncity, '_F', format(round(F,6),nsmall=6), 'R', R, 'J', J, 'S', S, 'T', ntimes, 'if_', ifelse(runif,'T','F'), '_wi_', ifelse(weighted,'T','F'), '_la', max_la, '_', n, sep='') 
thfilename <- paste('im_theta', filenamecore, '.txt', sep='')
theta <- read.table(thfilename, sep=' ')
trfun <- c(exp, exp)
parname <- c("sig_p", "sig_m")
param <- sapply(1:2, function(j) { sapply(theta[,j+1], trfun[[j]]) } )
par(mfrow=c(2,1), mar=c(0,2,0,1), oma=c(3,1,1,0))
plot(param[,1], type='l', xlab='', ylab='', xaxt='n')
axis(side=1, at=ntimes*1:M, labels=FALSE)
mtext(side=2,line=2, text=parname[1])
plot(param[,2], type='l', xlab='', ylab='', xaxt='n')
mtext(side=2, line=2, text=parname[2])
mtext(side=1, line=2, text='IF2 iterations')
axis(side=1, at=ntimes*1:M, labels=1:M)



## plot likelihood estimates ##
Fvalues <- seq(6,10,by=.5)
param4 <- cbind(rep(Fvalues, each=4), rep(0:3,length(Fvalues))) # four repetitions were made for each value of F
Frep2 <- rep(Fvalues, each=2)
profile <- FALSE; lls4 <- matrix(apply(param4, 1, function(x) { ll <- do.call(lth, as.list(x))$ll; ifelse(is.finite(ll), ll, NA)}), nrow=4) # slice
llstop2 <- apply(lls4, 2, function(x) which(rank(x, na.last=FALSE) %in% 3:4)) # choose top two points per a value of F
lls <- c(sapply(1:length(Fvalues), function(n) lls4[llstop2[,n],n]))
profile <- TRUE; llp4 <- matrix(apply(param4, 1, function(x) do.call(lth, as.list(x))$ll), nrow=4) #profile
llptop2 <- apply(llp4, 2, function(x) which(rank(x, na.last=FALSE) %in% 3:4)) # choose top two points per a value of F
llp <- c(sapply(1:length(Fvalues), function(n) llp4[llptop2[,n],n]))

pdf(paste('plots/', 'slice', '_dim', ncity, '.pdf', sep=''), width=6, height=5)

cex=1.35
par(mar=c(3+(cex-1)*2,3+(cex-1)*2+.1,1.5,0.5))
## slice likelihood
plot(param4[,1], c(lls4), pch=1, xlab='', ylab='', xaxt='n', yaxt='n')
abline(v=F, lty=2)
axis(side=1, at=Frep2, labels=NA)
axis(side=1, at=Frep2, labels=Frep2, tick=FALSE, line=-0.25+(cex-1), cex.axis=cex)
axis(side=2, labels=NA)
axis(side=2, tick=FALSE, line=-0.25+(cex-1), cex.axis=cex)
mtext(side=1, line=1.8+(cex-1)*2, text='F', cex=cex)
mtext(side=2, line=1.8+(cex-1)*2, text="log likelihood estimate", cex=cex)

dev.off()

pdf(paste('plots/', 'mcap_profile', '_dim', ncity, '.pdf', sep=''), width=6, height=4)

cex=1.35
par(mar=c(3+(cex-1)*2,3+(cex-1)*2+.1,0.5,0.5))
## profile likelihood
plot.profile(llp,Frep2,lambda=0.75, plotci=TRUE, pch=4) #ylim=range(c(llp,lls)), 
abline(v=F, lty=2)
axis(side=1, at=Frep2, labels=NA)
axis(side=1, at=Frep2, labels=Frep2, tick=FALSE, line=-0.25+(cex-1), cex.axis=cex)
axis(side=2, labels=NA)
axis(side=2, tick=FALSE, line=-0.25+(cex-1), cex.axis=cex)
mtext(side=1, line=1.8+(cex-1)*2, text='F', cex=cex)
mtext(side=2, line=1.8+(cex-1)*2, text="log likelihood estimate", cex=cex)

dev.off()

## slice likelihood
pdf(paste('plots/', 'mcap_slice', '_dim', ncity, '.pdf', sep=''), width=6, height=4)

sFpts <- 1:14 # for slice likelihood curve, we omit F=10 (for dim=50), F=9.5,10 (for dim=100)
par(mar=c(3+(cex-1)*2,3+(cex-1)*2+.1,0.5,0.5))
plot.profile(lp=lls[sFpts], parameter=Frep2[sFpts], paramgrid=seq(min(Fvalues), max(Fvalues), length.out=1000), lambda=0.75, color='blue', add=0, pch=1, plotci=TRUE)
abline(v=F, lty=2)
axis(side=1, at=Frep2, labels=NA)
axis(side=1, at=Frep2, labels=Frep2, tick=FALSE, line=-0.25+(cex-1), cex.axis=cex)
axis(side=2, labels=NA)
axis(side=2, tick=FALSE, line=-0.25+(cex-1), cex.axis=cex)
mtext(side=1, line=1.8+(cex-1)*2, text='F', cex=cex)
mtext(side=2, line=1.8+(cex-1)*2, text="log likelihood estimate", cex=cex)
##legend('topright', legend=c('profile', 'slice'), lty=1, col=c('red', 'blue'), pch=c(4,1), bg='white', cex=cex)

dev.off()

plot.profile <- function(lp,parameter, paramgrid=NULL, confidence=0.95,np=500, lambda, color='red', add=FALSE, pch=1, ylim=NULL, plotci=TRUE){
    mcap.obj <- mcap(lp,parameter,confidence,lambda)
    fit <- mcap.obj$fit

    if (!add) { plot(lp ~ parameter, ylab='', xlab='', xaxt='n', yaxt='n', pch=pch, ylim=ylim, col=color) }
    if (add) { points(lp ~ parameter, pch=pch,ylim=ylim, col=color) }

    if (is.null(paramgrid)) { lines(fit$parameter, fit$smoothed, col = color, lwd = 1.5) }
    if (!is.null(paramgrid)) { lines(paramgrid, predict(mcap.obj$smooth_fit, newdata=paramgrid), col = color, lwd = 1.5) }

    diff_max <- max(fit$smoothed,na.rm=T) - fit$smoothed
    lower <- min(fit$parameter[diff_max < mcap.obj$delta],na.rm=T)
    upper <- max(fit$parameter[diff_max < mcap.obj$delta],na.rm=T)
    if (plotci) {
        abline(v=c(lower,upper),col=color)
        abline(h=max(fit$smoothed,na.rm=T)-mcap.obj$delta,col=color)
    }
}

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




## plot MC-MLEs
runif <- TRUE; MCMLE <- t(sapply(1:length(Frep2), function(x) { th=lth(Frep2[x], c(llptop2)[x]-1)$th; unlist(th[dim(th)[1],]) }))
plot.smooth <- function(x, y, lambda=0.75, add=FALSE, color='black', pch=1, ylim=NULL, cex=1) {
    smooth_fit <- loess(y ~ x, span=lambda)
    x_grid <- seq(min(x), max(x), length.out = 1000)
    smoothed_y <- predict(smooth_fit,newdata=x_grid)
    if (!add) { plot(x, y, col=color, pch=pch, ylim=ylim, cex=cex, xaxt='n', yaxt='n', xlab='', ylab='') }
    if (add) { points(x, y, col=color, pch=pch, cex=cex) }
    lines(x_grid, smoothed_y, col = color, lwd = 1.5)
}

library('latex2exp')
    
pdf(paste('plots/', 'MCMLE', '_dim', ncity, '.pdf', sep=''), width=6, height=4)

cex=1.35
par(mar=c(3+(cex-1)*2,3+(cex-1)*2+.1,0.5+(cex-1),0.5))
plot.smooth(MCMLE[,1], exp(MCMLE[,2]), color='blue', cex=1.5, ylim=exp(range(MCMLE[,2:3]))) #ylim=c(min(exp(MCMLE[,2:3])), 2)) ## sig_p
plot.smooth(MCMLE[,1], exp(MCMLE[,3]), color='red', add=TRUE, pch=4, cex=1.5) ## sig_m
abline(h=1.0, lty=2); abline(v=8, lty=2)
axis(side=1, at=Frep2, labels=NA)
axis(side=1, at=Frep2, labels=Frep2, tick=FALSE, line=-0.25+(cex-1), cex.axis=cex)
axis(side=2, labels=NA)
axis(side=2, tick=FALSE, line=-0.25+(cex-1), cex.axis=cex)
mtext(side=1, line=1.8+(cex-1)*2, text='F', cex=cex)
mtext(side=2, line=1.8+(cex-1)*2, text="Monte Carlo MLE", cex=cex)
legend('topleft', legend=c(TeX('$\\sigma_p$'), TeX('$\\sigma_s$')), lty=1, pch=c(1,4), col=c('blue','red'), cex=cex, bg='white')

dev.off()
