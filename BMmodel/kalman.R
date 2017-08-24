## State inference
LMSE <- function(nvar, alpha, diagCov=FALSE) {
    ## alpha: non-diagonal entries of the jump covariance matrix
    ## diagCov: was diagonal covariance matrix used for computing predictive likelihood?
    ntimes = 50
    state0 = rep(0,nvar) #2 * (-(nvar-1)/2) : ((nvar-1)/2);

    s = 1.0 # std dev of the marginal distribution of jump in one dimension
    d = 1.0 # measurement error standard deviation

    Sig = matrix(s^2*alphaPt, nvar, nvar); diag(Sig) <- s^2
    D = d^2 * diag(rep(1, nvar))

    ofilename = paste('data/obs_K', nvar, 'alpha', format(alpha, nsmall=2), '.txt', sep='')
    Y = c(t(as.matrix(read.table(ofilename))))[1:(ntimes*nvar)]

    lll <- c()
    library(mvtnorm)
    Ymat <- t(matrix(Y, nvar, ntimes))
    ## initialization
    P <- Sig # covariance of X_t conditional on the data available so far
    mu <- state0 # mean of X_t conditional on the data available so far
    ll <- 0 # log likelihood of y
    EXFILT <- c() # filter mean (E(X_t|y_{1:t}))
    SXXFILT <- list() # posterior covariance (Var(X_t|y_{1:t}))
    for (t in 1:ntimes) {
        residual <- Ymat[t,] - mu
        S <- P + D # residual covariance
        K <- t(solve(S,P))
        ## compute likelihood (y_t | y_{1:(t-1)})
        ll <- ll + dmvnorm(Ymat[t,], mean=mu, sigma=S, log=TRUE)
        lll <- c(lll,dmvnorm(Ymat[t,], mean=mu, sigma=S, log=TRUE))
        ## update
        mu <- mu + K%*%residual
        P <- P - K%*%P
        EXFILT <- rbind(EXFILT, c(mu))
        SXXFILT[[t]] <- P
        ## predict
        P <- P + Sig
    }

    ll
    
    date = '170216'
    directory = date

    le <- numeric(); sum_X50_mean <- numeric(nvar); sum_X50_mean_sq <- numeric(nvar)
    MSEE <- numeric()
    rep <- 1
    for (n in 1:rep) {
        R = 60; J = 1000; N = nvar; options(scipen=7)
        commonString = paste("_K", nvar, "alpha", format(alpha, nsmall=2), "R", R, "J", J, "N", N, "T", ntimes, "if_", "F", "wi_", "T", "_la", 2, "_", (n-1), ifelse((!diagCov) || (alpha==0), ".txt", "_diag.txt"), sep='')
        ##commonString = paste("_K", nvar, "alpha", format(alpha, nsmall=2), "R", R, "J", J, "N", N, "_1_unsmoothed.txt", sep='')
        lestimate <- read.table(paste("data/",directory,"/im_est_likelihood", commonString, sep=''))
        lestimate <- lestimate[dim(lestimate)[1]-(ntimes-1):0,]
        le[n] <- sum(lestimate)
        paste('exact log likelihood =', ll, 'likelihood estimate =', le)
        ##
        state_mean= as.matrix(read.table(paste("data/",directory,"/im_est_state_mean", commonString, sep=''))); colnames(state_mean) <- NULL
        sum_X50_mean <- sum_X50_mean + state_mean[50,]
        sum_X50_mean_sq <- sum_X50_mean_sq+ state_mean[50,]^2
        ##f1_mean= as.matrix(read.table(paste("data/",directory,"/im_est_f1_mean", commonString, sep=''))); colnames(f1_mean) <- NULL
        ##state_var=(R*J)/(R*J-1)*(f1_mean - state_mean^2)
        MSEE[n] <- sum((EXFILT[50,]-state_mean[50,])^2)/nvar
    }

    MSE <- c(k=nvar, alpha=alpha, bias_sq = (bias_sq <- sum((EXFILT[50,]-sum_X50_mean/rep)^2)/nvar), var = (variance <- sum((sum_X50_mean_sq - sum_X50_mean^2/rep)/rep)/nvar), MSE = bias_sq+variance)
    L <- c(k=nvar, alpha=alpha, exact_l = ll, log_mean_est_l = mean(le)+log(mean(exp(le-mean(le)))), se_log_mean_est_l = sqrt(var(exp(le-mean(le)))/rep)/mean(exp(le-mean(le))), mean_log_est_l = mean(le), se_log_est_l = sqrt(var(le)/rep) )

    return(list(MSE = MSE, L = L))
}

MSE_dim <- sapply(c(20,50,100,200), function(d) LMSE(d, alpha=0.00, diagCov=FALSE)$MSE) ## saved as MSE_dim.txt 
L_dim <- sapply(c(20,50,100,200), function(d) LMSE(d, alpha=0.00, diagCov=FALSE)$L) 
MSE_alpha_exact20 <- sapply(0.1*0:5, function(a) LMSE(20, a, diagCov=FALSE)$MSE)
L_alpha_exact20 <- sapply(0.1*0:5, function(a) LMSE(20, a, diagCov=FALSE)$L) 
MSE_alpha_diag20 <- sapply(0.1*0:5, function(a) LMSE(20, a, diagCov=TRUE)$MSE)
L_alpha_diag20 <- sapply(0.1*0:5, function(a) LMSE(20, a, diagCov=TRUE)$L) 
MSE_alpha_exact50 <- sapply(0.1*0:5, function(a) LMSE(50, a, diagCov=FALSE)$MSE)
L_alpha_exact50 <- sapply(0.1*0:5, function(a) LMSE(50, a, diagCov=FALSE)$L) 
MSE_alpha_diag50 <- sapply(0.1*0:5, function(a) LMSE(50, a, diagCov=TRUE)$MSE)
L_alpha_diag50 <- sapply(0.1*0:5, function(a) LMSE(50, a, diagCov=TRUE)$L) 

MSE_dim <- read.table('data/161210/MSE_dim.txt', row.names=1)
L_dim <- read.table('data/161210/L_dim.txt', row.names=1)
MSE_alpha_exact20 <- read.table('data/161210/MSE_alpha_exact20.txt', row.names=1)
L_alpha_exact20 <- read.table('data/161210/L_alpha_exact20.txt', row.names=1)
MSE_alpha_diag20 <- read.table('data/161210/MSE_alpha_diag20.txt', row.names=1)
L_alpha_diag20 <- read.table('data/161210/L_alpha_diag20.txt', row.names=1)
MSE_alpha_exact50 <- read.table('data/161210/MSE_alpha_exact50.txt', row.names=1)
L_alpha_exact50 <- read.table('data/161210/L_alpha_exact50.txt', row.names=1)
MSE_alpha_diag50 <- read.table('data/161210/MSE_alpha_diag50.txt', row.names=1)
L_alpha_diag50 <- read.table('data/161210/L_alpha_diag50.txt', row.names=1)

##pdf('plots/161210/MSE_dim(alpha0.00R60J1000).pdf')
par(mar=c(3.5,3.9,1,1))
plot(MSE_dim["k",], MSE_dim["MSE",], type='b', ann=FALSE, xlim=c(20,200), ylim=c(0,0.01), xaxt='n')
points(MSE_dim["k",], MSE_dim["bias_sq",], type='b', lty=2, pch=2)
abline(h=0)
legend('topleft', legend=c('MSE', expression(bias^2)), lty=c(1,2), pch=c(1,2), cex=1.3)
axis(side=1, at=c(0,20,50,100,200))
mtext('Dimension', side=1, line=2.5, cex=1.3)
mtext('Mean squared error (average over all sites)', side=2, line=2.5, cex=1.3)
##dev.off()

##pdf('plots/161210/MSE_alpha(dim20R60J1000).pdf')
par(mar=c(3.5,3.9,1,1))
plot(MSE_alpha_diag20["alpha",], MSE_alpha_diag20["MSE",], type='b', lty=5, pch=4, ann=FALSE, xlim=c(0,0.5), ylim=c(0,0.0088))
points(MSE_alpha_exact20["alpha",], MSE_alpha_exact20["MSE",], type='b')
abline(h=0)
legend('topleft', legend=c('Diagonal covariance used', 'Exact covariance used'), lty=c(5,1), pch=c(4,1), cex=1.3)
mtext('Correlation coefficient', side=1, line=2.5, cex=1.3)
mtext('Mean squared error (average over all sites)', side=2, line=2.5, cex=1.3)
##dev.off()

##pdf('plots/161210/MSE_alpha(dim50R60J1000).pdf')
par(mar=c(3.5,3.9,1,1))
plot(MSE_alpha_diag50["alpha",], MSE_alpha_diag50["MSE",], type='b', lty=5, pch=4, ann=FALSE, xlim=c(0,0.5), ylim=c(0,0.046))
points(MSE_alpha_exact50["alpha",], MSE_alpha_exact50["MSE",], type='b')
abline(h=0)
legend('topleft', legend=c('Diagonal covariance used', 'Exact covariance used'), lty=c(5,1), pch=c(4,1), cex=1.3)
mtext('Correlation coefficient', side=1, line=2.5, cex=1.3)
mtext('Mean squared error (average over all sites)', side=2, line=2.5, cex=1.3)
##dev.off()


t=ntimes
##pdf(paste('plots/',directory,'/posterior_mean',substr(commonString,start=1,stop=nchar(commonString)-4),'.pdf',sep=''))
ui=state_mean[t,] + sqrt(state_var[t,])
li=state_mean[t,] - sqrt(state_var[t,])
plotCI(1:nvar + 0.2, state_mean[t,], ui=ui, li=li, ylim=range(c(ui,li)), slty=2, pch=4, xlab='', ylab='', main='')
if(direct) {SXXFILT_t <- SXXPOST[(t-1)*nvar+1:nvar,(t-1)*nvar+1:nvar]}
if(!direct) {SXXFILT_t <- SXXFILT[[t]]}
plotCI(1:nvar, EXFILT[t,], ui=EXFILT[t,] + sqrt(diag(SXXFILT_t)), li=EXFILT[t,] - sqrt(diag(SXXFILT_t)), add=TRUE)
legend('topright', legend=c('Exact', 'Estimated'), lty=c(1,2), pch=c(1,4), bg='white')
title(main=paste('t = ',t, ',  J = ', J, ', R = ', R,  sep=''), xlab='Sites', ylab='Posterior mean')
##points(1:nvar, latent_states[t,], col=2)
##points(1:nvar, Ymat[t,], col=3, pch=4)

dev.off()



## bias variance plot
## plotting ##
##n <- 5; R <- 60; J <- 1000; N <- nvar; date <- 161210
##ESSrawdata <- read.table(paste('data/', date, '/im_ESS_K', nvar, 'alpha', format(alpha, nsmall=2), "R", R, "J", J, "N", N, "T", ntimes, "if_", "F", "wi_", "T", "_la", 2, "_", (n-1), ".txt", sep='')) #The row of the ESS data file writes the ESS of the intermediate time steps for each island and each observation time. (the ESS data file is (ntimes*N)-by-R matrix. Note: Before Oct. 18, 2016, the ESS data file is (ntimes*R)-by-N matrix)
##ESS <- apply(ESSrawdata, 1, sum)
##ESS <- matrix(0, ntimes, N) # old version
##for(t in 1:ntimes) 
##    ESS[t,] <- apply(ESSrawdata[(t-1)*R+1:R,], 2, sum)

require(plotrix)
sfilename = paste('data/state_K', nvar, 'alpha', format(alpha, nsmall=2), '.txt', sep='')
latent_states = as.matrix(read.table(sfilename))



## IF2 run
directory <- 170311
alpha = .1; nvar = 20; ntimes = 50;
R = 60; J = 1000; N = nvar; options(scipen=7)
M <- 30
nvec <- 0
le <- matrix(0,M,length(nvec)); thList <- list()
for(n in nvec) {
    commonString = paste("_K", nvar, "alpha", format(alpha, nsmall=2), "R", R, "J", J, "N", N, "T", ntimes, "if_", "T", "wi_", "T", "_la", 2, "_", n, sep='')
    lestimate <- read.table(paste("data/",directory,"/im_est_likelihood", commonString, '.txt', sep=''))
    le[,n+1] <- apply(matrix(unlist(lestimate),nrow=ntimes),2,sum)
    th <- read.table(paste("data/",directory,"/im_theta", commonString, '.txt', sep=''))
    thList[[n+1]] <- th
}

logistic <- function(x) { 1/(1+exp(-x)) }

ftrans <- list(exp,exp,logistic)
ran <- sapply(1:3, function(comp) { thcomp <- c(); for(n in nvec) {thcomp <- c(thcomp,ftrans[[comp]](thList[[n+1]][,comp])) }; range(thcomp) })

plotcomp <- function(comp) {
    plot(ftrans[[comp]](thList[[1]][,comp]),type='l',ylim=ran[,comp])
    for(n in nvec) {
        if (n==0) next
        points(ftrans[[comp]](thList[[n+1]][,comp]),type='l',col=(n+1))
    }
}

par(mfrow=c(2,2))
sapply(1:3, plotcomp)
plot(le[,1],xlim=c(0,M),ylim=range(c(le)))
for(n in nvec) { if (n==0) next;
    points(le[,n+1], type='l', col=(n+1))
}
    
        

## direct way of computing the filter means
direct <- FALSE # should the filter mean and the likelihood be computed directly by inverting the covariance matrix? If FALSE, they are computed by Kalman filter.

if (direct) {
    TemporalCov = matrix(0, ntimes, ntimes)
    TemporalCov <- pmin(row(TemporalCov), col(TemporalCov))
    M <- kronecker(TemporalCov, Sig)
    COV = rbind(M,M)
    COV = cbind(COV,COV)
    for (t in 1:ntimes)
        COV[nvar*ntimes + (t-1)*nvar + (1:nvar), nvar*ntimes + (t-1)*nvar + (1:nvar)] = COV[nvar*ntimes + (t-1)*nvar + (1:nvar), nvar*ntimes + (t-1)*nvar + (1:nvar)] + D
    SYY = COV[nvar*ntimes + (1:(nvar*ntimes)), nvar*ntimes + (1:(nvar*ntimes))]
    SXXPOST = M - M %*% solve(SYY, M) # COV(X|Y)
    ##log likelihood
    lldirect <- -ntimes*nvar/2*log(2*pi)-.5*log(det(SYY))-.5*Y%*%solve(SYY,Y)
    EXlast = state0
    EX = rep(state0, ntimes)
    EY = EX
    EY_IND = rep(state0, ntimes)
    EX_IND = rep(state0, ntimes)
    ## marginal distributions
    MXCov <- matrix(numeric(ntimes*ntimes), ntimes, ntimes) # Cov of X_k,1:T under no interaction
    for (t in 1:ntimes)
        for (u in 1:ntimes)
            MXCov[t,u] <- min(t, u) * s^2
    MYCov <- MXCov + diag(rep(d^2,ntimes)) # Cov of Y_k,1:T under no interaction
    ## posterior
    EXPOST =  EX + M %*% solve(SYY, (Y - EY)) # E(X|Y)
    EXFILT = c() # E(X_t|Y_1:t), theoretical mean of filtered particles
    for (t in 1:ntimes) 
        EXFILT <- c(EXFILT, EX[(t-1)*nvar+1:nvar] + M[(t-1)*nvar+1:nvar, 1:(t*nvar)] %*% solve(SYY[1:(t*nvar), 1:(t*nvar)], (Y- EY)[1:(t*nvar)]))
    EXFILT <- t(matrix(EXFILT, nvar, ntimes))
    EXFILT_IND = numeric(ntimes*nvar) # EXFILT where interaction between cities have been ignored
    for (k in 1:nvar)
        for (t in 1:ntimes)
            EXFILT_IND[k + (t-1)*nvar] <- EX_IND[k + (t-1)*nvar] + MXCov[t, 1:t] %*% solve(MYCov[1:t, 1:t], (Y- EY_IND)[k + (1:t - 1)*nvar])
    EXFILT_IND <- t(matrix(EXFILT_IND, nvar, ntimes))
    ##save(SXXPOST, EXPOST, EXFILT, EXFILT_IND, file=paste('data/filtermean_K', nvar, 'alpha', format(alpha, nsmall=2), '.RData', sep=''))
}



date = '160829'; directory = date

##condition <- c('noweighted_noshuffle', 'noweighted_shuffle', 'weighted_noshuffle', 'weighted_shuffle')
##condition <- c('multinomial', 'systematic')
condition <- c('smoothed', 'unsmoothed')

R <- 40; J <- 800; N <- nvar

state_mean_avgBias2 <- numeric(length(condition)) # bias^2 of mean estimates (average over all sites & times)
state_mean_avgVariance <- numeric(length(condition)) # variance of mean estimates (avearge over all sites & times)
state_mean_avgMSE <- numeric(length(condition)) # mse of mean estimates (average over all sites & times)
f1_mean_avgBias2 <- numeric(length(condition)) # bias^2 of mean estimates (average over all sites & times)
f1_mean_avgVariance <- numeric(length(condition)) # variance of mean estimates (avearge over all sites & times)
f1_mean_avgMSE <- numeric(length(condition)) # mse of mean estimates (average over all sites & times)

for (i in 1:length(condition)) {
    ## each particle filter with a (R,J) pair was repeated 20 times
    lestimate_sum <- 0
    state_mean_sum <- matrix(0, nrow=ntimes, ncol=nvar)
    f1_mean_sum <- matrix(0, nrow=ntimes, ncol=nvar)
    state_mean_sqsum <- matrix(0, nrow=ntimes, ncol=nvar)
    f1_mean_sqsum <- matrix(0, nrow=ntimes, ncol=nvar)
    for(j in 1:20) {
        comstr <- paste("_K", nvar, "alpha", format(alpha, nsmall=2), "R", R, "J", J, "N", N, "_", j, "_", condition[i], ".txt", sep='')
        lestimate <- as.matrix(read.table(paste("data/",directory,"/im_est_likelihood", comstr, sep='')))[1:ntimes,]
        state_mean=  as.matrix(read.table(paste("data/",directory,"/im_est_state_mean", comstr, sep='')))[1:ntimes,]; colnames(state_mean) <- NULL
        f1_mean= as.matrix(read.table(paste("data/",directory,"/im_est_f1_mean", comstr, sep='')))[1:ntimes,]; colnames(f1_mean) <- NULL
        lestimate_sum <- lestimate_sum + sum(apply(lestimate,1,mean))
        state_mean_sum <- state_mean_sum + state_mean
        f1_mean_sum <- f1_mean_sum + f1_mean
        state_mean_sqsum <- state_mean_sqsum + state_mean^2
        f1_mean_sqsum <- f1_mean_sqsum + f1_mean^2
    }
    state_mean_bias = state_mean_sum / 20 - EXFILT
    f1_mean_bias = f1_mean_sum / 20 - (EXFILT^2+t(matrix(diag(SXXPOST), nrow=nvar)))
    state_mean_variance = state_mean_sqsum / 20 - (state_mean_sum / 20)^2
    f1_mean_variance = f1_mean_sqsum / 20 - (f1_mean_sum / 20)^2
    state_mean_avgBias2[i] = mean(state_mean_bias[ntimes,]^2)
    state_mean_avgVariance[i] = mean(state_mean_variance[ntimes,])
    state_mean_avgMSE[i] = state_mean_avgBias2[i] + state_mean_avgVariance[i]
    f1_mean_avgBias2[i] = mean(f1_mean_bias[ntimes,]^2)
    f1_mean_avgVariance[i] = mean(f1_mean_variance[ntimes,])
    f1_mean_avgMSE[i] = f1_mean_avgBias2[i] + f1_mean_avgVariance[i]
}


##pdf(paste('plots/',directory,'/MSE_state_mean_K',nvar,'_', "alpha", format(alpha, nsmall=2), "R", R, "J", J, "N", N, '_', date,'.pdf', sep=''))
par(mar=c(4,4,4,4)+.1)
xaxis = 1:length(condition)
plot(xaxis, state_mean_avgMSE[xaxis], type='b', xaxt='n', xlab='', ylab='', ylim=c(0,max(state_mean_avgMSE[xaxis])))
points(xaxis, state_mean_avgBias2[xaxis], type='b', lty=2)
points(xaxis, state_mean_avgVariance[xaxis], type='b', lty=3)
abline(h=0)
mtext('MSE of Mean Estimates\n(average over all sites at time 10)', side=2, line=2)
title(paste(main='Filter performance (dimension = ', nvar, ')\n Average over 20 filtering repititions', sep=''))
xaxislabel <- c('not weighted\n no shuffle', 'not weighted\n shuffle', 'weighted\n no shuffle', 'weighted\n shuffle')
axis(side=1, at=xaxis, labels=xaxislabel)
legend('topright', legend=c('MSE', 'Bias^2', 'Variance'), lty=1:3)

dev.off()

##pdf(paste('plots/',directory,'/MSE_f1_mean_K',nvar,'_', "alpha", format(alpha, nsmall=2), "R", R, "J", J, "N", N, '_', date,'.pdf', sep=''))
par(mar=c(4,4,4,4)+.1)
xaxis = 1:length(condition)
plot(xaxis, f1_mean_avgMSE[xaxis], type='b', xaxt='n', xlab='', ylab='', ylim=c(0,max(f1_mean_avgMSE[xaxis])))
points(xaxis, f1_mean_avgBias2[xaxis], type='b', lty=2)
points(xaxis, f1_mean_avgVariance[xaxis], type='b', lty=3)
abline(h=0)
mtext('MSE of X^2 Estimates \n(average over all sites and times)', side=2, line=2)
title(paste(main='Filter performance (dimension = ',nvar,')\n Average over 20 filtering repititions', sep=''))
axis(side=1, at=xaxis, labels=xaxislabel)
legend('topright', legend=c('MSE', 'Bias^2', 'Variance'), lty=1:3)

dev.off()




## bias variance plot (Feb 2016)
# R {400, 200, 80, 40, 20, 10, 2};
# R {800, 400, 160, 80, 40, 20, 4};
# J {80, 160, 400, 800, 1600, 3200, 16000};
# J {160, 320, 800, 1600, 3200, 6400, 32000};

R = c(400, 200, 80, 40, 20, 10, 2, 800, 400, 160, 80, 40, 20, 4)
J = c(80, 160, 400, 800, 1600, 3200, 16000, 160, 320, 800, 1600, 3200, 6400, 32000)

## all pairs (R,J) were repeated 20 times. Each iteration generated estimates of statistics (of state mean, f1 mean, f2 mean).
state_mean_avgBias2 <- numeric(length(R)) # bias^2 of mean estimates (average over all sites & times)
state_mean_avgVariance <- numeric(length(R)) # variance of mean estimates (avearge over all sites & times)
state_mean_avgMSE <- numeric(length(R)) # mse of mean estimates (average over all sites & times)
f1_mean_avgBias2 <- numeric(length(R)) # bias^2 of mean estimates (average over all sites & times)
f1_mean_avgVariance <- numeric(length(R)) # variance of mean estimates (avearge over all sites & times)
f1_mean_avgMSE <- numeric(length(R)) # mse of mean estimates (average over all sites & times)

time_avg <- numeric(length(R))

for (i in 1:length(R)) {
    ## each particle filter with a (R,J) pair was repeated 20 times
    state_mean_sum <- matrix(0, nrow=ntimes, ncol=nvar)
    f1_mean_sum <- matrix(0, nrow=ntimes, ncol=nvar)
    state_mean_sqsum <- matrix(0, nrow=ntimes, ncol=nvar)
    f1_mean_sqsum <- matrix(0, nrow=ntimes, ncol=nvar)
    time_sum <- 0
    for(j in 1:20) {
        comstr <- paste("_K", nvar, "alpha", format(alpha, nsmall=2), "R", R[i], "J", J[i], "N", N, "_", (j-1), "_mixFewerResample.txt", sep='')
        state_mean=  as.matrix(read.table(paste("data/",directory,"/im_est_state_mean", comstr, sep=''))); colnames(state_mean) <- NULL
        f1_mean= as.matrix(read.table(paste("data/",directory,"/im_est_f1_mean", comstr, sep=''))); colnames(f1_mean) <- NULL
        time= as.matrix(read.table(paste("data/",directory,"/im_time", comstr, sep=''))); colnames(time) <- NULL
        time_sum <- time_sum + time
        state_mean_sum <- state_mean_sum + state_mean
        f1_mean_sum <- f1_mean_sum + f1_mean
        state_mean_sqsum <- state_mean_sqsum + state_mean^2
        f1_mean_sqsum <- f1_mean_sqsum + f1_mean^2
    }
    state_mean_bias = state_mean_sum / 20 - EXFILT
    f1_mean_bias = f1_mean_sum / 20 - (EXFILT^2+t(matrix(diag(SXXPOST), nrow=nvar)))
    state_mean_variance = state_mean_sqsum / 20 - (state_mean_sum / 20)^2
    f1_mean_variance = f1_mean_sqsum / 20 - (f1_mean_sum / 20)^2
    state_mean_avgBias2[i] = mean(state_mean_bias^2)
    state_mean_avgVariance[i] = mean(state_mean_variance)
    state_mean_avgMSE[i] = state_mean_avgBias2[i] + state_mean_avgVariance[i]
    f1_mean_avgBias2[i] = mean(f1_mean_bias^2)
    f1_mean_avgVariance[i] = mean(f1_mean_variance)
    f1_mean_avgMSE[i] = f1_mean_avgBias2[i] + f1_mean_avgVariance[i]
}




for(i in 1:1){
##pdf(paste('plots/',directory,'/performance_state_mean_K',nvar,'_',i,'_',date,'.pdf',sep=''))
par(mar=c(4,4,4,4)+.1)
xaxis = 1:(length(R)/2) + (length(R)/2)*(i-1)
plot(xaxis, state_mean_avgMSE[xaxis], type='b', xaxt='n', xlab='', ylab='', ylim=c(0,max(state_mean_avgMSE[xaxis])))
points(xaxis, state_mean_avgBias2[xaxis], type='b', lty=2)
points(xaxis, state_mean_avgVariance[xaxis], type='b', lty=3)
abline(h=0)
mtext('(#Particles per filter, #Parallel filters)', side=1, line=2)
mtext('MSE of Mean Estimates \n(average over all sites and times)', side=2, line=2)
title(main='Filter performance (dimension = 30)\n Average over 20 filtering repititions')
axis(side=1, at=xaxis, labels=sapply(xaxis, function(ii) paste('(',J[ii],',',R[ii],')', sep='')), cex.axis=.8)
legend('topright', legend=c('MSE', 'Bias^2', 'Variance'), lty=1:3)
dev.off()
}

for(i in 1:1){
##pdf(paste('plots/',directory,'/performance_f1_mean_K',nvar,'_',i,'_',date,'.pdf',sep=''))
par(mar=c(4,4,4,4)+.1)
xaxis = 1:7 + 7*(i-1)
plot(xaxis, f1_mean_avgMSE[xaxis], type='b', xaxt='n', xlab='', ylab='', ylim=c(0,max(f1_mean_avgMSE[xaxis])))
points(xaxis, f1_mean_avgBias2[xaxis], type='b', lty=2)
points(xaxis, f1_mean_avgVariance[xaxis], type='b', lty=3)
abline(h=0)
mtext('(#Particles per filter, #Parallel filters)', side=1, line=2)
mtext('MSE of X^2 Estimates \n(average over all sites and times)', side=2, line=2)
title(main='Filter performance (dimension = 30)\n Average over 20 filtering repititions')
axis(side=1, at=xaxis, labels=sapply(xaxis, function(ii) paste('(',J[ii],',',R[ii],')', sep='')), cex.axis=.8)
legend('topright', legend=c('MSE', 'Bias^2', 'Variance'), lty=1:3)
dev.off()
}





## computation of likelihood
require(mvtnorm)
mllik <- function(k) { # marginal log likelihood
    Yk <- Y[k + (1:ntimes - 1) * nvar]
    return(dmvnorm(Yk, sigma = MYCov, log = TRUE))
}
ml <- sapply(1:nvar, mllik)

#log likelihood
-ntimes*nvar/2*log(2*pi)-.5*log(det(SYY))-.5*Y%*%solve(SYY,Y)



##filename = paste(paste('mean_stf_1_K', nvar, sep=''), '.txt', sep='')
##stf_1_result = c(as.matrix(read.table(filename)))

filename = paste(paste('mean_stf_2_K', nvar, sep=''), '_np5000_G125000_10.txt', sep='')
stf_2_result = c(as.matrix(read.table(filename)))

filename = paste(paste('cov_stf_2_K', nvar, sep=''), '_np5000_G125000_10.txt', sep='')
stf_2_v_result = as.matrix(read.table(filename))

filename = paste(paste('mean_pf_K', nvar, sep=''), '_np100000.txt', sep='')
pf_result = c(as.matrix(read.table(filename)))

filename = paste(paste('cov_pf_K', nvar, sep=''), '_np100000.txt', sep='')
pf_v_result = as.matrix(read.table(filename))

filename = paste(paste('mean_if_K', nvar, sep=''), '_np5000.txt', sep='')
if_result = c(as.matrix(read.table(filename)))

filename = paste(paste('var_if_K', nvar, sep=''), '_np5000.txt', sep='')
if_v_result = c(as.matrix(read.table(filename)))


# Bias plot
plot(pf_result - EXlastPOST, type = 'b', col = 'blue', pch = 3, xlab='k', ylab='Bias')
#points(stf_1_result - EXlastPOST, type = 'b', pch = 3)
#plot(stf_1_result - EXlastPOST, type = 'b', pch = 3, lwd = 2, ylim = c(-3,3))
points(stf_2_result - EXlastPOST, type = 'b', pch = 4, lwd = 2)
points(if_result - EXlastPOST, type = 'b', pch = 3, col = 'green')
points(Y[(ntimes-1)*nvar + (1:nvar)] - EXlastPOST, type = 'p', col = 'red')
legend('topright',legend=c('STPF','PF','Independent','Measurement'),col=c('black','blue','green','red'),lty=1)
abline(h = 0, lty = 2)

# Variance plot
plot(diag(pf_v_result), type = 'b', col = 'blue', pch = 3, xlab='k', ylab='Variance')
points(diag(SX10POST),type='l')
points(diag(stf_2_v_result), type = 'b', pch = 4, lwd = 2)
points(if_v_result, type = 'b', pch = 3, col = 'green')
legend('topleft',legend=c('STPF','PF','Independent'),col=c('black','blue','green'),lty=1)

# MSE
stf_2_MSE <- (stf_2_result - EXlastPOST)^2 + diag(stf_2_v_result)
pf_MSE <- (pf_result - EXlastPOST)^2 + diag(pf_v_result)
if_MSE <- (if_result - EXlastPOST)^2 + if_v_result

plot(pf_MSE, type = 'b', col = 'blue', pch = 3, xlab='k', ylab='Mean Squared Error')
points(stf_2_MSE, type = 'b', pch = 4, lwd = 2)
points(if_MSE, type = 'b', pch = 3, col = 'green')
legend('topleft',legend=c('STPF','PF','Independent'),col=c('black','blue','green'),lty=1)


par(ask = TRUE)
for (t in 1:ntimes)
    plot(Y[(t-1)*nvar + (1:nvar)], type = 'b', ylim = c(-5,5))







## Likelihood function
loglikelihood <- function(theta) {
    nvar = 1
    ntimes = 50
    state0 = rep(0,nvar)#2 * (-(nvar-1)/2) : ((nvar-1)/2);

    alpha <- 0.10 # non-diagonal entries of the jump covariance matrix

    s = theta[1] # std dev of the marginal distribution of jump in one dimension
    d = theta[2] # measurement error standard deviation

    Sig = matrix(s^2*alpha, nvar, nvar); diag(Sig) <- s^2
    D = d^2 * diag(rep(1, nvar))

    TemporalCov = matrix(0, ntimes, ntimes)
    TemporalCov <- pmin(row(TemporalCov), col(TemporalCov))

    M <- kronecker(TemporalCov, Sig)

    COV = rbind(M,M)
    COV = cbind(COV,COV)

    for (t in 1:ntimes)
        COV[nvar*ntimes + (t-1)*nvar + (1:nvar), nvar*ntimes + (t-1)*nvar + (1:nvar)] = COV[nvar*ntimes + (t-1)*nvar + (1:nvar), nvar*ntimes + (t-1)*nvar + (1:nvar)] + D

    SYY = COV[nvar*ntimes + (1:(nvar*ntimes)), nvar*ntimes + (1:(nvar*ntimes))]
    
    ofilename = paste('data/obs_K', nvar, 'alpha', format(0.10, nsmall=2), '.txt', sep='')
    Y = c(t(as.matrix(read.table(ofilename))))[1:(ntimes*nvar)]

    ##log likelihood
    -ntimes*nvar/2*log(2*pi)-.5*log(det(SYY))-.5*Y%*%solve(SYY,Y)
}

sdvalues <- expand.grid(s=0.1*1:20, d=0.1*1:20)
ll <- matrix(apply(sdvalues, 1, loglikelihood), 20, 20)
persp(0.1*1:20, 0.1*1:20, exp(ll))
arrayInd(which.max(ll), dim(ll))




## create MIF plot
plot(0:(length(theta_hat)-1), theta_hat, xlim=c(0,length(theta_hat)-1),ylim=c(min(theta.vec),max(theta.vec)),type = 'b', pch = 3, axes=FALSE, xlab="", ylab="theta")
axis(2, ylim=c(min(theta.vec),max(theta.vec)))
axis(1, xlim = c(0,length(theta_hat)-1))
mtext("MIF iteration", side=1, line = 2)
box()
par(new=TRUE)
plot(lik.vec, theta.vec, type = 'l', xlim=c(min(lik.vec),max(lik.vec)), ylim=c(min(theta.vec),max(theta.vec)), axes=FALSE,ylab="",xlab="")
axis(3, xlim = c(min(lik.vec),max(lik.vec)))
mtext("Likelihood", side=3, line = 2)



