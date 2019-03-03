## State inference
LMSE <- function(nvar, alpha, diagCov=FALSE) {
    ## for the linear Gaussian model, compute the exact likelihood and filtering distributions
    ## Then compare with the estimated values from the GIRF and compute the mean squared errors (MSE)
    ## alpha: non-diagonal entries of the jump covariance matrix
    ## diagCov: was diagonal covariance matrix used for computing predictive likelihood?
    ntimes = 50
    state0 = rep(0,nvar) #2 * (-(nvar-1)/2) : ((nvar-1)/2);

    s = 1.0 # std dev of the marginal distribution of jump in one dimension
    d = 1.0 # measurement error standard deviation

    Sig = matrix(s^2*alpha, nvar, nvar); diag(Sig) <- s^2
    D = d^2 * diag(rep(1, nvar))

    ofilename = paste('../obs_K', nvar, 'alpha', format(alpha, nsmall=2), '.txt', sep='')
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
    
    date = '171019'
    directory = date

    le <- numeric(); sum_X50_mean <- numeric(nvar); sum_X50_mean_sq <- numeric(nvar)
    MSEE <- numeric()
    rep <- 40
    for (n in 1:rep) {
        R = 60; J = 1000; N = nvar; options(scipen=7)
        commonString = paste("_K", nvar, "alpha", format(alpha, nsmall=2), "R", R, "J", J, "N", N, "T", ntimes, "if_", "F", "wi_", "T", "_la", 2, "_", (n-1), ifelse((!diagCov) || (alpha==0), ".txt", "_diag.txt"), sep='')
        ##commonString = paste("_K", nvar, "alpha", format(alpha, nsmall=2), "R", R, "J", J, "N", N, "_1_unsmoothed.txt", sep='')
        lestimate <- read.table(paste("im_est_likelihood", commonString, sep=''))
        lestimate <- lestimate[dim(lestimate)[1]-(ntimes-1):0,]
        le[n] <- sum(lestimate)
        paste('exact log likelihood =', ll, 'likelihood estimate =', le)
        ##
        state_mean= as.matrix(read.table(paste("im_est_state_mean", commonString, sep=''))); colnames(state_mean) <- NULL
        sum_X50_mean <- sum_X50_mean + state_mean[50,]
        sum_X50_mean_sq <- sum_X50_mean_sq+ state_mean[50,]^2
        ##f1_mean= as.matrix(read.table(paste("data/",directory,"/im_est_f1_mean", commonString, sep=''))); colnames(f1_mean) <- NULL
        ##state_var=(R*J)/(R*J-1)*(f1_mean - state_mean^2)
        MSEE[n] <- sum((EXFILT[50,]-state_mean[50,])^2)/nvar
    }

    MSE <- c(k=nvar, alpha=alpha, bias_sq = (bias_sq <- sum((EXFILT[50,]-sum_X50_mean/rep)^2)/nvar), var = (variance <- sum((sum_X50_mean_sq - sum_X50_mean^2/rep)/rep)/nvar), MSE = bias_sq+variance, seMSE = sqrt(var(MSEE)/rep)) # seMSE = standard error of the MSE
    L <- c(k=nvar, alpha=alpha, exact_l = ll, log_mean_est_l = mean(le)+log(mean(exp(le-mean(le)))), se_log_mean_est_l = sqrt(var(exp(le-mean(le)))/rep)/mean(exp(le-mean(le))), mean_log_est_l = mean(le), se_log_est_l = sqrt(var(le)/rep) )

    return(list(MSE = MSE, L = L))
}

#MSE_dim <- sapply(c(20,50,100,200), function(d) LMSE(d, alpha=0.00, diagCov=FALSE)$MSE) ## saved as MSE_dim.txt 
#L_dim <- sapply(c(20,50,100,200), function(d) LMSE(d, alpha=0.00, diagCov=FALSE)$L) 
MSE_alpha_exact20 <- sapply(0.1*0:5, function(a) LMSE(20, a, diagCov=FALSE)$MSE)
L_alpha_exact20 <- sapply(0.1*0:5, function(a) LMSE(20, a, diagCov=FALSE)$L) 
MSE_alpha_diag20 <- sapply(0.1*0:5, function(a) LMSE(20, a, diagCov=TRUE)$MSE)
L_alpha_diag20 <- sapply(0.1*0:5, function(a) LMSE(20, a, diagCov=TRUE)$L) 
MSE_alpha_exact50 <- sapply(0.1*0:5, function(a) LMSE(50, a, diagCov=FALSE)$MSE)
L_alpha_exact50 <- sapply(0.1*0:5, function(a) LMSE(50, a, diagCov=FALSE)$L) 
MSE_alpha_diag50 <- sapply(0.1*0:5, function(a) LMSE(50, a, diagCov=TRUE)$MSE)
L_alpha_diag50 <- sapply(0.1*0:5, function(a) LMSE(50, a, diagCov=TRUE)$L) 
#write.table(MSE_alpha_exact20, "MSE_alpha_exact20.txt")
#write.table(L_alpha_exact20, "L_alpha_exact20.txt")
#write.table(MSE_alpha_diag20, "MSE_alpha_diag20.txt")
#write.table(L_alpha_diag20, "L_alpha_diag20.txt")
#write.table(MSE_alpha_exact50, "MSE_alpha_exact50.txt")
#write.table(L_alpha_exact50, "L_alpha_exact50.txt")
#write.table(MSE_alpha_diag50, "MSE_alpha_diag50.txt")
#write.table(L_alpha_diag50, "L_alpha_diag50.txt")


## the results from above commands are saved as text files. 
#MSE_dim <- as.matrix(read.table('MSE_dim.txt', row.names=1))
#L_dim <- as.matrix(read.table('L_dim.txt', row.names=1))
MSE_alpha_exact20 <- as.matrix(read.table('MSE_alpha_exact20.txt', row.names=1))
L_alpha_exact20 <- as.matrix(read.table('L_alpha_exact20.txt', row.names=1))
MSE_alpha_diag20 <- as.matrix(read.table('MSE_alpha_diag20.txt', row.names=1))
L_alpha_diag20 <- as.matrix(read.table('L_alpha_diag20.txt', row.names=1))
MSE_alpha_exact50 <- as.matrix(read.table('MSE_alpha_exact50.txt', row.names=1))
L_alpha_exact50 <- as.matrix(read.table('L_alpha_exact50.txt', row.names=1))
MSE_alpha_diag50 <- as.matrix(read.table('MSE_alpha_diag50.txt', row.names=1))
L_alpha_diag50 <- as.matrix(read.table('L_alpha_diag50.txt', row.names=1))

##pdf('MSE_dim(alpha0.00R60J1000).pdf')
##pdf('../../report/JRSSB/figures/MSE_dim(alpha0.00R60J1000).pdf', width=8, height=5)
#par(mar=c(3.5,5.5,1,1))
#plot(MSE_dim["k",], MSE_dim["MSE",], type='b', ann=FALSE, xlim=c(20,200), ylim=c(0,0.01), xaxt='n', yaxt='n', cex=1.6)
#points(MSE_dim["k",], MSE_dim["bias_sq",], type='b', lty=2, pch=2, cex=1.6)
#abline(h=0)
##legend('topleft', legend=c('MSE', expression(bias^2)), lty=c(1,2), pch=c(1,2), cex=1.3)
#axis(side=1, at=c(0,20,50,100,200), cex.axis=1.6)
#axis(side=2, cex.axis=1.6)
#mtext('Dimension', side=1, line=2.5, cex=1.6)
#mtext('Mean squared error \n(average over all sites)', side=2, line=2.5, cex=1.6)
##dev.off()

##pdf('MSE_alpha(dim20R60J1000).pdf')
##pdf('../../report/JRSSB/figures/MSE_alpha(dim20R60J1000).pdf', width=8, height=5)
par(mar=c(3.5,6,1,1))
plot(MSE_alpha_diag20["alpha",], MSE_alpha_diag20["MSE",], type='b', lty=5, pch=4, ann=FALSE, xlim=c(0,0.5), ylim=c(0,0.028), cex.axis=1.8, cex=1.6)
arrows(MSE_alpha_diag20["alpha",], MSE_alpha_diag20["MSE",]-MSE_alpha_diag20["seMSE",], MSE_alpha_diag20["alpha",], MSE_alpha_diag20["MSE",]+MSE_alpha_diag20["seMSE",], length=0.05, angle=90, code=3, lty=5)
points(MSE_alpha_exact20["alpha",], MSE_alpha_exact20["MSE",], type='b', cex=1.6)
arrows(MSE_alpha_exact20["alpha",], MSE_alpha_exact20["MSE",]-MSE_alpha_exact20["seMSE",], MSE_alpha_exact20["alpha",], MSE_alpha_exact20["MSE",]+MSE_alpha_exact20["seMSE",], length=0.05, angle=90, code=3)
abline(h=0)
#legend('topleft', legend=c('Diagonal covariance used', 'Exact covariance used'), lty=c(5,1), pch=c(4,1), cex=1.3)
mtext('Correlation coefficient', side=1, line=2.5, cex=1.8)
mtext('Mean squared error\n (average over all sites)', side=2, line=2.8, cex=1.8)
##dev.off()

##pdf('MSE_alpha(dim50R60J1000).pdf')
##pdf('../../report/JRSSB/figures/MSE_alpha(dim50R60J1000).pdf', width=8, height=5)
par(mar=c(3.5,2.5,1,1))
plot(MSE_alpha_diag50["alpha",], MSE_alpha_diag50["MSE",], type='b', lty=5, pch=4, ann=FALSE, xlim=c(0,0.5), ylim=c(0,0.026), cex.axis=1.8, cex=1.6)
arrows(MSE_alpha_diag50["alpha",], MSE_alpha_diag50["MSE",]-MSE_alpha_diag50["seMSE",], MSE_alpha_diag50["alpha",], MSE_alpha_diag50["MSE",]+MSE_alpha_diag50["seMSE",], length=0.05, angle=90, code=3)
points(MSE_alpha_exact50["alpha",], MSE_alpha_exact50["MSE",], type='b', cex=1.6)
arrows(MSE_alpha_exact50["alpha",], MSE_alpha_exact50["MSE",]-MSE_alpha_exact50["seMSE",], MSE_alpha_exact50["alpha",], MSE_alpha_exact50["MSE",]+MSE_alpha_exact50["seMSE",], length=0.05, angle=90, code=3)
abline(h=0)
##legend('topleft', legend=c('Diagonal covariance used', 'Exact covariance used'), lty=c(5,1), pch=c(4,1), cex=1.3)
mtext('Correlation coefficient', side=1, line=2.5, cex=1.8)
##mtext('Mean squared error\n (average over all sites)', side=2, line=2.5, cex=1.8)
##dev.off()
