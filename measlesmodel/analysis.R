ntimes <- 832
ncity <- 20
cityname <- c("Halesworth","Lees","Mold","Dalton\n in.Furness","Oswestry","Northwich","Bedwellty","Consett","Hastings","Cardiff","Bradford","Hull","Nottingham","Bristol","Sheffield","Leeds","Manchester","Liverpool","Birmingham","London")
casedata.raw <- read.csv('20measles.csv', row.names=1)
startyear <- 1949
year <- (1:ntimes)*7/365.25 + startyear
startrow <- min(which(as.numeric(row.names(casedata.raw))>=1949))
casedata <- casedata.raw[startrow+0:(ntimes-1),]
thsize <- 17

dir <- 'data'; ncity=20
state <- read.table(paste(dir,'/state_ncity', ncity, '_ps9.txt', sep=''))[1+1:ntimes,]
ydata <- read.table(paste(dir,'/obs_ncity', ncity, '_ps9.txt', sep=''))[1:ntimes,]


date <- '170602'
ncity <- 20
R <- 5; J <- 4000; N <- ncity; options(scipen=6)
weighted <- TRUE; max_la <- 2; runif <- FALSE
case <- TRUE; # real data analysis?
M <- ifelse(runif, 8, 1)

lth <- function(G, n) {
    filenamecore <- paste('_K', ncity, ifelse(case, '_case_', ''), 'G', format(round(G,6),nsmall=6), ifelse(case, '', paste('rep', format(round(0.5,6),nsmall=6), sep='')), 'R', R, 'J', J, 'N', N, 'T', ntimes, 'if_', ifelse(runif,'T','F'), '_wi_', ifelse(weighted,'T','F'), '_la', max_la, '_', n, sep='') 
    ## compute likelihood
    lfilename <- paste('data/', date, '/im_est_likelihood', filenamecore, '.txt', sep='')
    if(file.exists(lfilename)) { lest <- as.matrix(read.table(lfilename))
    } else { print(paste("Warning: no lest file exits.", "G =",G, 'n =', n)); lest <- NA }
    ll <- apply(matrix(c(lest), nrow=ntimes), 2, mean, na.rm=TRUE) * ntimes # ll for each m
    nNaN <- apply(matrix(c(lest), nrow=ntimes), 2, function(x) sum(is.na(x)))
    if ((lenll <- length(ll)) < M) { ll[(lenll+1):M] <- NA }
    if ((lennn <- length(nNaN)) < M) { nNaN[(lennn+1):M] <- NA }
    if (runif) {
        thfilename <- paste('data/', date, '/im_theta', filenamecore, '.txt', sep='')
        if(file.exists(thfilename)) { th <- read.table(thfilename)
        } else { print(paste('Warning: no theta file exists.', "G =", G, 'n =', n)); th <- matrix(NA, ntimes*M, thsize) }
        if ((lenth <- dim(th)[1]) < M*ntimes) { th[(lenth+1):(M*ntimes),] <- NA }
        return(list(ll=ll, nNaN=nNaN, th=th))
    }
    return(list(ll=ll, nNaN=nNaN))
}


## IF2 plot (beta, G, dist_pow)
##pdf(paste('plots/', date, '/IF2', filenamecore, '.pdf', sep=''))

expand.grid2 <- function(...) { al <- list(...); y <- expand.grid(rev(al)); retval <- y[,dim(y)[2]:1]; colnames(retval) <- NULL; return(retval) }

nvec <- c(0:6); gvec <- c(50,100*c(1,3,5,7,9,11,15)); evalid <- c(0)

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

lasttheta <- t(sapply(1:length(thList), function(n) { lastline <- max(which(!is.na(thList[[n]][,1]))); sapply(estcomp, function(comp) ftrans[[comp]](thList[[n]][lastline,comp])) } ))
lasttheta <- cbind(lasttheta, sapply(1:length(thList), function(n) { lastline <- max(which(!is.na(ll[,n]))); ll[lastline,n] } ) )
colnames(lasttheta) <- c(compnames[estcomp], "ll")

pairs(lasttheta,col=color)



findRange <- function(comp,evall) { # find the range (ylim) for plotting
    concat <- c()
    for (n in evall) 
        concat <- c(concat, ftrans[[comp]](thList[[n]][,comp]))
    return(range(concat, na.rm=TRUE))
}

evall <- 1:length(thList)
ran <- sapply(1:thsize, findRange, evall=evall)

plotcomp <- function(comp, n, plot_add=FALSE) {
    par(new=plot_add)
    plot(NA, xlim=c(0,M), ylim=ran[,comp], type='l', ann=FALSE, xaxt='n')    
    points((1:(ntimes*M))/ntimes, ftrans[[comp]](thList[[n]][,comp]), type='l', col=color[n])
    mtext(text=compnames[comp], side=2, line=2, cex=1)
}

par(mfrow=c(ceiling((length(estcomp)+1)/2),2)); par(mar=c(0,3.5,0,1), oma=c(3.5,0,2,0), cex=1)
sapply(estcomp, function(comp) sapply((nn <- evall), function(n) plotcomp(comp,n,(n!=nn[1]))))

plot(NA, type='l', ann='FALSE', xlim=c(0,M), ylim=range(ll[,evall],na.rm=TRUE,finite=TRUE))#c(-52000,-49200))
for (n in evall)#1:(dim(ll)[2]))
    points(1:M, ll[,n], type='l', col=color[n])
mtext(text='likelihood', side=2, line=2, cex=1.3)
mtext(text='Iteration', side=1, line=2, cex=1.3)
mtext(text='Iterated filtering', outer=TRUE, cex=1.3)




## IVP estimation
case=TRUE; G=500; rep=0.5; R=1; J=4000; ncity=20; ntimes=3; N=ncity; runif=TRUE; weighted=TRUE; max_la = 2; M_ivp=60
filenamecore <- paste('_K', ncity, ifelse(case, '_case_', ''), 'G', format(round(G,6),nsmall=6), 'R', R, 'J', J, 'N', N, 'T', ntimes, 'if_', ifelse(runif,'T','F'), '_wi_', ifelse(weighted,'T','F'), '_la', max_la, '_', sep='')

ivpList <- list(); livp <- matrix(ncol=3, nrow=M_ivp)
for (n in 0:2) {
    ivpList[[n+1]] <- as.matrix(read.table(paste('data/170602/im_ivp', filenamecore, n, '.txt', sep='')))
    livp[,n+1] <- apply(matrix(unlist(read.table(paste('data/170602/im_ivp_est_likelihood', filenamecore, n, '.txt', sep=''))),nrow=3),2,sum)
}


plot(livp[,1],type='l')
points(livp[,2],type='l',col=2)
points(livp[,3],type='l',col=3)

plot(NA,xlim=c(0.00,0.002), ylim=c(0,20))
for (k in 1:ncity) {
    points(ivpList[[1]][ntimes*M_ivp,3*k-0],k,col=k)
    points(ivpList[[2]][ntimes*M_ivp,3*k-0],k,col=k)
    points(ivpList[[3]][ntimes*M_ivp,3*k-0],k,col=k)
    abline(h=k,lty=2)
}

plot(ivpList[[1]][,1],type='l',ylim=c(0,0.1))
for (k in 2:ncity)
    points(ivpList[[1]][,3*k-2],type='l',col=k)
for (k in 1:ncity)
    points(ivpList[[2]][,3*k-2],type='l',col=k,lty=1)
for (k in 1:ncity)
    points(ivpList[[3]][,3*k-2],type='l',col=k,lty=1)
    
plot(ivpList[[1]][,2],type='l',ylim=c(.0000,.002))
for (k in 2:ncity)
    points(ivpList[[1]][,3*k-1],type='l',col=k)
for (k in 1:ncity)
    points(ivpList[[2]][,3*k-1],type='l',col=k,lty=1)
for (k in 1:ncity)
    points(ivpList[[3]][,3*k-1],type='l',col=k,lty=1)
    
plot(ivpList[[1]][,3],type='l',ylim=c(.0000,.0020))
for (k in 2:ncity)
    points(ivpList[[1]][,3*k],type='l',col=k)
for (k in 1:ncity)
    points(ivpList[[2]][,3*k],type='l',col=k,lty=1)
for (k in 1:ncity)
    points(ivpList[[3]][,3*k],type='l',col=k,lty=1)
    




##legend('bottomleft', legend=gvec, col=1:length(gvec), lty=1)
##mtext(text=paste('K', ncity, 'T', ntimes, 'R', R, 'J', J, 'M', M), outer=TRUE, cex=1.3)



plot(param[,1],lasttheta[,dim(lasttheta)[2]],xlab='G',ylab='log lokelihood estimate with perturbed model')
mtext(side=3, text="mif (R0, nu_IR, nu_EI, alpha), No. city = 5 (4 iterations)")



####
ncity <- 5; N <- ncity; G <- 1500
thList <- list()
for(n in 0:2) {
    filenamecore <- paste('_K', ncity, ifelse(case, '_case_', ''), 'G', format(round(G,6),nsmall=6), 'mifR0nu_IR','R', R, 'J', J, 'N', N, 'T', ntimes, 'if_', ifelse(runif,'T','F'), '_wi_', ifelse(weighted,'T','F'), '_la', max_la, '_', n, sep='')
    lfilename <- paste('data/', date, '/im_est_likelihood', filenamecore, '.txt', sep='')
    lest <- as.matrix(read.table(lfilename))
    ll <- apply(matrix(c(lest), nrow=ntimes), 2, mean, na.rm=TRUE) * ntimes # ll for each m
    nNaN <- apply(matrix(c(lest), nrow=ntimes), 2, function(x) sum(is.na(x)))
    if ((lenll <- length(ll)) < M) { ll[(lenll+1):M] <- NA }
    if ((lennn <- length(nNaN)) < M) { nNaN[(lennn+1):M] <- NA }
    if (runif) {
        thfilename <- paste('data/', date, '/im_theta', filenamecore, '.txt', sep='')
        th <- read.table(thfilename)
        if ((lenth <- dim(th)[1]) < M*ntimes) { th[(lenth+1):(M*ntimes),] <- NA }
        thList[[n+1]] <- th
    }
}

concat <- c()
for (n in 1:(length(thList))) 
    concat <- c(concat, exp(thList[[n]][,1]))
ran <- range(concat)

thsname <- paste('_K', ncity, ifelse(case, '_case_', ''), 'G', format(round(G,6),nsmall=6), 'mifR0', 'R', R, 'J', J, 'N', N, 'T', ntimes, 'if_', ifelse(runif,'T','F'), '_wi_', ifelse(weighted,'T','F'), '_la', max_la, '_', sep='')

plot(NA, xlim=c(0,M), ylim=ran, type='l', xlab='Iteration', ylab='G')    
for (n in 0:2)
    points((1:(ntimes*M))/ntimes, exp(thList[[n+1]][,1]), type='l', col=n+1)
abline(h=20, lty=2)
plot(NA, xlim=c(0,832), ylim=c(10,2000), type='l', xlab='t', ylab='G', log='y')    
for (n in 0:4)
    points(1:ntimes, exp(thList[[n+1]][ntimes*3+1:ntimes,10]), type='l', col=n+1)
####




## expand.grid but first arg varying slowliest
expand.grid2 <- function(...) { al <- list(...); y <- expand.grid(rev(al)); retval <- y[,dim(y)[2]:1]; colnames(retval) <- NULL; return(retval) }
ll <- cbind((eg <- expand.grid2(gvec, reps, evalid)), apply(eg, 1, function(x) do.call("lth", as.list(x))$ll))


## plot likelihood evaluation
gs <- rep(gvec, each=length(reps))
profile <- cbind(param=gs, lp=unlist(ll.per.iter))
plot(gs, unlist(ll.per.iter))



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



## real data MCAP for G (take top 4 points for each G value)
ll4 <- c(apply(matrix(ll,nrow=7),2, function(x) sort(x, decreasing=TRUE)[1:4]))
gvec <- c(50,100*c(1,3,5,7,9,11,15)); reps <- 0:3
gs <- rep(gvec, each=length(reps))

##pdf('plots/170602/mcap_G.pdf')
cex=1.7
par(mar=c(3+(cex-1)*2,3+(cex-1)*2,0.5+(cex-1),0.5))
plot.profile(ll4,sqrt(gs),lambda=1.0)
axis(side=1, at=sqrt(gvec), labels=NA)
axis(side=1, at=sqrt(gvec), labels=gvec, tick=FALSE, line=-0.25+(cex-1), cex.axis=cex)
axis(side=2, labels=NA)
axis(side=2, tick=FALSE, line=-0.25+(cex-1), cex.axis=cex)
mtext(side=1, line=2+(cex-1), text='G', cex=cex)
mtext(side=2, line=2+(cex-1), text="log likelihood estimate", cex=cex)
##dev.off()



## simulated data MCAP for G (take top 3 points for each G value)
ll3 <- c(apply(matrix(ll,nrow=5),2, function(x) sort(x, decreasing=TRUE)[1:3]))
gvec <- c(50,100*c(1,3,5,7,9,11)); reps <- 0:2
gs <- rep(gvec, each=length(reps))
mcapresults <- mcap(ll3,gs,confidence=0.95,lambda=1.0,Ngrid=1000)

##pdf('plots/170502/mcap_sqrtG_lambda1.0.pdf')
par(mar=c(3,3,0.5,0.5))
plot.profile(ll3,sqrt(gs),lambda=1.0)
axis(side=1, at=sqrt(gvec), labels=NA)
axis(side=1, at=sqrt(gvec), labels=gvec, tick=FALSE, line=-0.25)
axis(side=2, labels=NA)
axis(side=2, tick=FALSE, line=-0.25)
mtext(side=1, line=2, text='G')
mtext(side=2, line=2, text="log likelihood estimate")
##dev.off()





##pdf(paste('plots/',date,'/prof_lik_G_case.pdf',sep=''))
plot.profile(lp=profile[,"lp"],parameter=profile[,"param"])
##dev.off()

##take 3 best points per g value
lp = c(apply(matrix(ll,nrow=5), 2, function(x) sort(x,decreasing=TRUE)[1:3]))
plot.profile(lp=lp, parameter=rep(gvec,each=3))

thetaswarm <- list()
for(n in reps) {
    thetaswarm[[n+1]] <- as.matrix(read.table(paste('data/', date, '/thetaswarm_m12', substr(filenamecore,start=1,stop=nchar(filenamecore)-2), '_', n, '.txt', sep='')))
}

plot(rep(1,60000), ftrans[[1]](thetaswarm[[4]][,14]))




### plot simulated data using the final MLE point found at data/170602/im_theta_K20_case_G300.000000R1J4000N20T832if_T_wi_T_la2_4.txt, but with varying values of G
dir <- 'data'; ncity=20
G <- c(0,100,321,1500)
state <- list(); ydata <- list()
for (i in 1:4) {
state[[i]] <- read.table(paste(dir,'/state_ncity', ncity, '_G', G[i], 'ps(G300_n4).txt', sep=''))[1+1:ntimes,]
ydata[[i]] <- read.table(paste(dir,'/obs_ncity', ncity, '_G', G[i], 'ps(G300_n4).txt', sep=''))[1:ntimes,]
}

##pdf('plots/170602/simulated_data_MLE_at_G300_n4.pdf', width=12, height=8.5)
cex=1
par(mfcol=c(4,3)); par(mar=c(1,2+(cex-1),0,0), oma=c(3+(cex-1),2.5+(cex-1),3+(cex-1),2))
color='red'
for(panel in 1:4) {
    city=20
    plot(year, casedata[,city], type='l', xaxt='n', yaxt='n', bty='n', col=color, lty=2)
    axis(side=1, labels=FALSE, at=c(1949,1950,1955,1960,1965))
    axis(side=2, labels=FALSE)
    axis(side=2, line=-.4+(cex-1), tick=FALSE, cex.axis=cex)
    points(year,ydata[[panel]][,city],type='l', lty=1)
    if (panel==1) {legend('topright', legend=c('simulated', 'real data') ,col=c(1, color), bg='white', lty=c(1,2), cex=cex) }
    if (panel==1) {mtext(side=3, line=.5+(cex-1), text=cityname[city], cex=cex)}
    if (panel==4) {axis(side=1, line=-.4+(cex-1), tick=FALSE, cex.axis=cex)}
    #mtext(side=4, line=3.8, text=paste('G=',G[panel],sep=''))
}
for(panel in 1:4) {
    city=10; col=4
    plot(year, casedata[,city], type='l', xaxt='n', yaxt='n', bty='n', col=color, lty=2)
    axis(side=1, labels=FALSE, at=c(1949,1950,1955,1960,1965))
    axis(side=2, labels=FALSE)
    axis(side=2, line=-.4+(cex-1), tick=FALSE, cex.axis=cex)
    points(year,ydata[[panel]][,city],type='l', lty=1)
    if (panel==1) {mtext(side=3, line=.5+(cex-1), text=cityname[city], cex=cex)}
    if (panel==4) {axis(side=1, line=-.4+(cex-1), tick=FALSE, cex.axis=cex); mtext(side=1, line=2+(cex-1)*2, text='year', cex=cex)}
}
for(panel in 1:4) {
    city=1; col=4
    plot(year, casedata[,city], type='l', xaxt='n', yaxt='n', bty='n', col=color, lty=2)
    text(x=par("usr")[2]+.1, y=mean(par("usr")[3:4]), srt = -90, xpd = NA, labels=paste('G=',G[panel],sep=''), cex=cex*1.5)
    axis(side=1, labels=FALSE, at=c(1949,1950,1955,1960,1965))
    axis(side=2, labels=FALSE)
    axis(side=2, line=-.4+(cex-1), tick=FALSE, cex.axis=cex)
    points(year,ydata[[panel]][,city],type='l', lty=1)
    if (panel==1) {mtext(side=3, line=.5+(cex-1), text=cityname[city], cex=cex)}
    if (panel==4) {axis(side=1, line=-.4+(cex-1), tick=FALSE, cex.axis=cex)}
}
mtext(side=2, line=0.5+(cex-1), outer=TRUE, text='weekly reported cases', cex=cex);
##dev.off()




###cityid <- 1
#par(mfrow=c(2,1))
#plot(year, state[1+1:ntimes, 9*5+2], type='l')
#plot(year, state[1+1:ntimes, 1*5+2], type='l', lty=2)

## likelihood estimates
G <- 0.00001
l <- rep(0,20)
for(n in 0:19) {
    lest <- read.table(paste('data/', date, '/im_est_likelihood_K', ncity, 'G', format(round(G,6),nsmall=6), 'R', R, 'J', J, 'N', N, 'T', ntimes, '_', n, '.txt', sep=''))
    l[n+1] <- sum(log(apply(lest, 1, mean)))
}



## plot estimated states ##
runno <- 1; # first, second, or ... mif run (runs can be continued)
G <- 500; no <- 0
filenamecore <- paste('_K20', ifelse(case,'_case_',''),'G',format(G,nsmall
=6), ifelse(case,'','rep0.500000'), 'exact','R5J4000N20T832if_F_wi_T_la2_', no, sep='')
date <- '170502'
dir <- 'data'; ncity=20
state <- read.table(paste(dir,'/state_ncity', ncity, '_ps9.txt', sep=''))[1+1:ntimes,]

ESSrawdata <- read.table(paste('data/', date, '/im_ESS', filenamecore, '.txt', sep=''))
filmean <- read.table(paste('data/', date, '/im_est_state_mean', filenamecore, '.txt', sep=''))[ntimes*(runno-1)+1:ntimes,]
filsqmean <- read.table(paste('data/', date, '/im_est_f1_mean', filenamecore, '.txt', sep=''))[ntimes*(runno-1)+1:ntimes,]
q10 <- read.table(paste('data/', date, '/im_est_q10', filenamecore, '.txt', sep=''))
q50 <- read.table(paste('data/', date, '/im_est_q50', filenamecore, '.txt', sep=''))
q90 <- read.table(paste('data/', date, '/im_est_q90', filenamecore, '.txt', sep=''))
filstd <- sqrt(filsqmean-filmean^2)
filstd <- apply(filstd, c(1,2), function(x) if(is.na(x)) 0 else x)
ui <- filmean+2*filstd
li <- filmean-2*filstd

 #The row of the ESS data file writes the ESS of the intermediate time steps for each island and each observation time. (ESS data file is (ntimes*R)-by-N matrix)
startrow <- (runno-1) * N*ntimes # if continued IF2 run, change startrow
ESS <- apply(ESSrawdata[startrow+1:(ntimes*N),], 1, sum)
plot(1:(N*ntimes), ESS, pch='.', log='y')

va <- 4; # variable to be plotted (1:S, 2:E, 3:I, 4:weekly cases, 5:pop)
maintitles <- c('Susceptible', 'Exposed', 'Infectious', 'Weekly Diagnosis/\nRecovery')
maintitle <- maintitles[va]
filetitles <- c('S', 'E', 'I', 'WeeklyIR')
filetitle <- filetitles[va]

#pdf(paste('plots/',date,'/',filetitle,'_K', ncity, 'G', format(round(G,6),nsmall=6), 'R', R, 'J', J, 'N', N, 'T', ntimes, '.pdf', sep=''), width=7.5, height=10)

##pdf('plots/170502/filterResults_trueParam.pdf', width=12, height=8.5)
par(mfcol=c(4,3), mar=c(1,2,0,0), oma=c(2.5,3,2,2))
cityids <- c(20,10,1)
for (cityid in cityids){
for (va in 1:4) {
column <- (cityid-1)*5+va # the column to be plotted
##plot(year, filmean[1:length(year),column], type='l', xaxt='n', xlab='', ylab='', main='')
plot(year, filmean[1:length(year),column], ylim=range(c(ui[,column], li[,column]), na.rm=TRUE), type='n', xaxt='n', yaxt='n', xlab='', ylab='', main='', bty='n')
axis(side=1, labels=FALSE, at=c(1949,1950,1955,1960,1965))
axis(side=2, labels=FALSE)
axis(side=2, line=-.4, tick=FALSE)
polygon(c(year, rev(year)), c(li[,column], rev(ui[,column])), col='grey80', border=NA)
points(year, filmean[1:length(year),column], type='l', lty=1)
points(year, q10[,column], type='l', lty=2, col=3)
points(year, q50[,column], type='l', lty=1, col=4)
points(year, q90[,column], type='l', lty=2, col=3)
points(year, state[1:length(year),column], type='l', col='red')
#points(year, casedata[1:length(year),cityid], type='l', col='red')
axis(side=1, labels=FALSE)
if(va==1) { mtext(text=cityname[20-ncity+cityid], side=3, line=.5) }
if(va==4 && cityid==1) {legend('topright', legend=c('Truth', 'Est. mean', 'Est. median', 'Est. 10%, 90% quantile'), bg='white', lty=c(1,1,1,2), col=c('red','black', 'blue', 'green'))}
if(va==4) {axis(side=1, line=-.4, tick=FALSE)}
if(cityid==20) { mtext(side=2, line=2, text=maintitles[va]) }
}
}
mtext(side=1, line=1, outer=TRUE, text='year');
##dev.off()



abline(v=year[588])
plot(c(t(ESS)), pch='.',  xaxt='n', xlab='', ylab='', ylim=c(1, R*J), log='y')
axis(side=1, labels=0:floor(ntimes/52), at=52*N*0:floor(ntimes/52), cex=1.3)
mtext(text='Year', side=1, line=2, cex=1.3)
mtext(text='ESS', side=2, line=2.5, cex=1.3)
mtext(text=maintitle, side=3, outer=TRUE, line=0.5, cex=1.3)

dev.off()



#### data plotting ###############
#pdf(paste('plots/',date,'/','ydata','_K', ncity, 'G', format(round(G,6),nsmall=6), 'T', ntimes, '.pdf', sep=''), width=7.5, height=3)

ydata <- read.table(paste('data/obs_ncity', ncity, '_ps5.txt', sep=''))[1:ntimes,]
par(mar=c(3.4,3.4,1,1))
plot(-1,1, xlim=c(0,max(year)), ylim=c(1,max(ydata)), xaxt='n', ann=FALSE, log='')
#plot(year, ydata[,20],type='l', col=1, xlim=c(-0.7,max(year)), xaxt='n', ann=FALSE, log='y')
for(cityid in 20:1){
    points(year, (ydata[,cityid]), type='l', col=21-cityid)
}
axis(side=1, at=0:floor(ntimes/52), cex=1.3)
mtext(side=1, line=2, text='year', cex=1.3)
mtext(side=2, line=2, text='Weekly cases', cex=1.3)
legend('topleft', legend=rev(cityname), col=1:20, lty=1, y.intersp=.5, ncol=2, cex=0.7)

dev.off()

######

state <- read.table(paste('data/state_ncity', ncity, '_ps5.txt', sep=''))
par(mar=c(3.4,3.4,1,1))
vc <- 1
plot(year, state[,5*20-5+vc],type='l', col=1, xlim=c(-0.7,max(year)), ylim=c(0,max(state[,5*20-5+vc])), xaxt='n', ann=FALSE, log='')
for(cityid in 19:1){
    #points(year, state[,5*cityid-5+vc], type='l', col=20-cityid+1)
}
axis(side=1, at=0:floor(ntimes/52), cex=1.3)
mtext(side=1, line=2, text='year', cex=1.3)
mtext(side=2, line=2, text='Weekly cases', cex=1.3)
legend('topleft', legend=rev(cityname), col=1:20, lty=1, y.intersp=.5, ncol=2, cex=0.7)


########



    plot(year, ydata[,12],type='l', col=1, xaxt='n', ann=FALSE, log='')
    points(year, ydata[,13],type='l', col=2, xaxt='n', ann=FALSE)


ydata <- read.table(paste('data/obs_ncity', ncity, '_ps4.txt', sep=''))[1:ntimes,]
par(mfrow=c(4,1))
par(mar=c(3.4,3.4,1,1))
for(i in 0:3){
    plot(year, ydata[,ncity-5*i],type='l', col=1, xaxt='n', ann=FALSE, log='')
    for(cityid in (ncity-5*i-1):(ncity-5*i-4)){
        points(year, ydata[,cityid], type='l', col=20-cityid)
    }
}


##pdf(paste('plots/',date,'/','I_Diagnosis_','_K', ncity, 'G', format(round(G,6),nsmall=6), 'R', R, 'J', J, 'N', N, 'T', ntimes, '.pdf', sep=''), width=7.5, height=10)

par(mfcol=c(6,2))
par(mar=c(0,0,0,0), oma=c(3,5,3,3), cex=1)
##
va <- 3; # variable to be plotted (1:S, 2:E, 3:I, 4:weekly cases, 5:pop)
maintitle <- ifelse(va==1, 'Susceptible Population', ifelse(va==2, 'Exposed Population', ifelse(va==3, 'Infectious Population', ifelse(va==4, 'Weekly Diagnosis/Recovery'))))
##
firstpanel <- TRUE
for(cityid in (ncity/5)*(5:1)) {
column <- (cityid-1)*5+va # the column to be plotted
plot(year, filmean[1:length(year),column], ylim=range(c(state[,column], ui[,column], li[,column])), type='n', xaxt='n', xlab='', ylab='', main='')
polygon(c(year, rev(year)), c(li[,column], rev(ui[,column])), col='grey80', border=NA)
points(year, filmean[1:length(year),column], type='l', lty=2)
points(year, state[1:length(year),column], type='l', col='red')
axis(side=1, labels=FALSE)
mtext(text=cityname[20-ncity+cityid], side=2, line=2.5, cex=1.3)
if(firstpanel) {legend('topleft', legend=c('Truth', 'Filter mean'), lty=c(1,2), col=c('red','black'))
mtext(text=maintitle, side=3, line=0.5, cex=1.3)
}
firstpanel <- FALSE
}
plot(c(t(ESS)), pch='.',  xaxt='n', xlab='', ylab='', ylim=c(1, R*J), log="y")
axis(side=1, labels=1:3, at=52*ncity*1:3, cex=1.3)
mtext(text='Year', side=1, line=2, cex=1.3)
mtext(text='Effective\n Sample Size', side=2, line=2.5, cex=1.3)
##
va <- 4; # variable to be plotted (1:S, 2:E, 3:I, 4:weekly cases, 5:pop)
maintitle <- ifelse(va==1, 'Susceptible Population', ifelse(va==2, 'Exposed Population', ifelse(va==3, 'Infectious Population', ifelse(va==4, 'Weekly Diagnosis/Recovery'))))
##
firstpanel <- TRUE
for(cityid in (ncity/5)*(5:1)) {
column <- (cityid-1)*5+va # the column to be plotted
plot(year, filmean[1:length(year),column], ylim=range(c(state[,column], ui[,column], li[,column])), type='n', xaxt='n', yaxt='n', xlab='', ylab='', main='')
polygon(c(year, rev(year)), c(li[,column], rev(ui[,column])), col='grey80', border=NA)
points(year, filmean[1:length(year),column], type='l', lty=2)
points(year, state[1:length(year),column], type='l', col='red')
axis(side=1, labels=FALSE)
axis(side=4)
if(firstpanel) {
mtext(text=maintitle, side=3, line=0.5, cex=1.3)
}
firstpanel <- FALSE
}
plot(c(t(ESS)), pch='.',  xaxt='n', yaxt='n', xlab='', ylab='', ylim=c(1, R*J), log="y")
axis(side=1, labels=1:3, at=52*ncity*1:3, cex=1.3)
axis(side=4)
mtext(text='Year', side=1, line=2, cex=1.3)



## plot likelihood estimate as a function of G
G <- 0.00001
l1 <- rep(0,20)
for(n in 0:19) {
    lest <- read.table(paste('data/', date, '/im_est_likelihood_K', ncity, 'G', format(round(G,6),nsmall=6), 'R', R, 'J', J, 'N', N, 'T', ntimes, '_', n, '.txt', sep=''))
    l1[n+1] <- sum(log(apply(lest, 1, mean)))
}
G <- 0.0001
l2 <- rep(0,20)
for(n in 0:19) {
    lest <- read.table(paste('data/', date, '/im_est_likelihood_K', ncity, 'G', format(round(G,6),nsmall=6), 'R', R, 'J', J, 'N', N, 'T', ntimes, '_', n, '.txt', sep=''))
    l2[n+1] <- sum(log(apply(lest, 1, mean)))
}
G <- 0.001
l3 <- rep(0,20)
for(n in 0:19) {
    lest <- read.table(paste('data/', date, '/im_est_likelihood_K', ncity, 'G', format(round(G,6),nsmall=6), 'R', R, 'J', J, 'N', N, 'T', ntimes, '_', n, '.txt', sep=''))
    l3[n+1] <- sum(log(apply(lest, 1, mean)))
}

#pdf(paste('plots/',date,'/', 'lest' ,'_K', ncity, 'G', format(round(G,6),nsmall=6), 'R', R, 'J', J, 'N', N, 'T', ntimes, '.pdf', sep=''))

par(mar=c(4,4,1,1))
plot(rep(1:3, each=20), c(l1,l2,l3), pch=3, ann=FALSE, xaxt='n', yaxt='n')
mtext(side=1, text='G', line=2.5, cex=1.3)
mtext(side=2, text='Likelihood estimates', line=2.5, cex=1.3)
axis(side=1, labels=c('1e-5', '1e-4', '1e-3'), at=c(1,2,3), cex.axis=1.3)
axis(side=2, cex.axis=1.3)
#text(x=2, y=min(l2)-100, labels='Success = 20/20', cex=1.3)
#text(x=1, y=min(l1,na.rm=TRUE)-100, labels='Success = 20/20', cex=1.3, adj=.1)
#text(x=3, y=max(l3,na.rm=TRUE)+100, labels='Success = 19/20', cex=1.3, adj=.9)

dev.off()


# read in old format ESS file
ESSfile <- file(paste('data/', date, '/im_ESS_K', ncity, 'G', format(round(G,6),nsmall=6), 'R', R, 'J', J, 'N', N, 'T', ntimes, '_', n, '.txt', sep=''))
open(ESSfile)
ESS <- matrix(0, nrow=(ntimes*N), ncol=R)
for (n in 1:(ntimes*N)) {
    for (r in 1:R) {
        line <- readLines(ESSfile, n=1)
        ESS[n,r] <- as.numeric(substr(line, start=regexpr("ESS : ", line)[1]+6, stop=regexpr(",   CV", line)[1]-1))
        if (is.na(ESS[n,r])) ESS[n,r] = 0
    }
}
close(ESSfile)
