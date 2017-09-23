ntimes <- 832
ncity <- 20
cityname <- c("Halesworth","Lees","Mold","Dalton\n in.Furness","Oswestry","Northwich","Bedwellty","Consett","Hastings","Cardiff","Bradford","Hull","Nottingham","Bristol","Sheffield","Leeds","Manchester","Liverpool","Birmingham","London")
casedata.raw <- read.csv('../../20measles.csv', row.names=1)
startyear <- 1949
year <- (1:ntimes)*7/365.25 + startyear
startrow <- min(which(as.numeric(row.names(casedata.raw))>=1949))
casedata <- casedata.raw[startrow+0:(ntimes-1),]
thsize <- 17

dir <- 'data'; ncity=20
state <- read.table(paste('../../',dir,'/state_ncity', ncity, '_ps9.txt', sep=''))[1+1:ntimes,]
ydata <- read.table(paste('../../',dir,'/obs_ncity', ncity, '_ps9.txt', sep=''))[1:ntimes,]


date <- '170502'
ncity <- 20
R <- 5; J <- 4000; N <- ncity; options(scipen=6)
weighted <- TRUE; max_la <- 2; runif <- FALSE
case <- FALSE; # real data analysis?
M <- ifelse(runif, 8, 1)

lth <- function(G, n) {
    filenamecore <- paste('_K', ncity, ifelse(case, '_case_', ''), 'G', format(round(G,6),nsmall=6), ifelse(case, '', paste('rep', format(round(0.5,6),nsmall=6), sep='')), 'R', R, 'J', J, 'N', N, 'T', ntimes, 'if_', ifelse(runif,'T','F'), '_wi_', ifelse(weighted,'T','F'), '_la', max_la, '_', n, sep='') 
    ## compute likelihood
    lfilename <- paste('../../data/', date, '/im_est_likelihood', filenamecore, '.txt', sep='')
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

nvec <- c(0:4); gvec <- c(50,100*c(1,3,5,7,9,11)); evalid <- c(0)

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


## simulated data MCAP for G (take top 3 points for each G value)
ll3 <- c(apply(matrix(ll,nrow=5),2, function(x) sort(x, decreasing=TRUE)[1:3]))
gvec <- c(50,100*c(1,3,5,7,9,11)); reps <- 0:2
gs <- rep(gvec, each=length(reps))
mcapresults <- mcap(ll3,gs,confidence=0.95,lambda=1.0,Ngrid=1000)

##pdf('plots/170502/mcap_sqrtG_lambda1.0.pdf')
pdf('../../../report/JRSSB/figures/simulated_data_mcap_sqrtG_lambda1.0.pdf', width=8, height=5)
cex=1.6
par(mar=c(3+(cex-1)*2,3+(cex-1)*2,0.5+(cex-1),0.5))
plot.profile(ll3,sqrt(gs),lambda=1.0)
axis(side=1, at=sqrt(gvec), labels=NA)
axis(side=1, at=sqrt(gvec), labels=gvec, tick=FALSE, line=-0.25+(cex-1), cex.axis=cex)
axis(side=2, labels=NA)
axis(side=2, tick=FALSE, line=-0.25+(cex-1), cex.axis=cex)
mtext(side=1, line=2+(cex-1)*1.5, cex=cex, text='G')
mtext(side=2, line=2+(cex-1)*1.5, cex=cex, text="log likelihood estimate")
dev.off()
##dev.off()


