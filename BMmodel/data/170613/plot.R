nvar <- 20
im <- read.table('im_states.txt')
obs <- read.table('../obs_K20alpha0.00.txt')

J <- 2000; draw <- 400 # how many particles will be plotted
d1 <- matrix(unlist(im[,1]), nrow=J)#[J/draw*1:draw,]  # the first component
d2 <- matrix(unlist(im[,2]), nrow=J)#[J/draw*1:draw,]  # the second component

library(plotrix)

S <- 20
##dev.new(width=9,height=3.3)
pdf(file='../../plots/170613/im_states.pdf', width=9, height=3.3)
par(mar=c(2,2,0,0.2), oma=c(0,2,0,0))
par(mfrow=c(1,3),pty='s')
for (s in c(4,12,20)) {
    plot(d1[,s],d2[,s],xlim=c(-4,2),ylim=c(-3,3), pch='.', xlab='', ylab='', yaxt='n', bty='n')
    if (s==4) { mtext(text=expression(x^2), side=2, line=2) }
    axis(side=2, labels=ifelse(s==4, TRUE, FALSE))
    mtext(text=expression(x^1), side=1, line=2.5)
    label = ifelse(s==4, 'A', ifelse(s==12, 'B', 'C'))
    mtext(text=paste(label,'.  s=',s,sep=''), side=3, line=0)
    mean = (1+s/S)/3*unlist(obs[1,1:2])
    var = 1+s/S - (1+s/S)^2/3 # 1/3*(1+s/S)*(2-s/S)
    points(0,0, pch='O', col='green', cex=2, lwd=2)
    points(mean[1],mean[2], pch=4, col='red', cex=2, lwd=3)
    points(obs[1,1], obs[1,2], pch=2, col='purple', cex=2, lwd=2)
    draw.circle(x=mean[1], y=mean[2], radius=sqrt(var*qchisq(0.95,df=2)), lty=2, border='blue', lwd=2)
    #segments(0,0,obs[1,1],obs[1,2], lty=2)
}
dev.off()

