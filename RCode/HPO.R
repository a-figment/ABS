setwd("~/Dropbox/ABS/Studies/S-HPO/Data/")


a  <- seq(pi/30, 30*pi/31, length=1000)
a  <- seq(pi/20, 20*pi/21, length=1000)
a  <- seq(pi/2, 30*pi/31, length=2000)
a  <- seq(2*pi/3, 20*pi/21, length=1000)
dy <- (1/2)*(1/tan(a/2))
d <- dy 
h0 <- ReadDataList(lf("Scatter"), idx=1)
vx <- ReadDataList(lf("Scatter"), idx=2)
vy <- ReadDataList(lf("Scatter"), idx=3)
a  <- ReadDataList(lf("Scatter"), idx=4)*180/pi
th <- atan2(vy,vx)
th <- th*180/pi

phi <- function(alpha, d, y0) { 
    id <- which(y0 > 0 & y0 < d*cos(alpha))
    res <- c()
    if (length(id) == 0) {
	res <- rep(0, length(alpha))
    } else {
	res[id] <- 180 - alpha[id]
	res[-id] <- 0;
    }
    res
#    for (i in seq(length(alpha))) {
#	if (y0[i] > 0 & y0[i] < d[i]*cos(alpha[i])) {
#	    result[i] <- 180 - alpha[i]
#	} else {
#	    result[i] <- 0
#	}
#    }
#    return(res)
}

mplot(h0[,], th, type='p', cex=0.2, pch='.',col='black') #ylim=c(0,1))

res <- list.cbind.na(sapply(1:ncol(a), function(x) phi(a[,x], d[,x], h0[,x])))
mplot(h0, (res), col='red', type='p', pch=19, cex=0.1, add=TRUE)
eplot(tpose(h0[,]), A, type='p', cex=0.2, pch='.',col='black') #ylim=c(0,1))
A <- a * 180/pi
lines(th,a)
mplot(th[,], a, type='p', cex=3.5, pch='.',col='black') #ylim=c(0,1))
mplot(h0[,], a, type='p', cex=3.5, pch='.',col='black') #ylim=c(0,1))


mplot(a, 90-atan2(vy,vx)*180/pi, type='p', cex=0.2, pch='.')
mplot(h0[,], atan(vy/vx)*180/pi, type='p', cex=0.2, pch='.', xlim=c(0,1))
mplot(th, h0[,], type='p', cex=0.1, pch='.', col='black')
mplot(a*180/pi, th[,], type='p', cex=0.1, pch='.', col='black')
lines(th[,1], pi-th[,1])

pos <- data.matrix(fread(lf("Pos")))
x <- pos[,1]
y <- pos[,2]

id <- 30
PCB(1/2, d[id], dy[id], from=-1, to=1, asp=FALSE)
n <- 15000:15201

mplot(-1/2, h0[1,1], type='l', asp=1, add=TRUE)
#points(x[n[1]], y[n[1]], pch=19, col='green')


t <- seq(0,1e5,length=1001)[-1]
#pdf("Vbal-Distribution.pdf", width=20, height=20)
par(mfrow=c(2,2), bg='grey', cex=2)
j <- 80
dx <- (ReadDataList(lf("DX"), colidx=j)[,])
hist(as.vector(abs(dx[]) / t[j]), breaks=100, freq=F, ylim=c(0,10), xlim=c(0,1), 
     main=TeX("$\\alpha = 2\\pi / 3$"), 
     xlab=TeX("$V_{bal}$"))

setwd("~/Dropbox/ABS/Studies/S-HPO/Data/")
t <- seq(1,1e5+1,length=20001)
t <- exp(seq(0,log(1e5),length=5000))
t <- seq(1,1e4,length=10000)
id <- 11
eps <- pi/(10000)
ang <- sort(c(0,seq(-eps, eps, length=10000)))
ratio <- seq(1/4, 1, length=21)
dy <- 1/2
d <- ratio * dy
vb <- function(n) { 1 - (1 / (1 + n*ratio)) } # IH
vb <- 1/(1 + ratio[id])

TSD <- tpose(fread("p0-a2-r11--TSD.dat"))
DT  <- tpose(fread("p0-a2-r10--DT.dat"))


spar(c(2,1))
s <- seq(150)
s <- seq(1)
k <- rindex(ncol(TSD), 0)
matplot(t[-s], abs(TSD)[-s,k] * t[-s],
	type='l', lty='solid', lwd=2, 
	#ylim=c(0.55,0.65), 
	#xlim=c(0,1000),
	#log='x',
	ylab=TeX("$\\Delta x / t$"), 
	xlab=TeX("$t$"),
	col=rainbow(10))

matplot(t[-s][-1], diff(abs(TSD)[-s,k] * t[-s]),
	type='l', lty='solid', lwd=2, 
	#ylim=c(0.55,0.65), 
	#xlim=c(0,1000),
	#log='x',
	ylab=TeX("$\\Delta x / t$"), 
	xlab=TeX("$t$"),
	col=rainbow(10))

thresh <- 0.005
abline(h=vb + c(thresh, -thresh))

idx <- apply(abs(TSD)[-s,], 2, function(x) min(which(x < (vb - thresh))))
abline(v=t[-s][idx])
abline(v=DT[,1], col='red')


#' Distribution of departure times
DT  <- ReadData(lf("DT.dat"))[,1,]
spar()
N <- 1-seq(nrow(DT))/nrow(DT)
S <- apply(DT, 2, sort)

pdf("DTDist.pdf", width=18, height=18) 
{
    spar()
    matplot(S, N, 
    	type='s', log='xy',
    	ylab=TeX("$1 - F(\\tau)$"),
    	xlab=TeX("$\\tau$ - Departure time"),
    	col=rainbow(10),
    	lty='solid',
    	lwd=3)
    LF <- LinFit(mat=N, time=S, range=seq(5000,9000), logs=TRUE)
    lines(S[,1], exp(1.1*LF$b[1])*S[,1]^(LF$ax[1]), lwd=2.2, lty='dashed')
    
    legend('bottomleft', 
            legend=TeX(paste0("$r = ", ratio[11:13], "$")), 
            col=rainbow(10), lty='solid', lwd=2)
    legend('topright', 
            legend=TeX(paste0("$1 - F(\\tau) \\sim \\tau^{", fm(LF$ax[1], 3), "}$")), 
            col='black', lty='dashed', lwd=2)
    
}; dev.off() 

#' Checking trajectories
{
    id <- 3
    a <- paste0("a", seq(6,8))
    d <- c(2.16505, .3, 1.21244)
    alpha <- c(pi/3, pi/2, 2*pi/3)
    dy <- (1/2) * (1/tan(alpha/2))
    V <- cbind(c(-1/2, 0, 1/2, 1/2, 0, -1/2), c(0, dy[id], 0, d[id], dy[id] + d[id], d[id]))
    from <- 5750
    to <- 5755
    PolygonBaseplate(V[,1], V[,2], from, to)
    PlotChannelBoundary(V[,1], V[,2], from, to)
    matplot(Traj.x, Traj.y, type='l', add=TRUE)
}
    Traj.x <- tpose(fread(lff(a[id], "xPos")))
    Traj.y <- tpose(fread(lff(a[id], "yPos")))
}
#fwrite(x=format(as.data.frame(cbind(Traj.x[,1], Traj.y[,1], Traj.x[,2], Traj.y[,2])),digits=5), 
       file="2pi3--HPO.tkz", sep=",", col.names=FALSE)


ct <- cumsum(sqrt(diff(Traj$V1)^2 + diff(Traj$V2)^2))
lines(ct, abs(Traj$V1[-1] / ct), col='red')



pdf("PPSample.pdf", width=20, height=15) 
spar()
{
    eps <- pi/(10000)
    ang <- sort(c(0,seq(-eps, eps, length=100)))
    K <- sort(c(seq(1,ncol(TSD),by=12), 51))
    matplot(t[-s], abs(TSD)[-s,K],
    	type='l', lty='solid', lwd=2, 
    	main=TeX(paste0("$r = ", ratio[id], ", \\alpha = \\pi/2$")),
    	log='x',
    	col=rainbow(10),
    	ylab=TeX("$|\\Delta x | / t$"), 
    	xlab=TeX("$t$"))
    legend('bottomleft', 
           legend=TeX(paste0("$\\theta_0 = ", fm(ang[K],3), "$")), 
           col=rainbow(10), lty='solid', lwd=2)
    thresh <- 0.005
    abline(h=vb -thresh, lty='dashed', lwd=2)
    
    legend(x=300, y=0.1, 
           legend="Threshold",
           col='black', lty='dashed', lwd=2)
}; dev.off() 



