
# Set the WD
setwd("~/Dropbox/ABS/CppCode/DATA/")
setwd("~/Dropbox/ABS/Studies/S-DiskChannel/Data/")
setwd("~/Dropbox/ABS/Studies/S-BallisticFront/Data/")
setwd("~/Dropbox/ABS/Studies/S-HPO/Data/")
setwd("~/Dropbox/PHD/Investigations/FrontSpeed/Data/PI4/FH")

Config <- function(idx,a) {
    alpha <<- a
    dy <<- 1/(2*tan(a/2))
    d <<- ratio[idx + 1] * dy
    V <<- getChannelVertices(dx = 1/2, alpha = a, d = d)
}

alpha <- pi/2
dy <- 1/(2*tan(alpha/2))
did <- 1
ratio <- seq(3/10, 1, length=11)    # pi/4 FH
ratio <- seq(1/4, 2, length=20)	    # pi/3
ratio <- seq(1/4, 1, length=10)

ratio <- 2^seq(log2(1/4), log2(10), length=20)
ratio <- 2^seq(log2(1/4), log2(2), length=20)

d <- ratio[did + 1] * dy
d <- 1/2


V <- getChannelVertices(dx = 1/2, alpha = alpha, d = 1/2)
        PolygonBaseplate(V[,1], V[,2], from, to)
        PlotChannelBoundary(V[,1], V[,2], from, to)

        lines(Traj.x[,nRows[i]], Traj.y[,nRows[i]], lwd=0.75, col='blue')

Traj.x <- ReadDataList(lf(paste0(did, "--ConTraj-")), colidx=1)
Traj.y <- ReadDataList(lf(paste0(did, "--ConTraj-")), colidx=2)

Traj.x <- ReadDataList(lf(paste0(did, "--ColTraj-")), colidx=1)
Traj.y <- ReadDataList(lf(paste0(did, "--ColTraj-")), colidx=2)

Traj.x <- ReadDataList(lf(paste0(did, "--DiscTraj_")), colidx=1)
Traj.y <- ReadDataList(lf(paste0(did, "--DiscTraj_")), colidx=2)


#' ---------------------------------------------------------------------
#' Plot all pdf's
from  <- 0
to    <- 2
nPerPage <- 3

pdf(paste0("d", did, "--Saturated-Traj.pdf"), width=20, height=10)
pdf(paste0("r", did, "--Traj.pdf"), width=20, height=10)
for (j in seq(floor(ncol(Traj.x) / nPerPage))) { 

    nRows <- seq((j-1)*nPerPage + 1, j * nPerPage)
    
    #pdf(paste0("d", j, "--SatTraj.pdf"), width=14, height=14)
    par(mfrow=c(length(nRows),1))

    for (i in seq(length(nRows))) {
    
        #s <- round(seq(which.min(Traj.x[,nRows[i]]), which.max(Traj.x[,nRows[i]])))

        PolygonBaseplate(V[,1], V[,2], from, to)
        PlotChannelBoundary(V[,1], V[,2], from, to)

        lines(Traj.x[,nRows[i]], Traj.y[,nRows[i]], lwd=0.75, col='blue')
        points(Traj.x[1,nRows[i]], Traj.y[1,nRows[i]], col='red', pch=19, cex=1)
    
        lines(Traj.x[1:10,], Traj.y[1:10,], lwd=0.75, col='blue')
        points(Traj.x[1,nRows[i]], Traj.y[1,nRows[i]], col='red', pch=19, cex=1)
    }
}



#' --------------------------- pANEL ---------------------------------


pdf(paste0("2PI3--Traj.pdf"), width=20, height=15)
par(mfrow=c(2,1), oma=rep(0,4), mar=rep(0,4))
dids <- 12 
Config(dids, 2*pi/3)
Config(dids, pi/4)
V

ratio <- 2^seq(log2(1/4), log2(10), length=20)
ratio <- 2^seq(log2(1/4), log2(2), length=20)

for (i in seq(length(dids))) {
    Traj.x <- ReadDataList(lf(paste0("a3-r", dids[i], "--ColTraj-")), colidx=1)
    Traj.y <- ReadDataList(lf(paste0("a3-r", dids[i], "--ColTraj-")), colidx=2)
    
    PlotChannelBoundary(V[,1], V[,2], -16, 16)
    
    s <- seq(1,ncol(Traj.x),by=4)
    s <- 50
    matplot(Traj.x[,s], Traj.y[,s], lwd=3, type='l', lty='solid', add=TRUE, col=rainbow(10))
    points(Traj.x[1,], Traj.y[1,], col='black', pch=19, cex=1)

    dev.off()
    matplot(Traj.x[,1] - floor(Traj.x[1,1]), Traj.y[,s], lwd=3, type='l', lty='solid', add=TRUE, col=rainbow(10))
    matplot(sweep(Traj.x[,s],2,floor(Traj.x[1,s]),'-'), Traj.y[,s], lwd=3, type='l', lty='solid', add=TRUE, col=rainbow(10))
    points(Traj.x[1,s], Traj.y[1,s], col='green', pch=19, cex=1.5)
}
#' --------------------------------------------------------------------------


LFP <- data.matrix(fread("p0-a3-r12--LFP.dat"))
TH  <- data.matrix(fread("p0-a3-r12--TH.dat"))
plot(TH[,2] * 180/pi, LFP[,2])








#' --------------------------- DISK CHANNEL ---------------------------------
{
    theta <- seq(0, 2*pi, length=200)
    alpha <- pi/3
    r <- 0
    d <- 1
    
    id <- 1
    Traj.x <- ReadDataList(lf(paste0("--DiscTraj_")), colidx=1)
    Traj.y <- ReadDataList(lf(paste0("--DiscTraj_")), colidx=2)
    
    
    
    Traj.x <- ReadDataList(lf(paste0("--ContTraj_")), colidx=1)
    Traj.y <- ReadDataList(lf(paste0("--ContTraj_")), colidx=2)
    
    from <- -3; to <- 3
    
    for (i in seq(nrow(Traj.x))) {
        png(paste0("FRAMES/", i, "-frame.png"), width=1440, height=1024)
        par(mar=c(0,0,0,0), bg='lightgray')
        
        PlotDiskChannel_Test(r, alpha, d=d, from=from, to=to)
        for (id in seq(ncol(Traj.x))) {
            polygon(r * cos(theta) + Traj.x[i,id], r * sin(theta) + Traj.y[i,id], col=adjustcolor(col="blue",alpha.f=0.5))
            points(Traj.x[i,id], Traj.y[i,id], pch=19, col='black')
        }
        dev.off()
    }
    
    
    
    id <- 1
    Traj.x <- ReadDataList(lf(paste0("--DiscTraj")), colidx=1)
    Traj.y <- ReadDataList(lf(paste0("--DiscTraj")), colidx=2)
    PlotDiskChannel_Test(r, alpha, d=d, from=0, to=1)
    lines(Traj.x[,id], Traj.y[,id], col='black')
    
    lines(Traj.x[20:23,id], Traj.y[20:23,id], col='black')
}
#' --------------------------------------------------------------------------


setwd("~/Dropbox/ABS/Studies/S-Scatter/Data/")

Traj.x <- data.matrix(fread(lf("--x")))
Traj.y <- data.matrix(fread(lf("--y")))

V <- getChannelVertices(dx = 1/2, alpha = alpha, d = 1/2)
from = 0
to=20
PolygonBaseplate(V[,1], V[,2], from, to)
PlotChannelBoundary(V[,1], V[,2], from, to)

n <- 100
lines(Traj.x[,]+1/2, Traj.y[,], lwd=0.75, col='red')

bid <- data.matrix(fread(lf("--BID")))
text(x=0.9+seq(length(bid[1,]))-2, y=1.2, bid[1,], col='blue')
text(x=1.1+seq(length(bid[1,]))-2, y=1.2, paste0(": ", bid[2,]), col='red')

text(x=1, y=1.4, "(h,theta) = (0.2, +/- 0.0713)", col='black')

paste(as.vector(bid[1,-1]), collapse=",")
paste(as.vector(bid[2,-1]), collapse=",")


# bids = 64664666366466636646
dev.copy2pdf(file="negtrajonstrip.pdf")

