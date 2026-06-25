
setwd("Data/")
{
    x <- data.matrix(fread(lf("x.dat")))
    #t <- data.matrix(fread(lf("t.dat")))
    y <- data.matrix(fread(lf("y.dat")))

    bp <- data.matrix(fread(lf("Boundary")))
    dy <- dx/tan(pi/4)
    d <- 2/10

    for (s in seq(555,562)) {
	# Plot channel as normal
	PlotChannelBoundary(xVert=c(-dx, 0, dx, dx, 0, -dx), yVert=c(0, dy, 0, d, dy+d,d), from=-1,to=1)

	# Add circles
	w <- 0.05
	r <- 0.035355
	DrawCircle(c(-dx, w),    r)
	DrawCircle(c(0, dy-w),   r)
	DrawCircle(c(dx, w),     r)
	DrawCircle(c(dx, d+w),   r)
	DrawCircle(c(0, dy+d-w), r)
	DrawCircle(c(-dx, d+w),  r)
    	points(bp[,1], bp[,2], pch=19, cex=0.5)
    	points(bp[,1]+1, bp[,2], pch=19, cex=0.5)

	#s <- rindex(nrow(x),1)[1]
	#s <- 1470
	n <- 1:3
	lines(x[s,], y[s,])
	points(x[s,n], y[s,n], pch=19, cex=0.5, col='red')
	points(x[,], y[,], pch=19, cex=0.5, col='red')

	Sys.sleep(0.5)
    }

    # Starting on periodic boundary whilst also being inside a scatterer
    # 
    s <- 1
    lines(x-X[s,1:2], y[s,1:2], pch=19, cex=0.5, col='black')
    points(x[s,1:5], y[s,1:5], pch=19, cex=0.5, col='red')
    arrows(x0=x[s,3], x1=dx, y0=y[s,3], y1=w, len=0.1) 

    ang <- atan(diff(y[s,])/diff(x[s,]))
    X <- x[s,3]
    ang <- atan(-0.8722/0.489)
    curve(expr=(x-X)*ang[1] + Y, from=-10, to=10, add=TRUE)
    y0 <- ((1/2)-X)*ang[3] + y[s,3]
    points(X, Y)
    curve(expr=(x+1/2)*ang[3] + y0, from=-10, to=10, add=TRUE)

    X <- 0.481214
    Y <- 0.020048
    lines(c(x[s,2],X), c(y[s,2], Y))

    X1 <- 0.525
    Y1 <- 0.03548
    lines(c(X, X1), c(Y, Y1))
    arrows(x0=X, x1=dx, y0=Y, y1=w, len=0.1) 
}


DrawRoundedCornerChannel <- function(bp, radius, from=-1, to=1) {

    d <- bp[6,2] - bp[5,2]
    a0 <- atan((bp[3,2] - bp[14,2]) / (bp[3,1] - bp[14,1]))
    a1 <- pi/2-atan((bp[2,2] - bp[14,2]) / (bp[2,1] - bp[14,1]))
    A <- seq(a1, a0, length=10000)

    #PolygonBaseplate(bp[,1], bp[,2], from=from-1, to=to+1)

    id1 <- which.max(bp[,2])
    id2 <- which.min(bp[,2])
    yv <- bp[,2]
    yv[id1] <- bp[id1,2] + radius
    yv[id2] <- bp[id2,2] - radius
    PolygonBaseplate(bp[-nrow(bp),1], yv, from=from-1, to=to+1)

    for (i in seq(from-1, to+1, by=1)) {

        polygon(c(bp[1:2,1], bp[9:10,1])+ 2*dx*i, c(bp[1:2,2], bp[9:10,2]), col='white', border="white")
        polygon(c(bp[3:4,1], bp[7:8,1]) + 2*dx*i, c(bp[3:4,2], bp[7:8,2]), col='white', border="white")

	lines(bp[1:2,1] + 2*dx*i, bp[1:2,2] , lwd=1.5, col='black')
	lines(bp[3:4,1] + 2*dx*i, bp[3:4,2] , lwd=1.5, col='black')
	lines(bp[7:8,1] + 2*dx*i, bp[7:8,2] , lwd=1.5, col='black')
	lines(bp[9:10,1]+ 2*dx*i, bp[9:10,2], lwd=1.5, col='black')

	xc <- bp[14,1] + radius*cos(A)
	yc <- bp[14,2] + radius*sin(A)
        polygon(c(bp[2,1], xc, bp[3,1], bp[8,1], xc, bp[9,1])+ 2*dx*i, 
		c(bp[2,2], yc, bp[3,2], bp[3,2]+d, yc+d, bp[3,2]+d), col='white', border="white")

	xc <- bp[15,1] - radius*cos(-A)
	yc <- bp[15,2] + radius*sin(-A)
        polygon(c(bp[4,1], xc, bp[5,1], bp[6,1], xc, bp[7,1])+ 2*dx*i,
		c(bp[4,2], yc, bp[5,2], bp[6,2], yc+d, bp[7,2]), col='white', border="white")

	A <- seq(a1, a0, length=10000)
	DrawCircle(center=c(bp[14,1] + 2*dx*i, bp[14,2]), radius, A, lwd=1.5)
	DrawCircle(center=c(bp[15,1] + 2*dx*i, bp[15,2]), radius,-A, lwd=1.5)
	DrawCircle(center=c(bp[16,1] + 2*dx*i, bp[16,2]), radius,-A, lwd=1.5)
	DrawCircle(center=c(bp[17,1] + 2*dx*i, bp[17,2]), radius, A, lwd=1.5)

    }

}

# Animation
{
    x <- data.matrix(fread(lf("x.dat")))
    y <- data.matrix(fread(lf("y.dat")))

    bp <- data.matrix(fread(lf("Boundary")))

    for (i in seq(ncol(x))) {
	png(paste0(i, "-frame.png"), width=2560, height=1440)#, res=1200)#, units='cm'Colorspace)
	par(mar=c(0,0,0,0), bg='gray')
	# Plot channel as normal
	DrawRoundedCornerChannel(bp, bp[nrow(bp),1], from=0,1)
	#PlotChannelBoundary(xVert=c(-dx, 0, dx, dx, 0, -dx), 
	#		    yVert=c(0, dy, 0, d, dy+d,d), from=-1,to=1)
	points(x[,i], y[,i], pch=19, cex=0.15)
	dev.off()
    }

}

    x <- tpose(fread(lf("x.dat")))
    t <- tpose(fread(lf("t.dat")))
    
    mplot(t, x, lwd=1, type='l')#, xlim=c(139800,139890), ylim=c(120,150))

#pdf("Ex2e.pdf", width=30, height=5)
{
    x <- data.matrix(fread(lf("x.dat")))
    y <- data.matrix(fread(lf("y.dat")))
    t <- data.matrix(fread(lf("t.dat")))
    
    mplot(t, x, lwd=0.5, type='l')
    bp <- data.matrix(fread(lf("Boundary")))

    DrawRoundedCornerChannel(bp, bp[nrow(bp),1], from=0,to=0)
    {
	s <- rindex(nrow(x), 1)[1]
	s <- 28
    	lines(x[s,], y[s,])
    	points(x[s,1], y[s,1], col='green', pch=19)
    	points(x[nrow(x),1], y[nrow(x),1], col='green', pch=19, cex=0.1)
    }
}
dev.off()


### DISPERSION FIGURE OF THESIS ###

# NON-PARALLEL CHANNEL
dx <- 1/2
n <- 10

{
    x <- data.matrix(fread(lf("x.dat")))
    y <- data.matrix(fread(lf("y.dat")))
    xVert <- c(-dx, 0, dx, dx, 0, -dx)
    #yVert <- c(0, 0.45, 0, 0.3, 1.07,0.3)
    yVert <- c(0, 0.45, 0, 0.2, 0.97,0.2)
    PlotChannelBoundary(xVert=xVert, yVert=yVert, from=-1,to=1)

    ts <- c(1, 12, 23, 38, 49, 60)
    j <- 1
    for (i in ts) {
	points(x[,i],y[,i], pch=19, cex=0.2, col=rbow()[j])
	#Sys.sleep(1)
	j <- j + 1
    }

    for (i in seq(floor(from), floor(to))) {
	lines(xVert[1:3] + i*1, yVert[1:3], col='black', lwd=1.5)
	lines(xVert[4:6] + i*1, yVert[4:6], col='black', lwd=1.5)
    }
}

# ROUNDED CHANNEL
{
    # Time recorded: seq(0, 2, length=201). Input params: 0 7 2 0
    x <- data.matrix(fread(lf("x1.dat")))
    y <- data.matrix(fread(lf("y1.dat")))
    bp <- data.matrix(fread(lf("Boundary")))
    DrawRoundedCornerChannel(bp, bp[nrow(bp),1], from=0,to=0)
    ts <- c(1, 12, 23, 38, 49, 60)
    j <- 1
    for (i in ts) {
	points(x[,i],y[,i], pch=19, cex=0.2, col=rbow()[j])
	#Sys.sleep(1)
	j <- j + 1
    }

    from <- -1
    to <- 1
    for (i in seq(from-1, to+1, by=1)) {
	lines(bp[1:2,1] + 2*dx*i, bp[1:2,2] , lwd=1.5, col='black')
	lines(bp[3:4,1] + 2*dx*i, bp[3:4,2] , lwd=1.5, col='black')
	lines(bp[7:8,1] + 2*dx*i, bp[7:8,2] , lwd=1.5, col='black')
	lines(bp[9:10,1]+ 2*dx*i, bp[9:10,2], lwd=1.5, col='black')
	A <- seq(a1, a0, length=5000)
	DrawCircle(center=c(bp[14,1] + 2*dx*i, bp[14,2]), radius, A, lwd=1.5)
	DrawCircle(center=c(bp[15,1] + 2*dx*i, bp[15,2]), radius,-A, lwd=1.5)
	DrawCircle(center=c(bp[16,1] + 2*dx*i, bp[16,2]), radius,-A, lwd=1.5)
	DrawCircle(center=c(bp[17,1] + 2*dx*i, bp[17,2]), radius, A, lwd=1.5)
    }
}

# Statistics of displacement
dx <- data.matrix(fread(lf("DX.dat")))
L <- getTSD_hist(dx=dx, n=500, t=c(100, 1000 ,1e4, 1e5))
mplot(L$x, L$y, type='s', log='y', xlim=c(-1, 1))#, xlim=c(-5000,5000))

L <- getTSD_ecdf(dx=abs(dx), n=1000, domain=exp.domain)
mplot(L$x, 1-L$y, type='s', log='xy', xlim=c(1e1, 1e4))

mom <- tpose(fread(lf("Moments")))
t <- seq(0, 1e5, length=1001)
mu <- getPowerScaling(t=t, M=mom[,2], z=c(0.7, 0.8, 0.9, 0.95))
mplot(t, mu[,], type='l')
mplot(t, mom[,2], log='xy')
LinFit(mat=mom[,2], time=t, plot = TRUE, range=200:1001)
abline(h=1)



# ROUNDED CHANNEL TRAJECTORIES
{
    dx <- 1/2
    bp <- data.matrix(fread(lf("Boundary")))
    DrawRoundedCornerChannel(bp, bp[nrow(bp),1], from=0,to=0)
    {
	x <- data.matrix(fread(lf("x1.dat")))
	y <- data.matrix(fread(lf("y1.dat")))
	lines(x,y, col='red', lwd=2)
	x <- data.matrix(fread(lf("x2.dat")))
	y <- data.matrix(fread(lf("y2.dat")))
	lines(x,y, col='blue', lwd=2)
	points(x[1], y[1], col='green', pch=19)
    }


    points(bp[-(1:12),1], bp[-(1:12),2])
    points(bp[13,1], bp[13,2])

    c1 <- cbind(c(bp[1:2,1]

    
    BP.skeleton <- bp[c(1:2,14,3:4,7:8,17, 9:10),]
    BP.centers  <- bp[13:18,]
    points(BP.skeleton[1,1], BP.skeleton[1,2])
    points(BP.skeleton[2,1], BP.skeleton[2,2])

    mplot(BP.skeleton[7:9,1], BP.skeleton[7:9,2], add=TRUE, pch=19, type='p')

    fwrite(x=cbind(BP.skeleton[,1] , BP.skeleton[,2]), file="BP.tkz", sep2 = c("(", ",", ")"), col.names=FALSE)

    points(BP.skeleton[1:4,1], BP.skeleton[1:4,2])
    points(BP.skeleton[4,1]-1, BP.skeleton[4,2], col='red')

    points(bp[-nrow(bp),1], bp[-nrow(bp),2], pch=19)
    abline(h=0.2)
    lines(x,y, lwd=2)
    x <- data.matrix(fread(lf("x1.dat")))
    y <- data.matrix(fread(lf("y1.dat")))
}


pdf("DisperseEx3--sharp.pdf", width=18, height=10)
{
    Vx <- c(-1/2, 0, 1/2, 1/2, 0, -1/2)

    Vy <- c(0, dyb, 0, d, dyt+d, d)
    Vy <- c(0, dyb, 0, .3, 1.07, 0.3)
    PlotChannelBoundary(xVert=Vx, yVert=Vy, from=-1, to=1, mult=1)
    x <- data.matrix(fread(lf("x2.dat")))
    y <- data.matrix(fread(lf("y2.dat")))
    lines(x, y, col='red')
    points(x[1], y[1], pch=19, col='green')

    x <- data.matrix(fread(lf("x1.dat")))
    y <- data.matrix(fread(lf("y1.dat")))
    lines(x, y, col='blue')
    points(x[1], y[1], pch=19, col='green')


    #x1 <- tpose(fread(lf("x1.dat")))
    #y1 <- tpose(fread(lf("y1.dat")))
    #x2 <- tpose(fread(lf("x2.dat")))
    #y2 <- tpose(fread(lf("y2.dat")))
   # fwrite(x=cbind(x1+0.5, y1, x2+0.5, y2), file="NPTraj.tkz", col.names=FALSE, sep=" ")

}
dev.off()

pdf("DisperseEx3--rounded.pdf", width=18, height=10)
dx <- 1/2
bp <- data.matrix(fread(lf("Boundary")))
DrawRoundedCornerChannel(bp, bp[nrow(bp),1], from=0,to=0)
{
    x <- data.matrix(fread(lf("x1.dat")))
    y <- data.matrix(fread(lf("y1.dat")))
    lines(x,y, col='red', lwd=2)
    x <- data.matrix(fread(lf("x2.dat")))
    y <- data.matrix(fread(lf("y2.dat")))
    lines(x,y, col='blue', lwd=2)
    points(x[1], y[1], col='green', pch=19)
}
dev.off()


