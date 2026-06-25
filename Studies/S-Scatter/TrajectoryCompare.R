# Compares scattering map trajectory with billiard trajectory

d <- 1

CompareTrajectory <- function(d) {
    dx <- dy <- 1/2
    PlotChannelBoundary(xVert=c(0, dx, 2*dx, 2*dx, dx, 0), yVert=c(0, dy, 0, d, dy+d,d), from=-2.5,to=2)
    i <- 0
    {
	# Billiard trajectory
	x <- data.matrix(fread(lff(paste0("p", i, "-"), "--x.dat")))
	y <- data.matrix(fread(lff(paste0("p", i, "-"), "--y.dat")))
	lines(x,y)

	# Heights from billiard trajectory
	h <- data.matrix(fread(lff(paste0("p", i, "-"), "--RawH.dat")))
	x <- data.matrix(fread(lff(paste0("p", i, "-"), "--RawX.dat")))
	points(x, h, cex=1.5, lwd=2)

	# Height from scattering map
	ST <- data.matrix(fread(lff(paste0("p", i, "-"), "--ST.dat")))
	points(ST[,3], ST[,1], pch=19, col='red')
    }
}

PlotTrajectory <- function(d, i=0, crossingPoint=FALSE) {
    dx <- dy <- 1/2
    PlotChannelBoundary(xVert=c(0, dx, 2*dx, 2*dx, dx, 0), yVert=c(0, dy, 0, d, dy+d,d), 
			from=-10,to=20)
    {
	# Billiard trajectory
	x <- data.matrix(fread(lff(paste0("p", i, "-"), "--x.dat")))
	y <- data.matrix(fread(lff(paste0("p", i, "-"), "--y.dat")))
	lines(x,y)

	if (crossingPoint) {
	    # Heights from billiard trajectory
	    h <- data.matrix(fread(lff(paste0("p", i, "-"), "--RawH.dat")))
	    x <- data.matrix(fread(lff(paste0("p", i, "-"), "--RawX.dat")))
	    #points(x, h, cex=1.5, lwd=2)

	    # Height from scattering map
	    ST <- data.matrix(fread(lff(paste0("p", i, "-"), "--ST.dat")))
	    points(ST[,3], ST[,1], pch=19, col='red')
	}
	points(x[1], y[1], col='green', pch=19, cex=2)
    }
}

PlotTrajectory(1, i=19, crossingPoint=TRUE)
# Persistance state
PS <- data.matrix(fread(lf("Persist")))
points(PS[,3], PS[,1], pch=19, cex=2, col='orange')

