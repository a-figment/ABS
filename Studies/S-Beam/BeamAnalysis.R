setwd("Data/")

#' (@x, @y) Position of particle initially on a beam
ShowBeamEvolution <- function(x, y, dx, dy, d) {
    xv <- c(-dx, 0, dx, dx, 0, -dx)
    yv <- c(0,  dy, 0,  d, d+dy, d)
    PlotChannelBoundary(xVert=xv, yVert=yv, from=-1, to=1)
    
    for (i in seq(nrow(x))) {
	png(paste0(i, "-frame.png"), width=1440, height=1024)
	PlotChannelBoundary(xVert=xv, yVert=yv, from=-1, to=1)
	points(x[i,], y[i,], pch=19, col='red', cex=0.5)
	dev.off()
    }
}

#' @brief Plots initial distribution of beams for parallel channel
ShowBeams.PC <- function(dx, dy, d, beams=NA, n=0.05) {
    x <- c(-dx, 0, dx, dx, 0, -dx)
    y <- c(0,  dy, 0,  d, d+dy, d)
    PlotChannelBoundary(xVert=x, yVert=y, from=0, to=0)
    for (i in seq(nrow(beams))) {
	lines(beams[i,c(1,3)], beams[i,c(2,4)], lwd=2, col='blue')
	mp <- (beams[i,1:2] + beams[i,3:4]) / 2
	va <- beams[i,5]
	arrows(x0=mp[1], y0=mp[2], 
	       x1=mp[1]+n*cos(va), y1=mp[2]+n*sin(va), 
	       angle=45, length=0.1)
    }
}


### 
##  PARALLEL CHANNELS -- 
### 
{
    dx <- 1/2
    varphi <- (sqrt(5) + 1)/2
    alpha <- c(pi/seq(9, 2, by=-1), seq(2,8)*pi/seq(3,9))
    alpha <- c(alpha, alpha * (varphi - 3/5))
    ratio <- seq(1/4, 10, length=40)
}

KeyID <- c(7,2,0)	# Input id 
{
    a <- alpha[KeyID[1]+1] 
    r <- ratio[KeyID[2]+1]
    dy <- 1/(2*tan(a/2))
    d <- dy * r
    PlotChannelBoundary(xVert=c(-dx, 0, dx, dx, 0, -dx), yVert=c(0,  dy, 0,  d, d+dy, d), from=0, to=0)
}

beam <- data.matrix(fread("Beams.dat"))
ShowBeams.PC(dx, dy, d, beam[1:20,])

x <- ReadDataList(lf("Traj"), idx=1)
y <- ReadDataList(lf("Traj"), idx=2)
ShowBeamEvolution(x, y, dx, dy, d)

id  <- tpose(data.matrix(fread(lf("--ID"))))
cid <- tpose(data.matrix(fread(lf("--CID"))))

np <- x <- list()
for (i in seq(nrow(cid))) {
    r <- rle(cid[i,])
    np[[i]] <- r$lengths
    x[[i]] <- r$values
}

mplot(x[[1]], np[[1]], type='p', cex=1.5, col='blue', log='',
      xlim=c(-10, 10), ylim=c(1,2000))
      xlim=c(min(unlist(x)), max(unlist(x))),
      ylim=c(0, ncol(id)))

for (i in seq(2, length(x)))
    mplot(x[[i]], np[[i]], type='p', cex=.7, add=TRUE)



tau <- data.matrix(fread("FJ.dat"))[,1]
xi <- data.matrix(fread("FJ.dat"))[,2]
ra <- data.matrix(fread("FJ.dat"))[,3]

L <- getTSD_hist(dx = tau, n=1000)
mplot(L$x, L$y, log='')


# Discretised position as a function of initial position
th  <- fread(lff("p0", "theta.dat"))
x   <- tpose(fread(lff("p0", "x.dat")))
y   <- tpose(fread(lff("p0", "y.dat")))
tau <- tpose(fread(lff("p0", "tt.dat")))

PlotChannelBoundary(xVert=c(-dx, 0, dx, dx, 0, -dx), yVert=c(0,  dy, 0,  d, d+dy, d), from=-2, to=2)
j <- 7
rgl::plot3d(x=x[,j], y=y[,j], z=tau[,j], cex=0.2, type='p', xlab="x(0)", ylab="y(0)", zlab="[x]")
id <- which(tau[,j] == -1)
mplot(x[id,j], y[id,j], type='p', cex=0.1, log='')

mplot(x[id,j] + cos(4.937265)*y[id,j], x[id,j] + 0.375 - y[id,j], type='p', cex=0.1, log='')

mplot(sort(x[id,j] + cos(4.937265)*y[id,j]), y[id,j], type='p', cex=0.1)



L <- getTSD_ecdf(dx = tau[,], n = 1000)
L <- getTSD_hist(dx = tau[,10], n = 200)

L$x <- sweep(L$x, 2, apply(L$x, 2, function(x) min(x)), '-')
mplot(L$x, L$y, type='s')#, xlim=c(0,5))
summary(L$x)

ntransfers <- seq(1, 1000, by=100)
mt <- apply(tau, 2, mean)
mt.sdv <- apply(tau, 2, function(x) sqrt(var(x)))
mplot(ntransfers, mt, type='o')
mplot(ntransfers, mt, type='o', add=TRUE)
ErrorBars.y(ntransfers, mt, mt.sdv)

mplot(x[,], tau[,], cex=0.2, type='p')




frac <- function(x) x - floor(abs(x)) * (x/abs(x))
###
#   Plot evolution of slicer
###
SliceChannel <- function(x, y, xv, yv) {
    th <- 0.97025
    mplot(x[,1], y[,1], type='p', col='black', cex=.01, xlab="Distance", ylab="Time", xlim=c(-4, 4), ylim=c(0, 8))
    abline(h=seq(0, 1000, by=max(yv)), v=seq(-1000,1000)+1/2, lty='dashed')
    curve(expr = max(yv)*x - max(yv)/2, from=1/2, to=100, add=TRUE)
    curve(expr = -max(yv)*x - max(yv)/2, from=-1/2, to=-100, add=TRUE)

    for (i in seq(2,ncol(x))) {
        cols <- ifelse(i %% 2 == 0, rbow()[1], rbow()[2])
	
        #mplot(x[,1] + round(x[,i]), y[,1]+(i-1)*max(yv), add=TRUE, type='p', col=cols, cex=.01)

	X <- x[,1] + round(x[,i])
	Y <- y[,1] + (i-1)*max(yv)

	left <- which(round(x[,2]) < 0)
        mplot(X[left],  Y[left], add=TRUE, type='p', col='red', cex=.01)
        mplot(X[-left], Y[-left], add=TRUE, type='p', col='blue', cex=.01)

	left <- which(round(x[,3]) < -1)
	right <- which(round(x[,3]) > 1)
        mplot(X[left],  Y[left], add=TRUE, type='p', col='green', cex=.01)
        mplot(X[right], Y[right], add=TRUE, type='p', col='orange', cex=.01)

	left <- which(round(x[,4]) < -2)
	right <- which(round(x[,4]) > 2)
        mplot(X[left],  Y[left], add=TRUE, type='p', col='black', cex=.01)
        mplot(X[right], Y[right], add=TRUE, type='p', col='pink', cex=.01)

	left <- which(round(x[,5]) < -3)
	right <- which(round(x[,5]) > 3)
        mplot(X[left],  Y[left], add=TRUE, type='p', col='yellow', cex=.01)
        mplot(X[right], Y[right], add=TRUE, type='p', col='turquoise2', cex=.01)
        
	#w <- which(X < round(x[,i]))
	#mplot(X[w]+1  - y[w,1], Y[w],    add=TRUE, type='p', col='blue', cex=.01)
	#mplot(X[-w] + y[-w,1] - 1, Y[-w], add=TRUE, type='p', col='magenta', cex=.01)

        #mplot(sort(x[,i]), y[,i] + (i-1)*max(yv), add=TRUE, type='p', col=cols, cex=.01)
        #mplot(X - cos(th)*y[,1], Y, add=TRUE, type='p', col='blue', cex=.01)

        #mplot(x[,1] + round(x[,i]), y[,1]+(i-1)*max(yv), add=TRUE, type='p', col=cols, cex=.01)
    }
}

# @brief Counts the number of particles to the right of sweeping lins
#	 defined by the end points (@A, @B)
# @x0, @y0 [vector] Points which are being counted
# @A, @B [nx2 matrix] x-y coordinates of start/end points of sweep line
count <- function(x0, y0, A, B) {
    counts <- c()
    for (i in seq(nrow(A))) {
	l <- isRight(A[i,], B[i,], cbind(x0, y0))
	wl <- which(l == TRUE)
	counts <- c(counts, length(wl))
    }
    return(counts)
}

# Edge detection
getEdge <- function(s, counts, thresh=2) {
    r <- rle(counts)
    horiz <- which(r$lengths > 2) # Repeated elements > 5
    id <- c(1, sapply(r$values[horiz], function(x) { w <- which(counts == x); c(w[1], w[length(w)]) }))
    return(id)
}


# Define elementary cell here
xv <- c(-dx, 0, dx, dx, 0, -dx)
yv <- c(0,  dy, 0,  d, d+dy, d)

#pdf(paste0("../Images/pi3--SlicedChannels.pdf"), height=12, width=18)
for (i in (c(0))) {
    x   <- tpose(fread(lff(paste0("p", i, "-a7"), "x.dat")))
    y   <- tpose(fread(lff(paste0("p", i, "-a7"), "y.dat")))
    SliceChannel(x[,], y[,], xv, yv)
}

# Setup
{
    th <- 2.54128
    i <- 8
    X <- x[,1] #+ round(x[,i])
    Y <- y[,1] #+ (i-1)*max(yv)
    w <- which(round(x[,i]) == -1)

    mplot(X[w], Y[w], type='p', cex=.01,asp=1, )#xlim=c(-3, 3), ylim=c(-2,2))
    polygon(xv,yv, border='blue')
    abline(h=seq(0, 1000, by=max(yv)), v=seq(-1000,1000)+1/2, lty='dashed')
}

# Sweep line
N <- 200 # Number of sweep lines
M <- 10  # Number of points to define the sweep line
s <- seq(tan(th)/2, max(yv) - tan(th)/2, length=N) # y-intercepts for sweep line
x.sweep <- seq(-2, 2, length=M)  # Domain of sweeping line
y.sweep <- sapply(1:length(s), function(x) { x.sweep*tan(th) + s[x] }) # The sweeping line
#mplot(x.sweep, y.sweep, add=TRUE)  # Overlay all sweeping lines

# Get the start and end point of each sweep line
A <- cbind(x.sweep[1], y.sweep[1,])
B <- cbind(x.sweep[M], y.sweep[M,])
#mplot(cbind(A[,1], B[,1]), cbind(A[,2], B[,2]), type='p', add=TRUE)


# Considering counts based on LHS and RHS particles
xx <- X[w]
yy <- Y[w]
counter <- count(xx, yy, A, B)			   # Ordinary counter
lhs.count <- count(xx[xx <= 0], yy[xx <= 0], A, B) # Count only LHS points
rhs.count <- count(xx[xx > 0], yy[xx > 0], A, B)   # Count only RHS points

# Side-by-side plot of sliced instance and sliced locations
#pdf("../Images/SlicingMechanismExamplePROBLEMprint.pdf", width=12, height=6)
{
    spar(mfrow=c(1,2), cx=1)
    mplot(X[w], Y[w], type='p', cex=.01,asp=1, spa=FALSE)
    k <- 94
    mplot(x.sweep, y.sweep[,k], add=TRUE, spa=FALSE, col='black')
    ymid <- (y.sweep[1,k]+y.sweep[nrow(y.sweep),k])/2 
    arrows(x0=0, x1=-cos(th-pi/2)/10, y0=ymid, y1=ymid-sin(th-pi/2)/10, lwd=2)
    text(x=0, y=0.1, labels=paste0("count = ", counter[k]), cex=1)
    text(x=-.425, y=ymid+0.35, labels=TeX("$y = x \\tan(\\theta_0) + s_i$"), cex = 1)
    polygon(xv,yv, border='blue')
    abline(h=seq(0, 1000, by=max(yv)), v=seq(-1000,1000)+1/2, lty='dashed')
    
    mplot(s, counter, xlab="s", ylab="Counts", spa=FALSE, )#xlim=c(5,6), ylim=c(0,1e4))
    abline(v=s[k], lty='dashed')
    lines(s, lhs.count, lwd=2, col=rbow()[2])
    lines(s, rhs.count, lwd=2, col=rbow()[3])
    legend('topleft', legend=c("LHS + RHS", "LHS", "RHS", tex("s_i")), 
	   lty=c(rep('solid',3),'dashed'), col=c(rbow()[1:3], 'black'), cex=1, lwd=3)
    
    id     <- getEdge(s, counter)
    id.lhs <- getEdge(s, lhs.count)
    id.rhs <- getEdge(s, rhs.count)
    points(s[id], counter[id], pch=19, col=rbow()[1], cex=1)
    points(s[id.lhs], lhs.count[id.lhs], pch=19, col=rbow()[2], cex=1)
    points(s[id.rhs], rhs.count[id.rhs], pch=19, col=rbow()[3], cex=1)
}
dev.off()


# Visualise the sweeping line
counter <- count(xx, yy, A, B)
for (i in seq(nrow(A))) {
    mplot(X[w], Y[w], type='p', cex=.01,asp=1)
    polygon(xv,yv, border='blue')
    abline(h=seq(0, 1000, by=max(yv)), v=seq(-1000,1000)+1/2, lty='dashed')
    mplot(x.sweep, y.sweep[,i], add=TRUE)
    points(X[w][wl], Y[w][wl], col='blue', cex=0.1, pch=19)
    text(x=-0.65, y=0.1, labels=paste0("count = ", counter[i]), cex=2)
    Sys.sleep(.1)
}


# Trying to apply shearing operations to 
j <- 2
X <- x[,1] + round(x[,j])
w <- which(X < round(x[,j]))
mplot(X, y[,1], type='p', col=cols, cex=.01, ylim=c(-2,2), xlim=c(-5, 5))

points(X[w], y[w,1], col='blue', pch=19, cex=0.01)
points(X[-w], y[-w,1], col='red', pch=19, cex=0.01)

points(X[w] + y[w,1], y[w,1] - X[w], col='blue', pch=19, cex=0.01)
points(X[-w] - y[-w,1], y[-w,1] + X[-w], col='red', pch=19, cex=0.01)

points(X[w] - y[w,1], y[w,1] + X[w], col='blue', pch=19, cex=0.01)
points(X[-w] + y[-w,1], y[-w,1] - X[-w], col='red', pch=19, cex=0.01)


# Testing points left/right 
{
    #pdf(paste0("../Images/SweepLineCounter.pdf"), width=18, height=9)
    x <- cumsum(runif(20))
    X <- Y <- list()
    j <- 1
    for (i in seq(1,length(x),by=3)) {
	X[[j]] <- runif(n=10000 * (x[i+1]-x[i]), min=x[i], x[i+1])
	Y[[j]] <- runif(10000 * (x[i+1]-x[i]))
	j <- j + 1
    }
    X <- as.vector(list.cbind.na(X))
    Y <- as.vector(list.cbind.na(Y))

    spar(mfrow=c(1,2))
    mplot(X, Y, type='p', spa=FALSE)

    s <- seq(min(x), max(x), length=1000)
    nleft <- nright <- c()
    for (i in seq(length(s))) {
	nleft[i] <- length(which(X < s[i])) #/ length(X)
	nright[i] <- length(which(X > s[i])) #/ length(X)
    }
    mplot(s, nleft, xlab="s", ylab="Num. points < s", spa=FALSE)
    r <- rle(nleft)
    horiz <- which(r$lengths > 2) # Repeated elements > 5
    id <- c(1, sapply(r$values[horiz], function(x) { w <- which(nleft == x); c(w[1], w[length(w)]) }))
    points(s[id], nleft[id], pch=19, col='blue', cex=2)
    #abline(v=s[id])
    #abline(h=nleft[id])
    #dev.off()
}

###
#  
###

T <- 101 # Time runs from [1, T-1]
t <- seq(-(T-1), (T-1), by=1)
counts <- data.matrix(fread(lf('counts')))

mplot(t/(T-1), counts, type='s', log='')

# De-symmetrify - 
posid <- seq(T+1, length(t))
negid <- seq(1, T-1)
zerid <- T
counts <- data.matrix(c(counts[T,], counts[posid,] + rev(counts[negid,])))
t <- seq(0, T-1)
mplot(t, counts, type='s', )#ylim=c(0,1/2), log='')
mplot(t+1, counts, type='s', col='blue', add=TRUE)


###
#   Discretised moments
###

setwd("Data/T0/")

p <- seq(2, 20, by=2)
t <- 1:10001
mom <- ReadData(lff("", "Moments"))[,,1]
mom <- data.matrix(fread(lff("a6", "Moments"), colClasses=c("numeric"), skip=9))
mplot(t, mom[1,], type='l', log='xy')

plf <- PowerLawFits(x=t, (mom), fitrange=1000:10001)

mplot(p, plf$scales, type='o', xaxt='n', xlab=tex("p"), ylab=tex("\\gamma_p"), xlim=c(0,max(p)), ylim=c(0,max(p)), xaxs='i', yaxs='i')
axis(side=1, at=c(0,p))

# Rough
LinFit(mat=plf$scales, time=p, range = 1:2, plot=TRUE)

# Fit through the origin
df <- data.frame(x=p, y=plf$scales[,2])

l <- lm(y ~ 0 + ., data=df, subset=1:2)
curve(x*coef(l)[1], lty='dashed', from=0, to=max(p), add=TRUE)
l <- lm(y ~ x, data=df, subset=2:10)
curve(x*coef(l)[2] + coef(l)[1], lty='dashed', from=0, to=max(p), add=TRUE)
