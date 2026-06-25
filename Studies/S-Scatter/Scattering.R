

# Real trajectories
#pdf("TrajEx6-shorter.pdf", width=15, height=8)
PlotTraj <- function(d, from=0,to=1) {
    x <- data.matrix(fread(lf("x.dat")))
    #t <- data.matrix(fread(lf("t.dat")))
    y <- data.matrix(fread(lf("y.dat")))

    dx <- dy <- 1/2
    dy <- dx
    #d <- 2.414
    #d <- dy/2
    #PlotChannelBoundary(xVert=c(-dx, 0, dx, dx, 0, -dx), yVert=c(0, dy, 0, d, dy+d,d), from=-2,to=2)
    PlotChannelBoundary(xVert=c(0, dx, 2*dx, 2*dx, dx, 0), yVert=c(0, dy, 0, d, dy+d,d), from=from,to=to)
    lines(x, y)
    lines(x[1:2], y[1:2], col='blue')
    

    atan(diff(y)/diff(x))

    #J <- 30
    #lines(x[1:J], y[1:J])
    #points(x[1], y[1], col='green', pch=19)
    #points(x[J], y[J], col='red', pch=19)
}
dev.off()

# Compare billiard trajectory with scattering map

    #d <- 1/2
    #h0 <- y[1]
    #a1 <- -atan(h0/d) + 2*pi
    #a2 <- -atan(h0/(2*d)) + 2*pi
    #phi0 <- -atan((h0+d)/(4*d)) + 2*pi
    #phi1 <- -atan((h0+d)/(5*d)) + 2*pi
    #curve(tan(a1)*x + h0, from=-10, to=10, add=TRUE)
    #curve(tan(a2)*x + h0, from=-10, to=10, add=TRUE)
    #curve(tan(phi0)*x + h0, from=-10, to=10, col='green', add=TRUE)
    #curve(tan(phi1)*x + h0, from=-10, to=10, col='green', add=TRUE)

    #abline(v=1/2, h=0)
    
h0 <- ReadData(lf("Branch-1"), REV = TRUE)[,1,1]
t0 <- ReadData(lf("Branch-1"), REV = TRUE)[,2,1]
t0 <- (((t0 -pi) %% (2*pi)) - pi) * 180/pi
mplot(h0, t0*180/pi, type='p', cex=0.1, 
      xlim=c(0,1/2), ylim=c(0,360))


# Initial conditions
#png(paste0("ICSpace/ICSpace-", 0, ".png"), height=1080, width=1080)
#pdf("ICSpace.pdf", height=18, width=18)
spar(mfrow=c(5,4))
par(mar=c(0,0,0,0))
mplot(0,0, col='white', spa=FALSE,
      xlab=tex("h_0"), ylab=tex("\\theta_0"),
      xlim=c(0,1/2), ylim=c(0,360))
for (k in seq(0,18)) {
    ht <- tpose(fread(paste0("p", k, "-a7-r0--H0T0.dat")))
    h1 <- ht[,1]
    t1 <- ht[,2]
    #mplot(h1[], t1[]*180/pi, type='p', cex=0.25, col=cols[k+1], spa=FALSE, add=TRUE)
    mplot(h1[], t1[]*180/pi, type='p', cex=0.25, col=cols[k+1], xlim=c(0,1/2), ylim=c(0,360), spa=FALSE, xaxt='n', yaxt='n', xlab='', ylab='')
}
dev.off()

N <- 19
for (j in seq(N)) {

    png(paste0("ICs-", j, ".png"), height=1440, width=1440)
    #spar(mfrow=c(5,4))
    #par(mar=c(0,0,0,0))
#    mplot(0,0, col='white', spa=FALSE,
#	      xlab=tex(paste0("h_{", j, "}")), ylab=tex(paste0("\\theta_{", j, "}")),
#	      xlim=c(0,1/2), ylim=c(0,360))
    h1 <- ReadDataList(lf("-a7-r0--H.dat"), idx=j)
    t1 <- ReadDataList(lf("-a7-r0--A.dat"), idx=j)
    for (i in seq(0,18)) {
	i <- 0
	#mplot(h1[,i+1], t1[,i+1]*180/pi, type='p', cex=0.1, col=cols[i+1], spa=FALSE, add=TRUE)
	mplot(h1[,i+1], t1[,i+1]*180/pi, type='p', cex=0.1, col=cols[i+1], xlim=c(0,1/2), ylim=c(0,360), xaxt='n', yaxt='n', xlab='', ylab='')#, spa=FALSE)
	abline(h=c(90,180,270))
    }
    dev.off()
}


# Branch ID's

BID <- data.matrix(fread(lf("BID")))
BID <- data.matrix(fread(lf("--L.dat")))
BID <- apply(BID, 2, diff)
for (i in seq(ncol(BID)-1)) {
    png(paste0(i, "-RegionCounts.png"), height=1080, width=1440)
    hist(BID[,i], xlab=tex("ID"), ylab="Frequency", main="",breaks=100)# xaxt='n', cex=3, ylim=c(0,20000), right=FALSE)
    axis(side=1, at=0:20+0.5, labels=0:20)
    dev.off()
}

    h <- tpose(fread(lf("-a7-r0--H.dat")))
    a <- tpose(fread(lf("-a7-r0--A.dat")))
    x <- tpose(fread(lf("-a7-r0--R.dat")))


    mplot(h, a*180/pi, type='p',
	  xlim=c(0,1/2), ylim=c(0,360))

    mplot(x, a*180/pi, type='p',
	  ylim=c(0,360))

    mplot(x, h, type='p', cex=0.5,
	  ylim=c(0,1/2))

    x <- ReadDataList(lf("_R.dat"), idx=1)
    n <- round(seq(1, 1e8, length=(1e9/1000)))
    n <- seq(1, 1e8+1, length=(1e8/100)+1)
    mplot(x[1:1000,], type='l', cex=0.5)
	  
{
    from=0
    to=1
    x <- data.matrix(fread(lf("x.dat")))
    y <- data.matrix(fread(lf("y.dat")))

    d <- 1/2
    dx <- dy <- 1/2
    dy <- dx

    # Height/position from the billiard simulation trajectory 
    BilliardH <- tpose(fread(lf("RawH")))
    BilliardX <- tpose(fread(lf("RawX")))
    SMH <- tpose(fread(lf("--H.dat")))
    SMX <- tpose(fread(lf("--R.dat")))
    SMA <- tpose(fread(lf("--A.dat")))

    #PlotChannelBoundary(xVert=c(0, dx, 2*dx, 2*dx, dx, 0), yVert=c(0, dy, 0, d, dy+d,d), from=from,to=to)
    plot(c(0,1), c(0,1), col='white', asp=1, xaxt='n', yaxt='n', xlab='', ylab='')
    polygon(x=c(0, dx, 2*dx, 2*dx, dx, 0), y=c(0, dy, 0, d, dy+d,d))
    sub <- 1:(min(which(x[,1] > 1)) - 1)
    x <- c(x[sub], 1)
    y <- c(y[sub], SMH[1])
    lines(x, y)
    lines(x[1:2], y[1:2], col='blue')
    idx <- 1:min(nrow(SMH), nrow(BilliardH))
    points(BilliardX[idx,], BilliardH[idx,], cex=2, lwd=2)
    points(SMX[idx,], SMH[idx,], pch=19, col='red')
}

s <- seq(200)
mplot(t[s], x[s])
Scatter <- function(h0, theta0, d) {
    k0 <- ceiling((d - h0 - tan(theta0))/(d*(tan(theta0)-1)))
    h1 <- h0 + k0*d*(tan(theta0) - 1) + tan(theta0)
    theta1 <- theta0
    tau <- (1 + k0*d) * 1/cos(theta0)
    return(c(h1, theta1, tau))
}
time <- c(0)
h0 <- 0.23073
th0 <- 0.51756
for (i in seq(100)) {
    S <- Scatter(h0, th0, d)
    h0 <- S[1]
    th0 <- S[2]
    time <- c(time, S[3])
}
lines(cumsum(time), seq(length(time)))

# Ballistic speeds for h_0 = d

d <- seq(0.05, 1/2-0.001, length=10000)
ths <- function(p, d) { atan(1/(1 + 1/(p*d))) }

#anu <- vector(mode='list', length=length(d))
ath <- ata <- vector(mode='list', length=length(d))
anu <- c() #vector(length=length(d))
for (i in seq(length(d))) {

    p.min <- ceiling((1 - 2*d[i])/(2*d[i]^2))
    p.max <- floor(1/(2*d[i]^2))

    thetas <- ths(p.min:p.max, d[i])

    nu.p <- 1/(p.min:p.max * sqrt(thetas^2 + 1))
    #anu[[i]] <- nu.p
    anu <- c(anu, nu.p)
    ath[[i]] <- thetas # Periodic angles
    ata[[i]] <- 1/nu.p # Period times
}
ata <- list.cbind.na(ata)
ath <- list.cbind.na(ath)

#pdf("h0eqd.pdf", width=18, height=9)
spar(mfrow=c(1,2))
mplot(d, tpose(ata), type='p', log='xy',cex=0.2, xlab="d", ylab=tex("\\tau_p"))
mplot(d, (ata[1,]), type='p', log='',cex=0.2, xlab="d", ylab=tex("\\tau_p"))
LinFit(mat=ata[1,], time=d, plot = TRUE, range=1:200)
lines(d, 100/d)
mplot(d, tpose(ata), type='p', cex=0.2, spa=FALSE, xlab="d", ylab=tex("\\tau_p"))
mplot(d, tpose(ath), type='p', cex=0.2, spa=FALSE, xlab="d", ylab=tex("\\theta_0^*"))
dev.off()

L <- getTSD_hist(dx = 1/as.vector(anu), n=10000)
L <- getTSD_ecdf(dx = 1/as.vector(anu), n=1000)#, domain=exp.domain)

mplot(L$x, L$y, type='s', log='', ylim=c(0,0.1))
mplot(L$x, L$y, type='s', log='', ylim=c(0,0.02))

X <- runif(min=0, pi/4, n=10000)

L <- getTSD_ecdf(dx = 1/sqrt(X^2 + 1), n=1000)
lines(L$x, 1-L$y)

# Angles in this branch
h0 <- runif(min=0, 1/2, n=10000)
tau <- c()
for (i in seq(length(h0))) {
    a0 <- runif(atan(-h0[i]/2), atan((1/2-h0[i])/2), n=1)
    lt0 <- which(a0 < 0)
    j1 <- ceiling(2*(1/2 - h0[i])/tan(a0[-lt0]))
    j2 <- floor(-2*h0[i]/tan(a0[lt0]))
    j <- c(j1, j2)
    tau <- c(tau, 2 * j / cos(a0))
}

L <- getTSD_ecdf(dx=tau, n=1000, domain=exp.domain)
L <- getTSD_hist(dx=tau, n=1000)#, domain=exp.domain)
mplot(L$x, 1-L$y, type='s', log='xy')
LinFit(mat=1-L$y, time=L$x, plot = TRUE, range=1:300)


###
#   COMPLETE SCATTERING MAP: 
###
In <- function(val, bounds) { (val > bounds[1] && val < bounds[2]) }

# Integer functions used
kappa0 <- function(d, h0, t0) { ceiling((d - h0 - tan(t0))/(d*(tan(t0) -1))) }
gamma0 <- function(d, h0, t0) { ceiling((h0 - 1 - d)/(d*(1 - tan(t0)))) }
alpha0 <- function(d, h0, t0) { ceiling((-h0 - tan(t0))/(d*(1 + tan(t0)))) }


#' Updates the position of particles as they cross boundary
UpdatePosition <- function(t0, t1, x0, sym) {
    print(paste0("t0 = ", t0 * 180/pi, ", t1 = ", t1*180/pi))
    # cos(t0) > 0 always - by the time we enter this function we have already moved into the symmetry branch if the direction was negative

    if (sym && cos(t1) > 0) { #) || (sym && cos(t0) > 0)) { # Moving left
        print("moving back")
        x1 <- x0 - 1	
    } else if (!sym && cos(t1) > 0) {
        print("moving forward")
        x1 <- x0 + 1
    } else {
        print("null")
        x1 <- x0
    }
}

cot <- function(a) { cos(a)/sin(a) }

{
    #' @brief Scattering map for a first collision on Gamma0
    #' @h0, @t0, @x0 Entrance height, angle and current position
    a0  <- function(d,h0) { -pi/2 }
    a1  <- function(d,h0) { -atan(h0/d) }
    a2  <- function(d,h0) { -atan((h0+d)/2) }
    a3  <- function(d,h0) { -atan(h0) }
    a4  <- function(d,h0) { -atan(2*h0/3) }
    a5  <- function(d,h0) { -atan(h0/2) }
    a6  <- function(d,h0) { +atan((d-h0)/2) }
    a7  <- function(d,h0) { +atan(2*(d-h0)/3) }
    a8  <- function(d,h0) { +atan(d-h0) }
    a9  <- function(d,h0) { +atan((2*d-h0)/2) }
    a10 <- function(d,h0) { +atan(2*(2*d-h0)/3) }
    a11 <- function(d,h0) { +atan(2*(3*d-h0)/3) } # Note: angle doesn't coincide with the lower bound for a first bounce on Gamma2

    # Gamma 2 angles
    a12 <- function(d, h0) { atan(1 - 2*h0) }
    a13 <- function(d, h0) { atan(1 - 2*h0/(1+2*d)) }
    a14 <- function(d, h0) { atan(1 - h0/(d + 1)) }
    a15 <- function(d, h0) { atan(1 + d - h0) }
    a16 <- function(d, h0) { atan(1 + 2*(d - h0)) }

    # Gamma 3 angles
    a17 <- function(d, h0) { atan(1 + 2*(d - h0)) }
    a18 <- function(d, h0) { atan(2*d + 1 - h0) }
    a19 <- function(d, h0) { atan(4*d + 1 - 2*h0) }
    #a20 <- function(d, h0) { atan((2*(d - h0) + 1)/(1-2*d)) }
    a20 <- function(d, h0) { pi/2 }
}

Sg0 <- function(d,h0,t0,x0, sym) { # BOTTOM LEFT
    X <- c(NA, NA) # = (h0, t0)
    if (In(t0, c(a0(d,h0),  a1(d,h0)))) {
	X <- c(-h0*cot(t0), -t0 + pi/2) # Note: Incoming angle is always negative
    } else if (In(t0, c(a1(d,h0),  a2(d,h0)))) {
	#X <- c(h0 + (alpha0(d,h0,t0) + 1)*(tan(t0) + d), t0)
	X <- c(h0 + d*(alpha0(d,h0,t0) + 3)*tan(t0) + d*(alpha0(d,h0,t0) + 1), t0)
    } else if (In(t0, c(a2(d,h0),  a3(d,h0))) && (h0 > 1/3)) {
	X <- c((h0+d)*cot(t0) + 5/2, 3*pi/2 - t0)
    } else if (In(t0, c(a2(d,h0),  a3(d,h0))) && (h0 < 1/3)) {
	X <- c(h0 + d + 5*tan(t0)/2, t0)
    } else if (In(t0, c(a3(d,h0),  a4(d,h0)))) {
	X <- c(h0+(1+d)*tan(t0) + d, t0)
    } else if (In(t0, c(a4(d,h0),  a5(d,h0)))) { # Checked int
	X <- c(h0*cot(t0) + 4*d, -t0 + 3*pi/2)	 # t0 always negative
    } else if (In(t0, c(a5(d,h0),  a6(d,h0)))) { 
	X <- c(h0 + 4*d*tan(t0), t0)
    } else if (In(t0, c(a6(d,h0),  a7(d,h0)))) {
	X <- c(h0 + (2 + d)*tan(t0) - d, t0)
    } else if (In(t0, c(a7(d,h0),  a8(d,h0)))) {
	X <- c(-h0 + 2*(d - tan(t0)), t0 + pi)
    } else if (In(t0, c(a8(d,h0),  a9(d,h0)))) {
	X <- c(-h0 + 1 - 2*tan(t0), t0 + pi)
    } else if (In(t0, c(a9(d,h0),  a10(d,h0)))) {
	X <- c(4*d - (1-h0)/tan(t0), t0 + pi/2)
    } else if (In(t0, c(a10(d,h0), a11(d,h0)))) {
	X <- c(-h0 + 3*d - (1 + d)*tan(t0), t0 + pi)
    } else {
	print("Impossible angle in Sg0.")
    }
    X <- c(X, UpdatePosition(t0, X[2], x0, sym))
    return(X)
}

Sg2 <- function(d, h0, t0, x0, sym) { # TOP RIGHT
    X <- c(NA, NA)
    if (In(t0, c(a12(d,h0),  a13(d,h0)))) {
	X <- c(h0 + d*kappa0(d,h0,t0)*(tan(t0)-1) + tan(t0), t0)
    } else if (In(t0, c(a13(t0,h0),  a14(d,h0)))) {
	X <- c(-h0 + (1 - tan(t0))*(1 + d), t0 + pi)
    } else if (In(t0, c(a14(d,h0),  a15(d,h0)))) {
	X <- c(cot(t0)*(h0 + (d + 1)*(tan(t0) - 1)), -t0-pi/2)#-(t0 + pi/2))
    } else if (In(t0, c(a15(d,h0),  a16(d,h0)))) {
	X <- c(-h0 + 2*d + 1 - tan(t0), t0 + pi)
    } else {
	print("Impossible angle in Sg2.")
    }
    X <- c(X, UpdatePosition(t0, X[2], x0, sym))
    return (X)
}

Sg3 <- function(d, h0, t0, x0, sym) { # TOP LEFT
    X <- c(NA, NA)
    if (In(t0, c(a17(d,h0),  a18(d,h0)))) {
	X <- c(2*d + 1 - tan(t0) - h0, t0 + pi)
    } else if (In(t0, c(a18(d,h0),  a19(d,h0)))) {
	X <- c(cot(t0)*(h0 - 1 - 2*d + tan(t0)), 3*pi/2 - t0)
    } else if (In(t0, c(a19(d,h0),  a20(d,h0)))) {
	X <- c(cot(t0) * (-h0 + 1 + d*gamma0(d,h0,t0)*(1-tan(t0)) + d*(1 + tan(t0))) , -t0 + pi/2)
    } else {
	print("Impossible angle in Sg3.")
    }
    X <- c(X, UpdatePosition(t0, X[2], x0, sym))
    return(X)
}

THETA0 <- function(d, h0, t0) { c(-pi/2, atan(1 - 2*h0)) }
THETA2 <- function(d, h0, t0) { c(atan(1 - 2*h0), atan(1 + 2*d - 2*h0)) }
THETA3 <- function(d, h0, t0) { c(atan(1 + 2*d - 2*h0), pi/2) }

#' @brief Global scattering function
Scatter <- function(d, h0, t0, x0, sym=FALSE) {
    if (In(t0, THETA0(d,h0))) {	       # BOTTOM LEFT
	print("First bounce on BOTTOM LEFT")
	Sg0(d,h0,t0,x0, sym)
    } else if (In(t0, THETA2(d,h0))) { # TOP RIGHT
	print("First bounce on TOP RIGHT")
	Sg2(d,h0,t0,x0, sym)
    } else if (In(t0, THETA3(d,h0))) { # TOP LEFT
	print("First bounce on TOP LEFT")
	Sg3(d,h0,t0,x0, sym)
    } else {
	#print("an impossible angle!?")
	print("Symmetry branch.")
	S <- Scatter(d,h0,pi-t0,x0,sym=TRUE)
	S[2] <- pi - S[2] # Revert angle
	S

    }
}



# Positions are correct but something is wrong with symmetry branch
{
    d  <- 1/2
    PlotTraj(d, from=1,to=5)
    ha <- data.matrix(fread(lf("Height")))
    h0 <- ha[1,1] 
    t0 <- ha[1,2]
    if (t0 > 3*pi/2 && t0 < 2*pi) { t0 <- t0 %% (-pi/2) }
    S <- Scatter(d, h0, t0, x0=0)
    points(S[3], S[1], col='blue', pch=19)

    h <- x <- 0
    for (i in seq(20)) {
	h0 <- S[1]
	t0 <- S[2]
	if (t0 > 3*pi/2 && t0 < 2*pi) { t0 <- t0 %% (-pi/2) }
	x0 <- S[3]
	S <- Scatter(d, h0, t0, x0)
	points(S[3], S[1], col='red', pch=19)
	h <- c(h, h0)
	x <- c(x, x0)
    }
    #points(-x[-1], h[-1], col='blue', pch=19)
}

curve(tan(pi-a4(d,h0))*(x-S[3]) + h0, from = 0, to=10, add=TRUE, lty='dashed', col='red', lwd=1)
curve(tan(pi-a5(d,h0))*(x-S[3]) + h0, from = 0, to=10, add=TRUE, lty='dashed', col='red', lwd=1)


################################### TESTING Sg0 #########################################
#########################################################################################

ilist <- list(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11)
d <- dy <- 1/2
for (i in seq(length(ilist)-1)) { # Loop over angle domains
    h0 <- runif(min=0, max=d,n=1)
    t0 <- runif(min=ilist[[i]](d,h0), max=ilist[[i+1]](d,h0),n=1)
}

# Checking entrance ranges
{
curve(tan(a0(d,h0))*x + h0, from = 0, to=10, add=TRUE, lty='dashed', col='red', lwd=1)
curve(tan(a1(d,h0))*x + h0, from = 0, to=10, add=TRUE, lty='dashed', col='red', lwd=1)
curve(tan(a2(d,h0))*x + h0, from = 0, to=10, add=TRUE, lty='dashed', col='red', lwd=1)
curve(tan(a3(d,h0))*x + h0, from = 0, to=10, add=TRUE, lty='dashed', col='red', lwd=1)
curve(tan(a4(d,h0))*x + h0, from = 0, to=10, add=TRUE, lty='dashed', col='red', lwd=1)
curve(tan(a5(d,h0))*x + h0, from = 0, to=10, add=TRUE, lty='dashed', col='red', lwd=1)
curve(tan(a6(d,h0))*x + h0, from = 0, to=10, add=TRUE, lty='dashed', col='red', lwd=1)
curve(tan(a7(d,h0))*x + h0, from = 0, to=10, add=TRUE, lty='dashed', col='red', lwd=1)
curve(tan(a8(d,h0))*x + h0, from = 0, to=10, add=TRUE, lty='dashed', col='red', lwd=1)
curve(tan(a9(d,h0))*x + h0, from = 0, to=10, add=TRUE, lty='dashed', col='red', lwd=1)
}

curve(tan(pi-a0(d,h0))*x + h0, from = -10, to=0, add=TRUE, lty='dashed', col='red', lwd=1)
curve(tan(pi-a7(d,h0))*x + h0, from = -10, to=0, add=TRUE, lty='dashed', col='red', lwd=1)

