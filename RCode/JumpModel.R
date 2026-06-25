###
#    Synthetic position evolution based on jump time definition 
###

rint <- function(n0, n1, n) floor(runif(n=n,min=n0,max=n1+1))
choose.speed <- function(m) { m[rint(1, n0=1, n1=length(m))] }
m <- c(0.6)#, 0.28) 
speeds <- c(m, -m)

#' @base.occ Number of jump time samples 
#' @base.repl Number of times to replicate @base.occ samples
#' @speeds Vector of jump speeds
getSyntheticPos <- function(base.occ=5, base.repl=2, speed=c(1/2, -1/2)) {

    #N <- 0
    #while (N < base.occ) {
	tau <- c(0, power.law(alpha=2, x0=0.1, n=base.occ))
	#tau <- tau[-which(tau > 1)]
	#N <- length(tau)
    #}
    tau <- rep(tau, base.repl)

    t <- cumsum(tau)
    xi <- c(); xi[1] <- 0
    for (i in seq(2, length(tau))) {
	#xi[i] <- xi[i-1] + choose.speed(speed) * tau[i]
	r <- runif(1)
	if (r < 1/2)
	    xi[i] <- xi[i-1] + speed[1] * tau[i]
	else 
	    xi[i] <- xi[i-1] + speed[2] * tau[i]
    }

    return(list(x=t, y=xi))
}

{
    P <- getSyntheticPos() 
    
    spar(cx=2.5, col='white')
    par(mar = c(8, 8, 1.3, 1.3), oma = c(2, 2, 2, 2), mgp=c(6,2,0))
    matplot(P$x, P$y, type='l', lwd=2, 
	    xlab=tex("t_k = \\sum_k \\tau_k"),
	    ylab=tex("x_k = x_{k-1} \\pm m \\tau_k"))
    points(P$x, P$y, pch=19, col='blue')

}

#' Interpolate and find jump times
t <- R$x
xi <- R$y
x <- seq(min(t), max(t), length=2*length(tau))
f <- approxfun(x = t, y = xi)(x)
    # If original data set included
    O <- order(c(x, t))
    x <- c(x, t)[O]
    f <- c(f, xi)[O]
points(x, f, pch=19, cex=1.5, col='blue')
points(t, xi, pch=19, cex=1.5, col='blue')

M <- getMM(x = x, y = f, w = 10)
M <- getMM(x = R$x, y = R$y, w = 1)
ShowJT(df = M$df)


#' @base.occ Number of jump time samples 
#' @base.repl Number of times to replicate @base.occ samples
#' @speeds Vector of jump speeds
getSyntheticPos <- function(base.occ=5, base.repl=2, speed=c(1/2, -1/2)) {

    #N <- 0
    #while (N < base.occ) {
	tau <- c(0, power.law(alpha=2, x0=0.1, n=base.occ))
	#tau <- tau[-which(tau > 1)]
	#N <- length(tau)
    #}
    tau <- rep(tau, base.repl)

    t <- cumsum(tau)
    xi <- c(); xi[1] <- 0
    for (i in seq(2, length(tau))) {
	#xi[i] <- xi[i-1] + choose.speed(speed) * tau[i]
	r <- runif(1)
	if (r < 1/2)
	    xi[i] <- xi[i-1] + speed[1] * tau[i]
	else 
	    xi[i] <- xi[i-1] + speed[2] * tau[i]
    }

    return(list(x=t, y=xi))
}

#' Broken
reflect <- function(x, y) {
    x.refl <- c(x, x[length(x)] + cumsum(rev(tau)))
    y.refl <- c(y, y[length(x)] + rev(y))

    rm <- length(x.refl)
    return(list(x=x.refl, y=y.refl))
}

flip <- function(x, y) {
    x.refl <- c(x, x[length(x)] + x)
    y.refl <- c(y, y[length(x)] - y)

    rm <- length(x.refl)
    return(list(x=x.refl, y=y.refl))
}


continue <- function(x, y) {
    x.cont <- c(x, x[length(x)] + x)
    y.cont <- c(y, y[length(x)] + y)
    return(list(x=x.cont, y=y.cont))
}

concat <- function(x0, x1, y0, y1) {
    x.concat <- c(x0, x0[length(x0)] + x1)
    y.concat <- c(y0, y0[length(x0)] + y1)
    return(list(x=x.concat, y=y.concat))
}

scale <- function(x, y, scale=1) {
    r <- runif(n=1, min=0.1, max=2)
    x.scale <- c(x, x[length(x)] + r*x)
    y.scale <- c(y, y[length(x)] + r*y)
    return(list(x=x.scale, y=y.scale))
}


Join <- function(R1, R2) {
    R <- concat(x0=R1$x, y0=R1$y, x1=R2$x, y1=R2$y)
    return(R)
}

getBallisticSegment <- function(base, ntimes=1) {
    R <- list(x=0, y=0)
    for (i in seq(ntimes)) {
	r <- runif(1)
	if (r < 0.8)
	    R <- Join(R, continue(base$x, -base$y))
	else 
	    R <- Join(R, continue(base$x, base$y))
    }
    return(R)
}

getWaitingSegment <- function(base, ntimes=1) {
    R <- list(x=0, y=0)
    for (j in seq(ntimes)) {
	r <- runif(1)
	if (r < .5)
	    R <- Join(R, flip(base$x, base$y))
	else
	    R <- Join(R, continue(base$x, base$y))
    }
    return(R)
}


N <- 10
runint <- function(n=1, a=1, b=7) { ceiling(runif(n=n, min=a, max=b)) }

base <- getSyntheticPos(base.occ = 10, 
		     speed=c(0.5, -0.5))

old.base <- base

base <- R
base <- old.base
R <- list(x=0, y=0)
{
    for (i in seq(1)) {
        r <- runif(1)
        rint <- plint(1, 10, 1)
        print(paste0(r, " ", rint))
        if (r < 0.85) {
	   R <- Join(R, getBallisticSegment(base, rint))
        } else {
	   R <- Join(R, getWaitingSegment(base, rint))
        }
    }
    
    spar(cx=2.5, col='white')
    matplot(R$x, R$y, type='l', 
    	xlab=tex("t_k = \\sum_k \\tau_k"),
    	ylab=tex("x_k = x_{k-1} \\pm m \\tau_k"))
}
s <- seq(100)
matplot(R$x[s], R$y[s], type='l', 
	xlab=tex("t_k = \\sum_k \\tau_k"),
	ylab=tex("x_k = x_{k-1} \\pm m \\tau_k"))



### Plot of synthetic data
{
img.path = "~/Dropbox/PHD/Papers/JumpStatistics/Images/"
pdf(paste0(img.path, "synthetic_pos.pdf"), width=18, height=16)
spar(cx=3, col='white', mfrow=c(2,2))
par(mar = c(8, 8, 1.3, 1.3), oma = c(2, 2, 2, 2), mgp=c(6,2,0))

s <- seq(1, length(R$x), by=1000)
matplot(R$x, R$y, type='l', 
	xlab=tex("t_k"),
	ylab=tex("x_k"))

s <- seq(1, min(which(R$x > 65000)), by=200)
matplot(R$x[s], R$y[s], type='l', 
	xlab=tex("t_k"),
	ylab=tex("x_k"))

s <- seq(1, min(which(R$x > 6780)), by=10)
matplot(R$x[s], R$y[s], type='l', 
	xlab=tex("t_k"),
	ylab=tex("x_k"))

s <- seq(min(which(R$x > 680)))
matplot(R$x[s], R$y[s], type='l', 
	xlab=tex("t_k"),
	ylab=tex("x_k"))

dev.off()
}



# ---- Jump speed  ----

jd <- ReadDataList(FileList = lf("jump"), idx=3)
jt <- ReadDataList(FileList = lf("jump"), idx=2)
L <- getTSD_hist(dx = jd / jt, n=200, plot=FALSE)

hist(velocities, freq=FALSE, breaks=200)
plot(L$x[,2], cumsum(L$y[,2])/sum(L$y[,2]))

np <- 1000
velocities <- approx(cumsum(L$y[,2])/sum(L$y[,2]), L$x[,2], runif(np))$y
X <- T <- vector(mode='list', length=np)
nsteps <- 10000
for (j in seq(length(velocities))) {
    
    t.k <- x.k <- c(0)
    for (i in seq(2, nsteps)) { 
	tau <- power.law(alpha = 1.7, x0=1, n = 1)
	t.k[i] <- t.k[i-1] + tau

	r <- runif(1)
	v <- ifelse(r < 1/2, velocities[j], -velocities[j])
	x.k[i] <- x.k[i-1] + v * tau
    }
    X[[j]] <- x.k
}
X <- list.cbind.na(X)
plot(t.k, x.k, type='l')

plot(rowMeans(X^2, na.rm=TRUE), type='l')

msd <- rowMeans(X^2, na.rm=TRUE)
mu <- getPowerScaling(t = 1:length(msd), M = msd, z = 0.8)

plot(1:length(mu), mu, type='l')

lines(x, S, type='s')


