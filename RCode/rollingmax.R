

# ===== UNDERSTANDING THE ROLLING MAX PROC. =====
setwd("../Studies/S-HPO/Data/")
setwd("../Studies/S-JumpStatistics/Data/")
X <- data.matrix(fread("p0-a5-r10-t1--Xt.dat"))
t <- seq(1, 1e5+1, length=1001)

s <- seq(1000)
X <- X[101,s]

# Animate windows
n <- length(s)
w <- 2
I <- c(1:w, n+1-1:w)
y.max <- rollmax(zoo(X), k=2*w+1, align="center")
y.min <- -rollmax(zoo(-X), k=2*w+1, align="center")

# MAIN
{
    xc.min <- xc.max <- yc.max <- yc.min <- c()
    for (i in seq(length(s))) {
        #png(paste0(i, "-frame.png"), width=1440, height=1024)
        spar(cx=2)
        matplot(t[s], X, type='o', pch=19, lwd=2, ylim=c(-50, 450), 
    	    xlab=tex("t_k"), 
    	    ylab=tex("x(t_k)"),)
        #rect(xleft=t[i], xright=t[i + 2*w+1], ytop=1e6, ybot=-1e6, col=rgb(0,0,1,alpha=0.2))
        S <- seq(i, i+2*w+1)
        k <- which.max(X[S])
        K <- which.min(X[S])
        xc.min <- c(xc.min, t[S[K]])
        xc.max <- c(xc.max, t[S[k]])
        yc.min <- c(yc.min, min(X[S]))
        yc.max <- c(yc.max, max(X[S]))

        points(t[S[k]], yc.max[i], pch=19, col='red', cex=1.2)
        arrows(x0=t[S[k]]-800, x1=t[S[k]]-200, y0=yc.max[i], y1=yc.max[i], angle=45, lwd=2.6, length = 0.1, col='red')
        text(x = t[S[k]]-3900, y=yc.max[i], labels = "WINDOW MAX.", col = 'red', cex = 1.5)
        points(xc.max, yc.max, col='red', cex=1.5, pch=19)
        points(t[S[K]], yc.min[i], pch=0, col='blue', cex=1.75, lwd=2.5)
        arrows(x0=t[S[K]]+800, x1=t[S[K]]+200, y0=yc.min[i], y1=yc.min[i], angle=45, lwd=2.6, length = 0.1, col='blue')
        text(x = t[S[K]]+3900, y=y.min[i], labels = "WINDOW MIN.", col = 'blue', cex = 1.5)
        points(xc.min, yc.min, col='blue', cex=1.75, pch=0, lwd=2.5)
        #if (i != 50)
        #	dev.off()
    }
}


w <- 100
legend(x=25200, y=460, 
       legend=c("Window maximum", "Window minimum", "Turning pt. maximum", "Turning pt. minimum"), 
       pch=c(19, 0, 17, 18), box.lwd=0, box.col='gray',
       col=c('red', 'blue', 'red', 'blue'), cex=2)
rect(xleft=t[24], xright=t[24 + 2*w+1], ytop=300, ybot=-1e6, col=rgb(0,0,1,alpha=0.1))#, border='gray')


# Turning point maximum
nf <- c(1, seq(3, nrow(MM$df), by=2))
points(MM$df$xl[nf], MM$df$yl[nf], col='blue', pch=18, cex=3.5)

# Turning point minimum
nf <- seq(1, nrow(MM$df), by=2)
points(MM$df$xr[nf], MM$df$yr[nf], col='red', pch=17, cex=3.5)


# INSETS
{
    par(fig=c(0.65, 0.96, 0.08, 0.5), new=TRUE)
    delta <- y.max - X[-I]
    matplot(xc.max[-I], delta, col='magenta', pch=19, type='o', lwd=2,
            xlab='', ylab=TeX(paste0("(WIN. MAX.) - ", "x(t_k)")))
    abline(v=0, h=0)
    z1 <- which(delta <= 0)
    points(xc.max[-I][z1], rep(0, length(z1)), col='red', pch=17, cex=1.7)

    par(fig=c(0.06, 0.35, 0.5, 0.92), new=TRUE)
    delta <- y.min - X[-I]
    matplot(xc.min[-I], y.min - X[-I], col='magenta', pch=19, type='o', lwd=2, 
            xlab='', ylab=TeX(paste0("(WIN. MIN.) - ", "x(t_k)")))
    abline(v=0, h=0)
    z2 <- which(abs(delta) <= 0)
    points(xc.min[-I][z2], rep(0, length(z2)), col='blue', pch=18, cex=2)
}

ShowJT <- function(df, text=TRUE) {

    N <- sum(abs(df$yr - df$yl))
    for (k in seq(nrow(df))) {
	
	yloc <- min(df$yl[k], df$yr[k])
   	arrows(x0 = df$xl[k], x1 = df$xr[k], 
	       y0 = yloc, y1 = yloc, 
	       angle = 45, code=3, col='magenta', length=0.1, lwd=2)

	lines(c(df$xr[k], df$xr[k]), 
	      c(yloc, df$yr[k]), col='black', lty='dashed', lwd=2.5)
	lines(c(df$xl[k], df$xl[k]), 
	      c(df$yl[k], yloc), col='black', lty='dashed', lwd=2.5)

	if (text) {
	    text(labels = TeX(paste0("$\\tau_{",k,"}$ = ", 
			      format(df$xr[k]-df$xl[k], digits=1))), 
		x = (df$xl[k] + df$xr[k])/2, 
		y=yloc - N*0.02, lwd=2, cex=1.75, col='purple')
	}
    }
}


#' @brief Applies a rolling window to a curve (x,y) 
#' @w Window half width
#' @span \in (0,1) is the 'degree of smoothing'. Small @span = tighter fit
#'
argmax <- function(x, y, w, span) {
    y.max <- zoo::rollmax(zoo::zoo(y), k=2*w+1, align="center")
    delta <- y.max - y[-c(1:w, length(y)+1-1:w)]
    idx <- which(delta <= 0) + w
    return(idx)
}

#' @brief Uses @argmax() to find the minima/maxima of a curve
#' @w Window half width
#' @span \in (0,1) is the 'degree of smoothing'. Small @span = tighter fit
#'
getMM <- function(x, y, w, span) {
    Min <- argmax(x=x, -y, w=w, span=span)
    Max <- argmax(x=x, y , w=w, span=span)
   
    if (length(Max) != length(Min)) {
	lb <- min(length(Max), length(Min))
	Max <- Max[1:lb]
	Min <- Min[1:lb]
    }

    tau.mins = x[Min]
    tau.maxs = x[Max]

    taus <- c(tau.maxs, tau.mins)
    SMM <- sort(taus)	       # MinMax in ascending order
    IDMM <- order(taus)	       # Indices which sort the min/max times

    IDs <- c(Max, Min)[IDMM] 
    Y <- y[IDs]


    # Consider alternating jumps
    df <- data.frame(xl=SMM[1:(length(SMM)-1)], xr=SMM[2:length(SMM)], 
		     yl=Y[1:(length(SMM)-1)], yr=Y[2:length(SMM)])

    # Consider pairs of times only
    #pos <- seq(3, length(SMM), by=2)
    #SMM <- insert(vect=SMM, pos=pos, elems=NA)

    #idx <- seq(1, length(SMM), by=2)
    #df <- data.frame(xl=SMM[idx], xr=SMM[idx], yl=Y[idx], yr=Y[idx])

    # Jump times - Difference between min/max pairs
    tau <- diff(SMM)
    tau <- tau[-which(is.na(tau))]

    # Central object for this study
    TP  <- list(x=SMM, 
		y=y,
		df=df,
		tau=tau,
		i=IDs,
		y=y)
    return(TP)
}

x <- seq(1,50000,by=10)
f <- splinefun(x=t[s], y = as.vector(X))(x)
mm <- getMM(t[s], X, w=4, span=1)
mm <- getMM(x, f, w=100, span=1)

spar(cx=2)
matplot(t[s], X, type='o', pch=19, lwd=2,
matplot(x, f, type='l', lwd=2,
    xlab=tex("t_k = 1, 1001, ..., 10^7 + 1"), 
    ylab=tex("x(t_k)"),)
ShowJT(mm$df)


#' Differently prepared time series x(t_k)
setwd("../Studies/S-HPO/Data/")
X  <- tpose(fread("p0-a5-r4-t1--Xt.dat"))
X  <- tpose(fread("p0-a5-r4-t1--Xc.dat"))
Yc <- tpose(fread("p0-a5-r4-t1--Yc.dat"))
FP <- tpose(fread("p0-a5-r4-t1--FP.dat"))

j <- 2384 #rindex(ncol(X), 0)
w <- 12 # coltime
w <- 9  # contime

t <- cumsum(FP[,j])
X <- X[,j]

t <- seq(1, 1e4+1, length=1001)[1:316]
X <- X[1:316,j]

pdf("test.pdf", width=12, height=16)
spar(mfrow=c(3,1), cx=2)
spar(cx=2)
matplot(t, X, type='l', pch=19, lwd=2.5, #xlim=c(0,2000), ylim=c(-100,0),
    	xlab=tex("t_k"), #add=TRUE, col='red',
    	ylab=tex("x(t_k)"),)

# Comparison with raw data
MM <- getMM(t, X, w=w, pairs=FALSE)
ShowJT(MM$df, text=FALSE)

text(x=max(t), y=max(X)-100, labels = "1", cex=5)
text(x=max(t), y=max(X)-100, labels = "2", cex=5)
text(x=max(t), y=max(X)-100, labels = "3", cex=5)

# Comparison with smoothing
y.smooth <- loess(X ~ t, span=0.2)$fitted  # Smoothing function 
lines(t, y.smooth, col='red', lwd=1.5)
MM <- getMM(t, y.smooth, w=w)
ShowJT(MM$df, text=FALSE)

# Comparison with interpolating
x <- seq(1,t[length(t)],length=5*length(t))
f <- approxfun(x=t, y = as.vector(X[,j]))(x)
lines(x, f, col='red', lwd=2)
MM <- getMM(x, f, w=w*10)
ShowJT(MM$df, text=FALSE)



