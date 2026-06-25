
setwd("../Studies/S-HPO/Data/")

t <- seq(1,1e7+1,length=10001)

files <- paste0("p0-a", seq(15,17), "-r1--X.dat")  # IRRATIONALS pi/3, pi/2, 2pi/3
files <- paste0("p0-a", seq(4,6), "-r1--X.dat")    # RATIONALS pi/3, pi/2, 2pi/3
pref <- "JT-r1-"

#' @brief 
#' @jdx Jump displacements computed up to tim
JumpDispMoments <- function(t, jdx, p) {

    mom <- list()
    for (i in seq(length(p))) {
	mom[[i]] <- apply(jdx, 2, function(x) mean(abs(x)^p[i], na.rm=TRUE))
    }

    return(list.cbind.na(mom))
}

getJT <- function(t, X, w) {
    TAU <- list()
    for (i in seq(length(w))) {
	MM <- getMM(x = t, y = X, w = w[i])
	TAU[[i]] <- MM$tau
    }
    return(TAU)
}

getJD <- function(t, X, w) {
    DY <- list()
    for (i in seq(length(w))) {
	MM <- getMM(x = t, y = X, w = w[i])
	DY[[i]] <- MM$df$yr - MM$df$yl
    }
    return(DY)
}

S <- seq(1000, 10000, by=500)
S <- c(100, 1000, 10000)
S <- 10001
eps <- 1

w <- c(10, 40, 160)

J <- list()
for (f in seq(length(files))) {
    X <- tpose(fread(files[f]))

    jump <- c()
    jumps <- vector(mode='list', length=length(w))
    for (i in seq(1000)) {
	tau <- getJT(t, X[,i], w=w)
	#dy  <- getJD(t, X[,i], w=w)

	if (length(tau) > 0) {
	    for (j in seq(length(tau))) {
	        jumps[[j]] <- c(jumps[[j]], tau[[j]]) 
	    }
	}
    }
    J[[f]] <- jump
}


matplot((dy)[,2] / tau[,2], pch=19, lty='solid', log='')

tau <- list.cbind.na(tau)
dy <- list.cbind.na(dy)

L <- getTSD_hist(dx = abs(dy) / tau, n = 10, prune = FALSE)

L <- getTSD_hist(dx = list.cbind.na(jumps), n = 200, prune = FALSE, domain=exp.domain)

matplot(L$x, L$y, type='s', lty='solid', lwd=1, log='xy')




source("~/Dropbox/Rlibs.R")
r <- i <- rindex(ncol(X), 0)
{    
    s <- seq(10000)
    eps <- 1
    MM <- getMM(x = t, y = X[,r], w = w[3], eps=eps)
    #MM$df <- ifelse(length(MM$wait.id) == 0, MM$df, MM$df[-MM$wait.id,])
    spar(cx=2)
    matplot(t[s], X[s,r], type='l', pch=19, lwd=1.5, #ylim=c(-5000,5000),
        xlab=tex("t"), 
        ylab=tex("x(t)"),)
    ShowJT(MM$df, text=FALSE, eps=eps)

    par(fig=c(0.075, 0.42, 0.075, 0.6), new=TRUE)
    s <- seq(2000,3500)
    matplot(t[s], X[s,r], type='l', pch=19, lwd=1.5, #ylim=c(-5000,5000),
        xlab="", yaxt='n',
        ylab="")
    ShowJT(MM$df, text=FALSE, eps=eps)
}


#' PLOT 0 - ILLUSTRATION OF LONG JUMPS IN TIME SERIES
files <- c(paste0("p0-a", seq(15,17), "-r1--X.dat"), paste0("p0-a", seq(4,6), "-r1--X.dat"))
r <- rindex(100, 0)
X <- cbind(tpose(fread(files[3]))[,584], tpose(fread(files[6]))[,125])
X <- cbind(tpose(fread(files[3]))[,r], tpose(fread(files[6]))[,r])

pdf("longjump-ex.pdf", width=24, height=12)
spar(mfrow=c(1,2), cx=2)
matplot(t, X[,1], type='l', pch=19, lwd=2, 
	xlab=tex("t"), 
        ylab=tex("x(t)"),)
matplot(t, X[,2], type='l', pch=19, lwd=2, 
	xlab=tex("t"), 
        ylab=tex("x(t)"),)

#' PLOT 1 - ILLUSTRATING DEP. OF WINDOW WIDTH
X <- tpose(fread(files[3]))
r <- 2268
eps <- 1
w <- c(10, 40, 160, 200)

spar(mfrow=c(2,2), cx=2)
spar(mfrow=c(1,2), cx=2)
r <- rindex(ncol(X), 0)
for (i in seq(length(w),1)) {
    ifelse(i == 1, s <- seq(1000, 2000), s <- seq(nrow(X)))
    MM <- getMM(x = t[s], y = X[s,r], w = w[i], eps=eps)
    matplot(t[s], X[s,r], type='l', pch=19, lwd=1.5, 
        xlab=tex("t"), 
        ylab=tex("x(t)"),)
    ShowJT(MM$df, text=FALSE, eps=eps, len=0.05, wid=1.5)
    dy <- max(X[s,r]) - min(X[s,r])
    text(x=t[s[length(s)]], y=max(X[s,r])-0.025*dy, labels = seq(4)[length(w)-i+1], cex=4, col='blue')
}

#' PLOT 2 - WIDTH DEP. FOR MANY RANDOM TRAJECTORIES
X <- tpose(fread(files[3]))
pdf(file="RSample--i2PI3.pdf", width=15, height=15)
eps <- 1
w <- c(10, 40, 160, 200)
for (p in seq(10)) {
    spar(mfrow=c(2,2), cx=2)
    r <- rindex(ncol(X), 0)
    for (i in seq(length(w),1)) {
        ifelse(i == 1, s <- seq(1000, 2000), s <- seq(nrow(X)))
        MM <- getMM(x = t[s], y = X[s,r], w = w[i], eps=eps)
        matplot(t[s], X[s,r], type='l', pch=19, lwd=1.5, 
            xlab=tex("t"), 
	    main=paste0("w = ", w[i]),
            ylab=tex("x(t)"),)
        ShowJT(MM$df, text=FALSE, eps=eps, len=0.05, wid=1.5)
	dy <- max(X[s,r]) - min(X[s,r])
        text(x=t[s[length(s)]], y=max(X[s,r])-0.025*dy, labels = seq(4)[length(w)-i+1], cex=4, col='blue')
    }
}
dev.off()


J <- jump.ti
J <- list.cbind.na(list(readRDS("JDX_pi3_w10_e1000.RData")[,1], 
       readRDS("JDX_pi2_w10_e1000.RData")[,1], 
       readRDS("JDX_2pi3_w10_e1000.RData")[,1]))
J <- list.cbind.na(J)
L <- getTSD_hist(dx = abs(J[,seq(1,19,by=2)]), n = 100, domain = exp.domain, prune = FALSE)
L <- getTSD_hist(dx = J, n = 500, prune = FALSE)
L <- getTSD_hist(dx = data.matrix(J), n = 200, domain = exp.domain, prune = FALSE)

LF <- LinFit(mat = L$y[,1], time = L$x[,1], range = 50:85+90, logs = TRUE, show.range=TRUE)

spar(cx=2)
matplot(L$x, L$y, type='s', log='xy', #xlim=c(-5e5, 5e5),
	xlab=tex("\\delta y"),
	ylab=tex("\\psi(\\delta y, s)"),
	col=rainbow(10), lwd=2, lty='solid')
lines(L$x[,1], exp(LF$b[1]) * L$x[,1]^LF$ax[1])

legend('topright', legend=paste0("s = ", S[seq(1,19,by=2)]*10^3), col=rainbow(10), lty='solid', lwd=2, cex=2)
legend('topright', legend=c(tex("\\alpha = \\pi/3"), 
			    tex("\\alpha = \\pi/2"), 
			    tex("\\alpha = 2\\pi/3")),
			    col=rainbow(10), lty='solid', lwd=2, cex=2)



par(fig=c(0.075, 0.5, 0.075, 0.5), new=TRUE)
plot(P, spec$ax, pch=19, col='blue', cex=1.5, xlab="p", ylab=tex("M_p"))
lines(P[], P[] - 1.67, col='black', lty='dashed', lwd=2)

