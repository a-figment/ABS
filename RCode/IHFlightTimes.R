#' Looking at the free flight distribution for infinite horizon channels
setwd("~/Dropbox/ABS/Studies/S-HPO/Data/IHCompare")


Xc <- tpose(fread("p0-a5-r38-t1--Xc.dat"))
FP <- tpose(fread("p0-a5-r38-t1--FP.dat"))

#' ========= INTERCOLLISIONAL DISTANCE DISTRIBUTION =========

FP <- as.vector(FP[,])
FP <- FP[-which(FP < 0.1)]
L <- getTSD_hist(dx = data.matrix(FP[-which(FP < 0.1)]), n = 100, prune = FALSE, domain=exp.domain)
LF <- LinFit(mat = L$y, time = L$x, range = 80:98,)#show.range=TRUE)

spar(cx=1.5)
matplot(L$x, L$y, type='s', log='xy', add=TRUE,
	col='blue', lwd=2,
	xlab=tex("\\tau"),
	ylab=tex("\\Psi(\\tau)"))
lines(L$x, exp(LF$b) * L$x^(LF$ax), lty='dashed', lwd=2)
#' ==========================================================


#' ========= JUMP TIME ALGORITHM FOR COLLISION TIME =========
t <- apply(FP,2,cumsum)

w <- c(10) #, 250)
eps <- c(5)
atau <- list()
for (j in seq(length(w))) {
    tau <- c()
    for (i in seq(length(lj))) {
        mm <- getMM(x = t[,lj[i]], y = Xc[,lj[i]], 
		      w = w[j], eps=eps[1])

	temp <- mm$tau
	if (length(mm$wait.id) > 0)
	    temp <- temp[-mm$wait.id]
        tau <- c(tau, temp)
    }
    atau[[j]] <- tau
}


#' Try sample IH trajectories with long ballistic jumps
lj <- apply(FP, 2, function(x) ifelse(max(x) > 500, TRUE, FALSE))
lj <- seq(1, ncol(FP))[lj]
lj <- seq(1, ncol(FP))

#' Check jump times appropriately computed
j <- lj[rindex(length(lj), 0)]
matplot(t[,j], Xc[,j], type='o', pch=19,)
mm <- getMM(x = t[,j], y = Xc[,j], w = 10)
ShowJT(mm$df, text=TRUE)

#' ========= JUMP TIME ALGORITHM FOR CONTINUOUS TIME =========
t <- seq(1, 1e5+1, length=1001)
Xt <- tpose(fread("p0-a5-r39-t1--Xt.dat"))
matplot(t, Xt[,j], type='o', pch=19)
mm <- getMM(x = t, y = Xt[,j], w = 1)
ShowJT(mm$df, text=FALSE)

w <- 1
tau <- c()
for (i in seq(ncol(Xt))) {
    temp <- getMM(x = t, y = Xt[,i], w = w)$tau
    tau <- c(tau, temp)
}
#' ===========================================================

# ================= DISTRIBUTION OF TAU ======================
pdf("ihgoodcompare.pdf", width=12, height=12)
tau <- list.cbind.na(atau)
tau <- unique(round(tau, digits=3))
L <- getTSD_hist(dx = data.matrix(tau), n = 100, prune = FALSE, domain=exp.domain)
LF <- LinFit(mat = L$y, time = L$x, range = 44:95, show.range=TRUE)

spar(cx=2)
matplot(L$x, L$y, type='s', log='xy', xlim=c(500,10000), ylim=c(1e-9, 1e-3),
	col=rainbow(10), lwd=3, lty='solid',
	xlab=tex("\\tau"),
	ylab=tex("\\Psi(\\tau)"))
lines(L$x[,1], exp(LF$b[1]) * L$x[,1]^(LF$ax[1]), lty='dashed', lwd=2)

legend('topright', 
       legend=c(tex("\\tau^{-3.00}"),
		"Intercollision distance", 
		"Coarse graining: w = 10", 
		"Coarse graining: w = 250"),
       lty=c('dashed', 'solid', 'solid', 'solid'), 
       col=c('black', 'blue', 'red', 'orange'),
       cex=2, lwd=2)


#' ================ W = 1, 2 comparison ================
FP <- tpose(fread("p0-a5-r39-t1--FP.dat"))
FP <- as.vector(FP[,])
FP <- FP[-which(FP < 0.1)]
w1 <- readRDS("w1tau.RData")
w2 <- readRDS("w2tau.RData")

pdf("IHCompare.pdf", width=12, height=12)
L <- getTSD_hist(dx = list.cbind.na(list(FP, w1, w2)), n = 100, prune = FALSE, domain=exp.domain)

L <- getTSD_hist(dx = data.matrix(tau), n = 50, prune = FALSE, domain=exp.domain)
spar(cx=2)
matplot(L$x, L$y, type='s', log='xy', xlim=c(5,1e4), 
	col=rainbow(10), lwd=3, lty='solid',
	xlab=tex("\\tau"),
	ylab=tex("\\Psi(\\tau)"))

LF <- LinFit(mat = L$y[,1], time = L$x[,1], range = 33:65,)# show.range=TRUE)

lines(L$x[,1], exp(LF$b[1]) * L$x[,1]^(LF$ax[1]), lty='dashed', lwd=2)
lines(L$x[,3], exp(LF$b[1]) * L$x[,3]^(LF$ax[1]), lty='dashed', lwd=2)

legend('topright', 
       legend=c(tex("\\tau^{-3.01}"),
		"Intercollision distance", 
		"Coarse graining: w = 1", 
		"Coarse graining: w = 2"),
       lty=c('dashed', 'solid', 'solid', 'solid'), 
       col=c('black', rainbow(10)[1:3]),
       cex=2, lwd=2)

legend('topright', 
       legend=c(tex("\\tau^{-3.01}"), 
		"Intercollision distance: w = 1", 
		"Coarse graining: w = 1", "Coarse graining: w = 2"),
       lty=c('dashed', 'solid', 'solid', 'solid'), 
       col=c('black', rainbow(10)[1:3]),
       cex=2, lwd=2)


#' DISTRIBUTION OF WINDOW TIMES
m <- apply(FP, 2, mean)
rollm <- rollmean(x=FP[,8], k=5)

hist(rollm, breaks=10)





